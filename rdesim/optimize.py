"""Geometry optimization driver.

Each evaluation writes a case and runs ``optimization.solverCommand`` in that
case directory. The command is expected to leave
``postProcessing/rdeObjectives/objectives.dat`` with one ``name value`` pair
per line. ``rdeFluid`` writes that file. A short run can use any executable
that writes the same file.

Evaluations are kept. Every evaluation lives in ``eval-<fingerprint>`` under
``optimization.workDirectory``, where the fingerprint covers the whole case
(design values, geometry, fuel, reaction model, numerics). A finished
evaluation stores its objective quantities in ``evaluation.json`` and is reused
by later runs without starting the solver. An interrupted evaluation has no
such file, so the next run starts the solver in the same folder; with
``startFrom latestTime`` OpenFOAM continues from the last time it wrote.
``optimization.maxNewEvaluations`` or ``--budget`` limits how many solver runs
one call may start. Objective and constraint expressions are evaluated from the
stored quantities, so editing them does not re-run anything.
"""

from __future__ import annotations

import copy
import itertools
import json
import shutil
import subprocess
from pathlib import Path

from .casewriter import prepare_case
from .config import ConfigError, as_bool, get_path, load_config, set_path
from .expr import ExpressionError, evaluate_with_names
from .runs import fingerprint, latest_time


class _Budget:
    def __init__(self, limit):
        self.limit = limit
        self.used = 0
        self.exhausted = False

    def available(self) -> bool:
        return self.limit is None or self.used < self.limit

    def spend(self) -> None:
        self.used += 1


def optimize(
    config_path: str | Path,
    budget: int | None = None,
    fresh: bool = False,
    model: str | None = None,
) -> dict:
    config_path = Path(config_path).resolve()
    base = load_config(config_path, model)
    opt = base["raw"]["optimization"]
    work = Path(str(opt["workDirectory"]))
    if not work.is_absolute():
        work = config_path.parent / work
    if "resume" in opt and not as_bool(opt["resume"]):
        fresh = True
    if fresh and work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    variables = list(opt["designVariables"])
    samples = [_samples(variable) for variable in variables]
    limit = int(float(opt["maxEvaluations"]))
    timeout = float(opt["timeoutSeconds"])
    if budget is None and "maxNewEvaluations" in opt:
        budget = int(float(opt["maxNewEvaluations"]))
    spend = _Budget(budget)
    records = []
    best = None
    deferred = []

    def run(values):
        nonlocal best
        if len(records) >= limit:
            return None
        record = _evaluate(base, variables, values, work, opt, timeout, spend)
        if record["status"] == "deferred":
            deferred.append(record)
            spend.exhausted = True
            return None
        records.append(record)
        best = _better(best, record, opt["objective"]["sense"])
        _summary(work, records, best, opt, spend, deferred)
        return record

    if opt["method"] == "grid":
        for values in itertools.product(*samples):
            if run(values) is None and (len(records) >= limit or spend.exhausted):
                break
    elif opt["method"] == "coordinate":
        levels = int(float(opt["coordinateLevels"]))
        if levels < 1:
            raise ConfigError("'optimization.coordinateLevels' must be at least 1")
        current = [column[0] for column in samples]
        run(tuple(current))
        for _level in range(levels):
            for index, column in enumerate(samples):
                level_best = None
                level_value = None
                for value in column:
                    if spend.exhausted:
                        break
                    trial = list(current)
                    trial[index] = value
                    record = run(tuple(trial))
                    if record is None:
                        break
                    if record["feasible"] and (
                        level_best is None
                        or (
                            opt["objective"]["sense"] == "maximize"
                            and record["objective"] > level_best
                        )
                        or (
                            opt["objective"]["sense"] == "minimize"
                            and record["objective"] < level_best
                        )
                    ):
                        level_best = record["objective"]
                        level_value = value
                if level_value is not None:
                    current[index] = level_value
                if len(records) >= limit or spend.exhausted:
                    break
            if spend.exhausted:
                break
    else:
        raise ConfigError(f"Unknown optimization.method '{opt['method']}'")
    summary = _summary(work, records, best, opt, spend, deferred)
    summary["budgetExhausted"] = spend.exhausted
    summary["newEvaluations"] = spend.used
    return summary


def _samples(variable) -> list[float]:
    count = int(float(variable["samples"]))
    low = float(variable["min"])
    high = float(variable["max"])
    if high < low:
        raise ConfigError(
            f"Design variable '{variable['name']}' has max below min"
        )
    if count == 1:
        return [low]
    return [low + (high - low) * i / (count - 1) for i in range(count)]


def _evaluate(base, variables, values, work: Path, opt, timeout: float, spend: _Budget) -> dict:
    raw = _clone(base["raw"])
    names = []
    for variable, value in zip(variables, values):
        set_path(raw, variable["name"], value)
        names.append(f"{variable['name']}={value:.8g}")
    record = {
        "directory": "",
        "values": names,
        "status": "pending",
        "objective": None,
        "feasible": False,
        "message": "",
        "quantities": {},
        "reused": False,
    }
    violation = _design_constraints(raw, opt["constraints"])
    if violation:
        record["status"] = "infeasible"
        record["message"] = violation
        return record
    config = {
        "path": base["path"],
        "raw": raw,
        "fuel": base["fuel"],
        "fuelPath": base["fuelPath"],
    }
    case_dir = work / f"eval-{fingerprint(config)}"
    record["directory"] = str(case_dir)
    stored = _read_record(case_dir)
    if stored is not None and stored.get("solved"):
        record["reused"] = True
        record["quantities"] = {k: float(v) for k, v in stored["quantities"].items()}
        return _score(record, opt)
    # Every evaluation folder resumes from its latest time after an interruption.
    raw["output"].setdefault("startFrom", "latestTime")
    try:
        prepared = prepare_case(config, case_dir)
        if not prepared["hyperbolicity"].hyperbolic:
            record["status"] = "infeasible"
            record["message"] = prepared["hyperbolicity"].text()
            return record
    except (ConfigError, OSError) as exc:
        record["status"] = "error"
        record["message"] = str(exc)
        return record
    if not spend.available():
        record["status"] = "deferred"
        record["message"] = "evaluation budget used up; run again to continue"
        return record
    spend.spend()
    resumed = latest_time(case_dir)
    if resumed is not None:
        record["message"] = f"resumed from time {resumed}. "
    command = str(opt["solverCommand"]).split()
    try:
        completed = subprocess.run(
            command,
            cwd=case_dir,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        record["status"] = "error"
        record["message"] += str(exc)
        return record
    if completed.returncode != 0:
        record["status"] = "error"
        tail = (completed.stderr or completed.stdout or "").strip().splitlines()
        record["message"] += "solver exit " + str(completed.returncode)
        if tail:
            record["message"] += ": " + tail[-1]
        return record
    path = case_dir / "postProcessing" / "rdeObjectives" / "objectives.dat"
    if not path.is_file():
        record["status"] = "error"
        record["message"] += f"Solver did not write {path}"
        return record
    quantities = _read_objectives(path)
    record["quantities"] = quantities
    _write_record(case_dir, {"solved": True, "values": names, "quantities": quantities})
    return _score(record, opt)


def _score(record, opt) -> dict:
    quantities = record["quantities"]
    try:
        value = evaluate_with_names(opt["objective"]["expression"], quantities)
    except ExpressionError as exc:
        record["status"] = "error"
        record["message"] += str(exc)
        return record
    violation = _result_constraints(quantities, opt["constraints"])
    record["objective"] = value
    if violation:
        record["status"] = "infeasible"
        record["message"] += violation
        return record
    record["status"] = "ok"
    record["feasible"] = True
    return record


def _read_record(case_dir: Path):
    path = case_dir / "evaluation.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _write_record(case_dir: Path, data: dict) -> None:
    (case_dir / "evaluation.json").write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _design_constraints(raw, constraints) -> str:
    for constraint in constraints:
        quantity = str(constraint["quantity"])
        if not _is_config_path(quantity):
            continue
        try:
            value = float(get_path(raw, quantity))
        except (ConfigError, TypeError, ValueError):
            return f"Cannot read design constraint '{quantity}'"
        message = _bound(quantity, value, constraint)
        if message:
            return message
    return ""


def _result_constraints(quantities, constraints) -> str:
    for constraint in constraints:
        quantity = str(constraint["quantity"])
        if _is_config_path(quantity):
            continue
        if quantity not in quantities:
            return f"Constraint quantity '{quantity}' was not written by the solver"
        message = _bound(quantity, quantities[quantity], constraint)
        if message:
            return message
    return ""


def _bound(name, value, constraint) -> str:
    if "min" in constraint and value < float(constraint["min"]):
        return f"{name}={value:.8g} is below min {constraint['min']}"
    if "max" in constraint and value > float(constraint["max"]):
        return f"{name}={value:.8g} is above max {constraint['max']}"
    return ""


def _is_config_path(quantity: str) -> bool:
    return "." in quantity or "[" in quantity


def _better(best, record, sense: str):
    if not record["feasible"]:
        return best
    if best is None:
        return record
    if sense == "maximize":
        return record if record["objective"] > best["objective"] else best
    return record if record["objective"] < best["objective"] else best


def _read_objectives(path: Path) -> dict[str, float]:
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ConfigError(f"Bad objectives line in {path}: {line}")
        values[parts[0]] = float(parts[1])
    return values


def _clone(raw):
    return copy.deepcopy(raw)


def _summary(work: Path, records, best, opt, spend: _Budget, deferred) -> dict:
    reused = sum(1 for record in records if record["reused"])
    lines = [
        f"method {opt['method']}",
        f"objective {opt['objective']['name']} {opt['objective']['sense']} {opt['objective']['expression']}",
        f"evaluations {len(records)}",
        f"reused {reused}",
        f"newSolverRuns {spend.used}",
    ]
    for index, record in enumerate(records):
        tag = " reused" if record["reused"] else ""
        lines.append(
            f"eval-{index:03d} {record['status']}{tag} {record['objective']} "
            f"{' '.join(record['values'])} {record['message']}".rstrip()
        )
    if deferred:
        lines.append(
            f"budgetExhausted yes. {len(deferred)} evaluation(s) not started; "
            "run again to continue. Finished evaluations are reused."
        )
    if best is None:
        lines.append("best none")
    else:
        lines.append(
            "best "
            + str(best["objective"])
            + " "
            + best["directory"]
            + " "
            + " ".join(best["values"])
        )
    (work / "summary").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"work": work, "records": records, "best": best, "text": "\n".join(lines)}

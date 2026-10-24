"""Geometry optimization driver.

Each evaluation writes a case and runs ``optimization.solverCommand`` in that
case directory. The command is expected to leave
``postProcessing/rdeObjectives/objectives.dat`` with one ``name value`` pair
per line. ``rdeFluid`` writes that file. A short run can use any executable
that writes the same file.
"""

from __future__ import annotations

import copy
import itertools
import shutil
import subprocess
from pathlib import Path

from .casewriter import prepare_case
from .config import ConfigError, get_path, load_config, set_path
from .expr import ExpressionError, evaluate_with_names


def optimize(config_path: str | Path) -> dict:
    config_path = Path(config_path).resolve()
    base = load_config(config_path)
    opt = base["raw"]["optimization"]
    work = Path(str(opt["workDirectory"]))
    if not work.is_absolute():
        work = config_path.parent / work
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    variables = list(opt["designVariables"])
    samples = [_samples(variable) for variable in variables]
    limit = int(float(opt["maxEvaluations"]))
    timeout = float(opt["timeoutSeconds"])
    records = []
    best = None

    def run(values):
        nonlocal best
        if len(records) >= limit:
            return None
        record = _evaluate(
            base,
            variables,
            values,
            work / f"eval-{len(records):03d}",
            opt,
            timeout,
        )
        records.append(record)
        best = _better(best, record, opt["objective"]["sense"])
        return record

    if opt["method"] == "grid":
        for values in itertools.product(*samples):
            if run(values) is None and len(records) >= limit:
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
                if len(records) >= limit:
                    break
    else:
        raise ConfigError(f"Unknown optimization.method '{opt['method']}'")
    summary = _summary(work, records, best, opt)
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


def _evaluate(base, variables, values, case_dir: Path, opt, timeout: float) -> dict:
    raw = _clone(base["raw"])
    names = []
    for variable, value in zip(variables, values):
        set_path(raw, variable["name"], value)
        names.append(f"{variable['name']}={value:.8g}")
    record = {
        "directory": str(case_dir),
        "values": names,
        "status": "pending",
        "objective": None,
        "feasible": False,
        "message": "",
        "quantities": {},
    }
    violation = _design_constraints(raw, opt["constraints"])
    if violation:
        record["status"] = "infeasible"
        record["message"] = violation
        return record
    try:
        config = {
            "path": base["path"],
            "raw": raw,
            "fuel": base["fuel"],
            "fuelPath": base["fuelPath"],
        }
        # Fuel overrides inside raw are already merged into base["fuel"].
        # Design variables do not change the fuel. Re-merge is unnecessary.
        prepared = prepare_case(config, case_dir)
        if not prepared["hyperbolicity"].hyperbolic:
            record["status"] = "infeasible"
            record["message"] = prepared["hyperbolicity"].text()
            return record
    except (ConfigError, OSError) as exc:
        record["status"] = "error"
        record["message"] = str(exc)
        return record
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
        record["message"] = str(exc)
        return record
    if completed.returncode != 0:
        record["status"] = "error"
        tail = (completed.stderr or completed.stdout or "").strip().splitlines()
        record["message"] = "solver exit " + str(completed.returncode)
        if tail:
            record["message"] += ": " + tail[-1]
        return record
    path = case_dir / "postProcessing" / "rdeObjectives" / "objectives.dat"
    if not path.is_file():
        record["status"] = "error"
        record["message"] = f"Solver did not write {path}"
        return record
    quantities = _read_objectives(path)
    try:
        value = evaluate_with_names(opt["objective"]["expression"], quantities)
    except ExpressionError as exc:
        record["status"] = "error"
        record["message"] = str(exc)
        record["quantities"] = quantities
        return record
    violation = _result_constraints(quantities, opt["constraints"])
    record["quantities"] = quantities
    record["objective"] = value
    if violation:
        record["status"] = "infeasible"
        record["message"] = violation
        return record
    record["status"] = "ok"
    record["feasible"] = True
    return record


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


def _summary(work: Path, records, best, opt) -> dict:
    lines = [
        f"method {opt['method']}",
        f"objective {opt['objective']['name']} {opt['objective']['sense']} {opt['objective']['expression']}",
        f"evaluations {len(records)}",
    ]
    for index, record in enumerate(records):
        lines.append(
            f"eval-{index:03d} {record['status']} {record['objective']} {' '.join(record['values'])} {record['message']}"
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

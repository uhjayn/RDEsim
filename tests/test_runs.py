"""Run folders, restart settings, and the resumable optimization driver."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from rdesim.casewriter import prepare_case
from rdesim.config import load_config
from rdesim.optimize import optimize
from rdesim.runs import fingerprint, run_directory, run_name
from tests.test_rdesim import EXAMPLE, ROOT, _raw_text

COUNTING_SOLVER = """\
import sys
from pathlib import Path
sys.path.insert(0, {root!r})
Path({log!r}).open("a").write(Path.cwd().name + "\\n")
fail = Path({flag!r})
if fail.exists():
    sys.exit(3)
runpy_target = {stub!r}
import runpy
runpy.run_path(runpy_target, run_name="__main__")
"""


def _optimizer_config(folder: Path, samples: int, **extra) -> Path:
    config = load_config(EXAMPLE)
    opt = config["raw"]["optimization"]
    opt["method"] = "grid"
    opt["maxEvaluations"] = 10
    opt["designVariables"] = [
        {"name": "geometry.steps[0].height", "min": 0.001, "max": 0.003, "samples": samples}
    ]
    opt["constraints"] = []
    opt["workDirectory"] = "optimization"
    opt["solverCommand"] = f"python3 {folder / 'solver.py'}"
    opt.update(extra)
    config["raw"]["fuel"]["library"] = str(ROOT / "fuels")
    path = folder / "rdeConfig"
    path.write_text(_raw_text(config["raw"]), encoding="utf-8")
    return path


def _write_solver(folder: Path, fail: bool = False) -> Path:
    log = folder / "solver.log"
    flag = folder / "fail"
    if fail:
        flag.write_text("x", encoding="utf-8")
    elif flag.exists():
        flag.unlink()
    (folder / "solver.py").write_text(
        COUNTING_SOLVER.format(
            root=str(ROOT), log=str(log), flag=str(flag), stub=str(ROOT / "tests" / "stub_solver.py")
        ),
        encoding="utf-8",
    )
    return log


def _runs(log: Path) -> list[str]:
    return log.read_text(encoding="utf-8").split() if log.exists() else []


class RunFolderTests(unittest.TestCase):
    def test_same_physics_same_folder_and_endtime_ignored(self):
        first = load_config(EXAMPLE)
        second = load_config(EXAMPLE)
        second["raw"]["numerics"]["endTime"] = 1.0e-3
        second["raw"]["output"]["writeInterval"] = 5.0e-5
        self.assertEqual(fingerprint(first), fingerprint(second))
        self.assertEqual(run_directory(first), run_directory(second))

    def test_geometry_and_model_open_a_new_folder(self):
        first = load_config(EXAMPLE)
        second = load_config(EXAMPLE)
        second["raw"]["geometry"]["steps"][0]["height"] = 0.0021
        self.assertNotEqual(fingerprint(first), fingerprint(second))
        self.assertTrue(run_name(first).startswith("annulus2D-oneStepProgress-"))

    def test_restart_and_checkpoint_entries(self):
        config = load_config(EXAMPLE)
        config["raw"]["output"]["startFrom"] = "latestTime"
        config["raw"]["output"]["checkpoint"] = {
            "enabled": "yes",
            "writeControl": "clockTime",
            "writeInterval": 600,
        }
        with tempfile.TemporaryDirectory() as tmp:
            prepare_case(config, Path(tmp) / "case")
            control = (Path(tmp) / "case" / "system" / "controlDict").read_text(encoding="utf-8")
        self.assertIn("startFrom       latestTime;", control)
        self.assertIn("type            writeObjects;", control)
        self.assertIn("writeControl    clockTime;", control)
        self.assertIn("writeInterval   600;", control)

    def test_default_start_is_unchanged(self):
        config = load_config(EXAMPLE)
        with tempfile.TemporaryDirectory() as tmp:
            prepare_case(config, Path(tmp) / "case")
            control = (Path(tmp) / "case" / "system" / "controlDict").read_text(encoding="utf-8")
        self.assertIn("startFrom       startTime;", control)
        self.assertNotIn("writeObjects", control)


class ResumableOptimizerTests(unittest.TestCase):
    def test_budget_then_resume_reuses_finished_evaluations(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            log = _write_solver(folder)
            path = _optimizer_config(folder, samples=3)
            first = optimize(path, budget=2)
            self.assertEqual(len(_runs(log)), 2)
            self.assertTrue(first["budgetExhausted"])
            self.assertEqual(len(first["records"]), 2)
            self.assertIn("budgetExhausted yes", first["text"])

            second = optimize(path, budget=5)
            self.assertEqual(len(_runs(log)), 3, "only the missing evaluation runs")
            self.assertEqual(len(second["records"]), 3)
            self.assertEqual(sum(r["reused"] for r in second["records"]), 2)
            self.assertFalse(second["budgetExhausted"])

            third = optimize(path)
            self.assertEqual(len(_runs(log)), 3, "nothing is re-run")
            self.assertEqual(third["newEvaluations"], 0)
            self.assertEqual(
                third["best"]["objective"], second["best"]["objective"]
            )

    def test_config_budget_and_fresh(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            log = _write_solver(folder)
            path = _optimizer_config(folder, samples=3, maxNewEvaluations=1)
            optimize(path)
            self.assertEqual(len(_runs(log)), 1)
            optimize(path, fresh=True)
            self.assertEqual(len(_runs(log)), 2, "fresh discards the stored evaluation")

    def test_failed_evaluation_is_retried_and_objective_edit_is_free(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            log = _write_solver(folder, fail=True)
            path = _optimizer_config(folder, samples=2)
            broken = optimize(path)
            self.assertIsNone(broken["best"])
            self.assertEqual({r["status"] for r in broken["records"]}, {"error"})
            self.assertEqual(len(_runs(log)), 2)

            _write_solver(folder, fail=False)
            fixed = optimize(path)
            self.assertIsNotNone(fixed["best"])
            self.assertEqual(len(_runs(log)), 4, "both failed points ran again")

            config = load_config(path)
            config["raw"]["optimization"]["objective"]["expression"] = "pMeanOutlet"
            path.write_text(_raw_text(config["raw"]), encoding="utf-8")
            edited = optimize(path)
            self.assertEqual(len(_runs(log)), 4, "a new objective needs no solver run")
            self.assertGreater(edited["best"]["objective"], 1.0e5)

    def test_interrupted_evaluation_resumes_in_the_same_folder(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            log = _write_solver(folder, fail=True)
            path = _optimizer_config(folder, samples=1)
            optimize(path)
            case = next((folder / "optimization").glob("eval-*"))
            (case / "0.001").mkdir()
            _write_solver(folder, fail=False)
            result = optimize(path)
            self.assertEqual(result["records"][0]["status"], "ok")
            self.assertIn("resumed from time 0.001", result["records"][0]["message"])
            self.assertEqual(set(_runs(log)), {case.name})
            control = (case / "system" / "controlDict").read_text(encoding="utf-8")
            self.assertIn("startFrom       latestTime;", control)


if __name__ == "__main__":
    unittest.main()

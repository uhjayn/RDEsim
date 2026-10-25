"""Smoke tests that do not need OpenFOAM."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path

from rdesim.acoustics import evaluate_acoustics, helmholtz_resonance_frequency
from rdesim.casewriter import prepare_case
from rdesim.config import ConfigError, get_path, load_config, set_path
from rdesim.expr import evaluate_with_names
from rdesim.foamdict import parse
from rdesim.geometry import build_profile
from rdesim.hyperbolicity import check_hyperbolicity
from rdesim.mesh import build_mesh
from rdesim.mixture import mixture_from_config
from rdesim.optimize import optimize

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "cases" / "annulus2D" / "rdeConfig"


class DictTests(unittest.TestCase):
    def test_nested_comments_and_include(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "piece").write_text("includedValue  3;\n", encoding="utf-8")
            text = textwrap.dedent(
                """
                // comment
                outer
                {
                    a 1;
                    b (2 3);
                    items
                    (
                        { name one; value 4; }
                    );
                    #include "piece"
                }
                """
            )
            data = parse(text, source="sample", base_dir=folder)
            self.assertEqual(data["outer"]["a"], 1)
            self.assertEqual(list(data["outer"]["b"]), [2, 3])
            self.assertEqual(data["outer"]["items"][0]["name"], "one")
            self.assertEqual(data["outer"]["includedValue"], 3)

    def test_example_config_loads(self):
        config = load_config(EXAMPLE)
        self.assertEqual(config["raw"]["fuel"]["name"], "hydrogen-air")
        self.assertEqual(config["fuel"]["fuel"]["name"], "H2")

    def test_missing_key_is_an_error(self):
        config = load_config(EXAMPLE)
        del config["raw"]["domain"]["axialLength"]
        text = _raw_text(config["raw"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "rdeConfig"
            path.write_text(text, encoding="utf-8")
            # The fuel library walk still needs to find fuels/. Put the file
            # next to a copy of the library path by writing under the repo.
        with self.assertRaises(ConfigError):
            load_config(_write_broken())


class FuelTests(unittest.TestCase):
    def test_two_fuels_differ_and_phi_scales_heat_release(self):
        hydrogen = load_config(EXAMPLE)
        state = mixture_from_config(hydrogen)
        self.assertAlmostEqual(state.phi, 1.0)
        self.assertAlmostEqual(state.y_fuel, state.y_fuel_stoich)
        self.assertGreater(state.heat_release, 0.0)
        self.assertGreater(state.sound_speed_squared, 0.0)
        ethylene = load_config(EXAMPLE)
        ethylene["raw"]["fuel"]["name"] = "ethylene-air"
        ethylene = load_config(_materialise(ethylene["raw"]))
        other = mixture_from_config(ethylene)
        self.assertNotEqual(state.mol_weight, other.mol_weight)
        self.assertNotEqual(state.gamma, other.gamma)
        self.assertNotEqual(state.heat_release, other.heat_release)

    def test_library_numbers_are_not_in_the_solver_sources(self):
        needles = ["0.002016", "0.028054", "1.2e8", "4.72e7"]
        roots = [ROOT / "rdesim", ROOT / "applications"]
        for root in roots:
            for path in root.rglob("*"):
                if path.suffix not in {".py", ".C", ".H"}:
                    continue
                text = path.read_text(encoding="utf-8")
                for needle in needles:
                    self.assertNotIn(needle, text, f"{needle} hardcoded in {path}")


class GeometryTests(unittest.TestCase):
    def test_step_angle_and_nozzle(self):
        config = load_config(EXAMPLE)
        config["raw"]["geometry"]["steps"] = []
        config["raw"]["geometry"]["nozzle"]["enabled"] = "no"
        profile = build_profile(config)
        self.assertAlmostEqual(profile.segments[0].r_out0, 0.05)
        config["raw"]["geometry"]["outerWallAngle"] = 5
        tapered = build_profile(config)
        expected = 0.05 + 0.12 * __import__("math").tan(__import__("math").radians(5))
        self.assertAlmostEqual(tapered.segments[-1].r_out1, expected, places=6)

    def test_vertical_step_changes_downstream_radius(self):
        config = load_config(EXAMPLE)
        config["raw"]["geometry"]["nozzle"]["enabled"] = "no"
        profile = build_profile(config)
        upstream = _radius_before(profile, 0.04)
        downstream = _radius_after(profile, 0.04)
        self.assertAlmostEqual(upstream[1], 0.05)
        self.assertAlmostEqual(downstream[1], 0.048)
        self.assertEqual(profile.step_faces[0].side, "outer")
        self.assertFalse(profile.radial_steps_resolved)

    def test_meridional_mesh_contains_the_step(self):
        config = load_config(EXAMPLE)
        config["raw"]["domain"]["topology"] = "meridional"
        config["raw"]["mesh"]["nCircumferential"] = 1
        config["raw"]["mesh"]["nRadial"] = 8
        config["raw"]["geometry"]["nozzle"]["enabled"] = "no"
        config["raw"]["geometry"]["injector"]["recessDepth"] = 0
        config["raw"]["boundaries"]["acoustics"]["targetPatches"] = ["outerWall", "outerStep"]
        profile = build_profile(config)
        mesh = build_mesh(config, profile)
        text = mesh.block_mesh_dict()
        self.assertIn("outerStep", text)
        self.assertIn("0.048", text)
        self.assertIn("0.05", text)
        self.assertTrue(profile.radial_steps_resolved)

    def test_unwrapped_slots_change_with_width(self):
        config = load_config(EXAMPLE)
        profile = build_profile(config)
        mesh = build_mesh(config, profile)
        text = mesh.block_mesh_dict()
        self.assertIn("injector", text)
        self.assertIn("periodicLeft", text)
        self.assertNotIn("outerStep", text)
        config["raw"]["geometry"]["injector"]["slotWidth"] = 0.003
        wider = build_mesh(config, build_profile(config)).block_mesh_dict()
        self.assertNotEqual(text, wider)


class HyperbolicityTests(unittest.TestCase):
    def test_example_is_hyperbolic(self):
        config = load_config(EXAMPLE)
        report = check_hyperbolicity(config, mixture_from_config(config))
        self.assertTrue(report.hyperbolic, report.text())

    def test_threats_are_reported(self):
        config = load_config(EXAMPLE)
        config["fuel"]["thermo"]["gamma"] = 0.9
        report = check_hyperbolicity(config, _state_or_skip(config))
        if report.gamma <= 1:
            self.assertFalse(report.hyperbolic)
        config = load_config(EXAMPLE)
        config["raw"]["numerics"]["fluxScheme"] = "ellipticPressure"
        report = check_hyperbolicity(config, mixture_from_config(config))
        self.assertFalse(report.hyperbolic)
        config = load_config(EXAMPLE)
        config["raw"]["hyperbolicity"]["sourceDependsOnGradients"] = "yes"
        report = check_hyperbolicity(config, mixture_from_config(config))
        self.assertFalse(report.hyperbolic)
        config = load_config(EXAMPLE)
        config["raw"]["numerics"]["flowModel"] = "navierStokes"
        report = check_hyperbolicity(config, mixture_from_config(config))
        self.assertFalse(report.hyperbolic)
        config = load_config(EXAMPLE)
        config["raw"]["hyperbolicity"]["boundaryCoupling"] = "elliptic"
        report = check_hyperbolicity(config, mixture_from_config(config))
        self.assertFalse(report.hyperbolic)


class AcousticTests(unittest.TestCase):
    def test_limits_and_helmholtz_resonance(self):
        config = load_config(EXAMPLE)
        density = 1.2
        speed = 340.0
        acoustics = config["raw"]["boundaries"]["acoustics"]
        acoustics["mode"] = "impedance"
        acoustics["incidenceAngleDeg"] = 0
        acoustics["ZReal"] = density * speed
        acoustics["ZImag"] = 0
        matched = evaluate_acoustics(config, density, speed)
        self.assertAlmostEqual(matched.r_real, 0.0, places=8)
        acoustics["ZReal"] = 0
        acoustics["ZImag"] = 0
        soft = evaluate_acoustics(config, density, speed)
        self.assertAlmostEqual(soft.r_real, -1.0, places=8)
        acoustics["ZReal"] = 1e15
        hard = evaluate_acoustics(config, density, speed)
        self.assertAlmostEqual(hard.magnitude, 1.0, places=6)
        self.assertGreater(hard.r_real, 0.9)
        acoustics["ZReal"] = density * speed
        acoustics["incidenceAngleDeg"] = 60
        oblique = evaluate_acoustics(config, density, speed)
        self.assertAlmostEqual(oblique.r_real, -1.0 / 3.0, places=6)
        acoustics["mode"] = "microstructure"
        acoustics["incidenceAngleDeg"] = 0
        acoustics["microstructure"]["type"] = "helmholtz"
        acoustics["microstructure"]["normalizedResistance"] = 0
        acoustics["frequency"] = helmholtz_resonance_frequency(config, speed)
        resonant = evaluate_acoustics(config, density, speed)
        self.assertAlmostEqual(resonant.magnitude, 1.0, places=5)
        self.assertAlmostEqual(resonant.r_real, -1.0, places=4)
        acoustics["microstructure"]["type"] = "groove"
        acoustics["microstructure"]["grooveWidth"] = 1
        acoustics["microstructure"]["grooveLand"] = 1
        acoustics["microstructure"]["cavityDepth"] = speed / (4.0 * acoustics["frequency"])
        acoustics["microstructure"]["viscousResistance"] = 0
        quarter = evaluate_acoustics(config, density, speed)
        self.assertAlmostEqual(abs(quarter.z_real), 0.0, places=6)


class CaseTests(unittest.TestCase):
    def test_prepare_writes_adjustable_entries(self):
        config = load_config(EXAMPLE)
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp) / "annulus"
            prepare_case(config, case)
            control = (case / "system" / "controlDict").read_text(encoding="utf-8")
            self.assertIn("solver          rdeFluid;", control)
            self.assertIn("maxCo           0.2;", control)
            self.assertIn("lambda", (case / "0" / "lambda").read_text(encoding="utf-8"))
            thermo = (case / "constant" / "physicalProperties").read_text(encoding="utf-8")
            self.assertIn("perfectGas", thermo)
            hydrogen_weight = _mol_weight(thermo)
            config["raw"]["fuel"]["name"] = "ethylene-air"
            switched = load_config(_materialise(config["raw"]))
            prepare_case(switched, Path(tmp) / "ethylene")
            other = _mol_weight(
                (Path(tmp) / "ethylene" / "constant" / "physicalProperties").read_text(encoding="utf-8")
            )
            self.assertNotEqual(hydrogen_weight, other)
            schemes = (case / "system" / "fvSchemes").read_text(encoding="utf-8")
            self.assertIn("Kurganov", schemes)
            self.assertTrue((case / "annulus.foam").is_file())
            vtk = (case / "VTK" / "geometryBlocks.vtk").read_text(encoding="utf-8")
            self.assertIn("UNSTRUCTURED_GRID", vtk)
            config["raw"]["numerics"]["CFL"] = 0.4
            prepare_case(load_config(_materialise(config["raw"])), Path(tmp) / "cfl")
            # The in-memory config still has the original fuel merge. Changing
            # CFL on the loaded object and preparing it is enough.
        config = load_config(EXAMPLE)
        config["raw"]["numerics"]["CFL"] = 0.4
        with tempfile.TemporaryDirectory() as tmp:
            prepare_case(config, Path(tmp) / "cfl")
            control = (Path(tmp) / "cfl" / "system" / "controlDict").read_text(encoding="utf-8")
            self.assertIn("maxCo           0.4;", control)


class OptimizerTests(unittest.TestCase):
    def test_grid_varies_step_height(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            config = load_config(EXAMPLE)
            config["raw"]["optimization"]["method"] = "grid"
            config["raw"]["optimization"]["maxEvaluations"] = 4
            config["raw"]["optimization"]["designVariables"] = [
                {
                    "name": "geometry.steps[0].height",
                    "min": 0.001,
                    "max": 0.002,
                    "samples": 2,
                }
            ]
            config["raw"]["optimization"]["constraints"] = []
            config["raw"]["fuel"]["library"] = str(ROOT / "fuels")
            config["raw"]["optimization"]["solverCommand"] = (
                f"python3 {ROOT / 'tests' / 'stub_solver.py'}"
            )
            config["raw"]["optimization"]["workDirectory"] = "optimization"
            path = folder / "rdeConfig"
            path.write_text(_raw_text(config["raw"]), encoding="utf-8")
            summary = optimize(path)
            self.assertIsNotNone(summary["best"])
            self.assertEqual(len(summary["records"]), 2)
            heights = []
            for record in summary["records"]:
                self.assertEqual(record["status"], "ok", record["message"])
                heights.append(record["objective"])
            self.assertGreater(max(heights), min(heights))
            self.assertIn("height=0.002", " ".join(summary["best"]["values"]))

    def test_expression(self):
        self.assertAlmostEqual(evaluate_with_names("pMeanOutlet/pMeanInlet", {
            "pMeanOutlet": 2.0,
            "pMeanInlet": 4.0,
        }), 0.5)


def _radius_before(profile, y):
    segment = max(
        (item for item in profile.segments if item.y1 <= y + 1e-9),
        key=lambda item: item.y1,
    )
    return segment.r_in1, segment.r_out1


def _radius_after(profile, y):
    segment = min(
        (item for item in profile.segments if item.y0 >= y - 1e-9),
        key=lambda item: item.y0,
    )
    return segment.r_in0, segment.r_out0


def _mol_weight(text: str) -> str:
    for line in text.splitlines():
        if "molWeight" in line:
            return line.split()[-1].rstrip(";")
    raise AssertionError("molWeight missing")


def _state_or_skip(config):
    try:
        return mixture_from_config(config)
    except ConfigError:
        from rdesim.mixture import Mixture
        return Mixture(
            gamma=0.9,
            density=1.0,
            pressure=1.0e5,
            temperature=300.0,
            sound_speed_squared=-1.0,
            equation_of_state="perfectGas",
        )


def _materialise(raw) -> Path:
    path = ROOT / "cases" / "annulus2D" / ".tmp-rdeConfig"
    path.write_text(_raw_text(raw), encoding="utf-8")
    return path


def _write_broken() -> Path:
    config = load_config(EXAMPLE)
    del config["raw"]["operating"]["pInlet"]
    path = ROOT / "cases" / "annulus2D" / ".tmp-broken"
    path.write_text(_raw_text(config["raw"]), encoding="utf-8")
    return path


def _raw_text(raw) -> str:
    from rdesim.foamdict import render_entries
    return render_entries(raw)


def tearDownModule():
    for name in (".tmp-rdeConfig", ".tmp-broken"):
        path = ROOT / "cases" / "annulus2D" / name
        if path.exists():
            path.unlink()


if __name__ == "__main__":
    unittest.main()

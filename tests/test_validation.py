"""Tests for the H2-air validation tooling. They do not need OpenFOAM or Cantera."""

from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

from rdesim.detonation import cj_state, heat_for_cj_speed, von_neumann
from rdesim.validation import (
    Model,
    Setup,
    ValidationError,
    _ref,
    read_field,
    write_case,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "cases" / "detonationTubeH2" / "validationConfig"


def _percent(value: float, reference: float) -> float:
    return 100.0 * (value - reference) / reference


class ReferenceDataTests(unittest.TestCase):
    def setUp(self):
        self.setup = Setup(CONFIG)
        self.ref = _ref(self.setup)

    def test_reference_is_stoichiometric_h2_air_cj(self):
        # Published stoichiometric H2-air CJ speed at 1 atm, 300 K is about 1.97 km/s.
        self.assertGreater(self.ref["D"], 1950.0)
        self.assertLess(self.ref["D"], 1990.0)
        self.assertGreater(self.ref["p"] / 101325.0, 15.0)
        self.assertLess(self.ref["p"] / 101325.0, 18.0)
        self.assertGreater(self.ref["T"], 2800.0)
        self.assertLess(self.ref["T"], 3100.0)

    def test_independent_cea_cj_agrees_with_cantera(self):
        computed = self.setup.reference["computed"]
        cea = float(computed["cea"]["DCJ"])
        self.assertLess(abs(_percent(cea, self.ref["D"])), 0.5)

    def test_cell_width_literature_entries_have_sources(self):
        entries = self.setup.reference["literature"]["cellWidth"]
        self.assertGreaterEqual(len(entries), 3)
        for entry in entries:
            self.assertGreater(float(entry["width"]), 0.0)
            self.assertTrue(str(entry["reference"]).strip())


class TheoryTests(unittest.TestCase):
    def setUp(self):
        self.setup = Setup(CONFIG)

    def test_original_fuel_keeps_its_published_numbers(self):
        # The original fuel is not edited by the calibration.
        text = (ROOT / "fuels" / "hydrogen-air" / "fuelDict").read_text(encoding="utf-8")
        self.assertIn("LHV         1.2e8", text.replace("1.2e+08", "1.2e8"))

    def test_cj_state_is_consistent(self):
        model = Model(self.setup, self.setup.fuel_name("original"))
        cj = cj_state(model.gas)
        # CJ speed gives a sonic burnt gas behind the wave.
        burnt_sound = math.sqrt(model.gas.gamma * cj["p"] / cj["rho"])
        self.assertAlmostEqual((cj["D"] - cj["u"]) / burnt_sound, 1.0, places=6)
        vn = von_neumann(model.gas, cj["D"])
        self.assertGreater(vn["p"], cj["p"])
        self.assertGreater(vn["T"], model.cj["T"] * 0.5)

    def test_heat_for_cj_speed_inverts_cj_state(self):
        model = Model(self.setup, self.setup.fuel_name("original"))
        q = heat_for_cj_speed(model.gas, 1950.0)
        self.assertAlmostEqual(cj_state(model.gas, q)["D"], 1950.0, places=4)

    def test_calibrated_fuel_improves_wave_speed_and_half_release_length(self):
        ref = _ref(self.setup)
        before = Model(self.setup, self.setup.fuel_name("original"))
        after = Model(self.setup, self.setup.fuel_name("calibrated"))
        self.assertLess(abs(_percent(after.cj["D"], ref["D"])), 0.1)
        self.assertLess(abs(_percent(after.znd["x50"], ref["x50"])), 1.0)
        self.assertLess(abs(_percent(after.cj["D"], ref["D"])), abs(_percent(before.cj["D"], ref["D"])))
        self.assertLess(abs(_percent(after.znd["x50"], ref["x50"])), abs(_percent(before.znd["x50"], ref["x50"])))

    def test_calibration_leaves_original_values_in_the_record(self):
        text = (ROOT / "fuels" / "hydrogen-air-validation" / "fuelDict").read_text(encoding="utf-8")
        self.assertIn("calibration", text)
        self.assertIn("original", text)

    def test_single_gamma_pressure_and_temperature_mismatch_is_the_known_one(self):
        ref = _ref(self.setup)
        after = Model(self.setup, self.setup.fuel_name("calibrated"))
        self.assertLess(_percent(after.cj["p"], ref["p"]), -5.0)
        self.assertLess(_percent(after.cj["T"], ref["T"]), -5.0)


class FieldReaderTests(unittest.TestCase):
    def _write(self, text: str) -> Path:
        handle = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8")
        handle.write(text)
        handle.close()
        self.addCleanup(Path(handle.name).unlink)
        return Path(handle.name)

    def test_scalar_list(self):
        path = self._write("internalField   nonuniform List<scalar> 3(1 2.5 3e2);\nboundaryField {}\n")
        self.assertEqual(read_field(path, 3), [1.0, 2.5, 300.0])

    def test_vector_list_returns_x(self):
        path = self._write("internalField   nonuniform List<vector> 2((1 0 0)(-2.5 0 0));\n")
        self.assertEqual(read_field(path, 2), [1.0, -2.5])

    def test_uniform(self):
        path = self._write("internalField   uniform 5;\n")
        self.assertEqual(read_field(path, 4), [5.0] * 4)

    def test_wrong_count_is_an_error(self):
        path = self._write("internalField   nonuniform List<scalar> 2(1 2);\n")
        with self.assertRaises(ValidationError):
            read_field(path, 3)


class CaseWriterTests(unittest.TestCase):
    def test_case_has_resolved_zone_and_required_files(self):
        setup = Setup(CONFIG)
        run = next(r for r in setup.raw["runs"] if str(r["name"]) == "calibrated-N10")
        model = Model(setup, setup.fuel_name("calibrated"))
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "case"
            meta = write_case(setup, run, model, out)
            for name in (
                "system/controlDict", "system/blockMeshDict", "system/fvSchemes", "system/fvSolution",
                "constant/physicalProperties", "constant/rdeProperties", "0/p", "0/T", "0/U", "0/lambda",
            ):
                self.assertTrue((out / name).is_file(), name)
            self.assertAlmostEqual(model.x_half / meta["dx"], 10.0, places=6)
            self.assertEqual(meta["cells"], round(meta["length"] / meta["dx"]))


if __name__ == "__main__":
    unittest.main()

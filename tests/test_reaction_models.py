"""Reaction model selection, the two-step case files, and the 1D CJ check helpers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from rdesim.casewriter import prepare_case
from rdesim.cjcheck import (
    cj_perfect_gas,
    fitted_speed,
    front_position,
    read_scalar_field,
    tube_block_mesh,
)
from rdesim.config import ConfigError, load_config, reaction_model
from rdesim.runs import fingerprint
from tests.test_rdesim import EXAMPLE


class ModelSelectionTest(unittest.TestCase):
    def test_default_is_one_step(self):
        self.assertEqual(reaction_model(load_config(EXAMPLE)), "oneStepProgress")

    def test_model_override_selects_two_step(self):
        config = load_config(EXAMPLE, "twoStep")
        self.assertEqual(reaction_model(config), "twoStep")

    def test_unknown_model_is_rejected(self):
        with self.assertRaises(ConfigError) as caught:
            load_config(EXAMPLE, "fiveStep")
        self.assertIn("reaction.model", str(caught.exception))

    def test_fingerprint_separates_models(self):
        self.assertNotEqual(
            fingerprint(load_config(EXAMPLE)), fingerprint(load_config(EXAMPLE, "twoStep"))
        )

    def test_two_step_needs_positive_rates(self):
        config = load_config(EXAMPLE, "twoStep")
        self.assertEqual(float(config["fuel"]["twoStep"]["releaseA"]), 8.0e6)


class TwoStepCaseTest(unittest.TestCase):
    def test_one_step_properties_keep_their_entries(self):
        with tempfile.TemporaryDirectory() as folder:
            prepare_case(load_config(EXAMPLE), Path(folder) / "case")
            text = (Path(folder) / "case" / "constant" / "rdeProperties").read_text()
            self.assertIn("reactionModel      oneStepProgress;", text)
            for key in ("A ", "beta", "Ta ", "TRef", "ignitionTemperature", "TFloor"):
                self.assertIn(f"    {key}", text)
            self.assertNotIn("inductionA", text)
            self.assertFalse((Path(folder) / "case" / "0" / "xi").exists())

    def test_two_step_writes_xi_and_its_settings(self):
        with tempfile.TemporaryDirectory() as folder:
            case = Path(folder) / "case"
            prepare_case(load_config(EXAMPLE, "twoStep"), case)
            properties = (case / "constant" / "rdeProperties").read_text()
            self.assertIn("reactionModel      twoStep;", properties)
            for key in ("inductionA", "inductionTa", "releaseA", "releaseTa", "releaseOrder"):
                self.assertIn(key, properties)
            self.assertNotIn("ignitionTemperature", properties)
            self.assertTrue((case / "0" / "xi").is_file())
            self.assertIn("reconstruct(xi)", (case / "system" / "fvSchemes").read_text())
            solution = (case / "system" / "fvSolution").read_text()
            self.assertIn("xiFinal", solution)
            report = (case / "constant" / "hyperbolicityReport").read_text()
            self.assertIn("reaction model        twoStep", report)

    def test_two_step_checkpoint_keeps_xi(self):
        config = load_config(EXAMPLE, "twoStep")
        config["raw"]["output"]["checkpoint"] = {
            "enabled": True,
            "writeControl": "timeStep",
            "writeInterval": 10,
        }
        with tempfile.TemporaryDirectory() as folder:
            prepare_case(config, Path(folder) / "case")
            control = (Path(folder) / "case" / "system" / "controlDict").read_text()
            self.assertIn("lambda xi", control)


class CjHelpersTest(unittest.TestCase):
    def test_cj_for_no_heat_release_is_sound_speed(self):
        result = cj_perfect_gas(1.4, 0.0, 340.0, 101325.0)
        self.assertAlmostEqual(result["mach"], 1.0)
        self.assertAlmostEqual(result["pressureRatio"], 1.0)

    def test_cj_known_value(self):
        # H = (gamma^2 - 1) q / (2 c0^2) = 1 gives M = sqrt(2) + 1.
        gamma, c0 = 1.4, 400.0
        q = 2.0 * c0 * c0 / (gamma * gamma - 1.0)
        result = cj_perfect_gas(gamma, q, c0, 1.0e5)
        self.assertAlmostEqual(result["mach"], 2.0 ** 0.5 + 1.0, places=12)
        self.assertAlmostEqual(
            result["pressureRatio"], (1 + gamma * result["mach"] ** 2) / (gamma + 1), places=12
        )

    def test_cj_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            cj_perfect_gas(1.0, 1.0, 340.0, 1.0e5)

    def test_front_position_interpolates_the_last_crossing(self):
        values = [10.0, 10.0, 8.0, 2.0, 0.0, 0.0]
        self.assertAlmostEqual(front_position(values, 1.0, 5.0), 3.0, places=12)
        self.assertIsNone(front_position([0.0, 0.0], 1.0, 5.0))

    def test_fitted_speed_is_the_slope(self):
        times = [0.0, 1.0, 2.0, 3.0]
        self.assertAlmostEqual(fitted_speed(times, [5.0 + 7.0 * t for t in times]), 7.0)

    def test_read_scalar_field_reads_lists_and_uniform(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "f"
            path.write_text("internalField   nonuniform List<scalar>\n3\n(\n1\n2\n3.5\n)\n;\n")
            self.assertEqual(read_scalar_field(path), [1.0, 2.0, 3.5])
            path.write_text("internalField   uniform 4.5;\n")
            self.assertEqual(read_scalar_field(path), 4.5)

    def test_tube_mesh_names_patches_the_fields_know(self):
        text = tube_block_mesh(0.1, 100)
        for name in ("inletWall", "outlet", "front", "back"):
            self.assertIn(name, text)
        self.assertIn("(100 1 1)", text)


if __name__ == "__main__":
    unittest.main()

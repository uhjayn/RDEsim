"""Detailed multi-species chemistry: configuration, mechanism data, and case files."""

from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

from rdesim.casewriter import prepare_case
from rdesim.config import ConfigError, load_config, reaction_model
from rdesim.hyperbolicity import check_hyperbolicity
from rdesim.mixture import mixture_from_config
from rdesim.runs import fingerprint
from rdesim.species import load_mechanism, mechanism_directory, unburned_mass_fractions
from tests.test_rdesim import EXAMPLE, ROOT, _materialise

SPECIES = ["H", "H2", "O", "OH", "H2O", "O2", "HO2", "H2O2", "N2"]
MECHANISM = ROOT / "fuels" / "hydrogen-air" / "detailed" / "burke2012"


def _detailed():
    return load_config(EXAMPLE, "detailed")


class MechanismDataTest(unittest.TestCase):
    def test_species_and_ranges(self):
        mechanism = load_mechanism(MECHANISM)
        self.assertEqual(mechanism.names, SPECIES)
        for specie in mechanism.species:
            self.assertEqual(specie.t_low, 200.0, specie.name)
            self.assertGreaterEqual(specie.t_high, 3500.0, specie.name)

    def test_reactions_file_has_the_23_published_reactions(self):
        text = (MECHANISM / "reactions").read_text()
        reversible = re.findall(r"type\s+reversible\w+;", text)
        # Falloff and third-body reactions use their own types; count them all.
        types = re.findall(r"\n\s+type\s+(\w+);", text)
        self.assertTrue(reversible)
        self.assertGreaterEqual(len(types), 23)
        self.assertNotIn("AR", text)
        self.assertNotIn("HE", text)

    def test_molecular_weights_are_standard(self):
        mechanism = load_mechanism(MECHANISM)
        self.assertAlmostEqual(mechanism.specie("H2").weight, 2.016, places=3)
        self.assertAlmostEqual(mechanism.specie("O2").weight, 31.9988, places=3)
        self.assertAlmostEqual(mechanism.specie("N2").weight, 28.0134, places=4)

    def test_unknown_specie_is_named(self):
        with self.assertRaises(ConfigError) as caught:
            load_mechanism(MECHANISM).specie("CH4")
        self.assertIn("CH4", str(caught.exception))


class DetailedConfigTest(unittest.TestCase):
    def test_model_is_selectable(self):
        self.assertEqual(reaction_model(_detailed()), "detailed")

    def test_unburned_stoichiometric_mixture(self):
        config = _detailed()
        fractions = unburned_mass_fractions(config, load_mechanism(mechanism_directory(config)))
        self.assertAlmostEqual(sum(fractions.values()), 1.0, places=12)
        self.assertAlmostEqual(fractions["H2"], 0.02851, places=4)
        self.assertEqual(fractions["H"], 0.0)

    def test_frozen_mixture_properties(self):
        state = mixture_from_config(_detailed())
        self.assertTrue(1.38 < state.gamma < 1.42)
        self.assertTrue(395.0 < state.sound_speed < 420.0, state.sound_speed)

    def test_fingerprint_separates_the_detailed_model(self):
        self.assertNotEqual(fingerprint(load_config(EXAMPLE)), fingerprint(_detailed()))

    def _refused(self, edit) -> str:
        config = load_config(EXAMPLE)
        edit(config["raw"])
        path = _materialise(config["raw"])
        try:
            with self.assertRaises(ConfigError) as caught:
                load_config(path, "detailed")
        finally:
            path.unlink()
        return str(caught.exception)

    def test_navier_stokes_is_refused(self):
        message = self._refused(lambda raw: raw["numerics"].update(flowModel="navierStokes"))
        self.assertIn("euler only", message)

    def test_progress_kernel_is_refused(self):
        message = self._refused(lambda raw: raw["ignition"].update(type="progressKernel"))
        self.assertIn("no progress variable", message)

    def test_hyperbolicity_check_runs_for_species(self):
        config = _detailed()
        report = check_hyperbolicity(config, mixture_from_config(config))
        self.assertTrue(report.hyperbolic, report.text())
        self.assertIn("detailed", report.text())
        self.assertIn("not known before the run", report.text())

    def test_kernel_above_janaf_range_fails_the_check(self):
        config = _detailed()
        config["raw"]["ignition"]["TKernel"] = 6000
        report = check_hyperbolicity(config, mixture_from_config(config))
        self.assertFalse(report.hyperbolic)
        self.assertIn("janaf range", report.text())


class DetailedCaseTest(unittest.TestCase):
    def test_case_files(self):
        with tempfile.TemporaryDirectory() as folder:
            case = Path(folder) / "case"
            prepare_case(_detailed(), case)
            for name in SPECIES:
                self.assertTrue((case / "0" / name).is_file(), name)
            properties = (case / "constant" / "rdeProperties").read_text()
            self.assertIn("reactionModel      detailed;", properties)
            thermo = (case / "constant" / "physicalProperties").read_text()
            self.assertIn("multicomponentMixture", thermo)
            self.assertIn("janaf", thermo)
            self.assertTrue((case / "constant" / "thermo").is_file())
            self.assertTrue((case / "constant" / "reactions").is_file())
            chemistry = (case / "constant" / "chemistryProperties").read_text()
            for key in ("TMin", "initialChemicalTimeStep", "fuelSpecie", "seulex"):
                self.assertIn(key, chemistry)
            self.assertIn('#include "reactions"', chemistry)
            h2 = (case / "0" / "H2").read_text()
            self.assertIn("uniform 0.0285", h2)

    def test_checkpoint_keeps_the_species(self):
        config = _detailed()
        config["raw"]["output"]["checkpoint"] = {
            "enabled": True,
            "writeControl": "timeStep",
            "writeInterval": 10,
        }
        with tempfile.TemporaryDirectory() as folder:
            prepare_case(config, Path(folder) / "case")
            control = (Path(folder) / "case" / "system" / "controlDict").read_text()
            self.assertIn("H2O2", control)

    def test_one_step_case_has_no_species_files(self):
        with tempfile.TemporaryDirectory() as folder:
            case = Path(folder) / "case"
            prepare_case(load_config(EXAMPLE), case)
            self.assertFalse((case / "0" / "H2").exists())
            self.assertFalse((case / "constant" / "chemistryProperties").exists())


class ChemcheckHelperTest(unittest.TestCase):
    def test_ignition_delay_is_the_steepest_interval(self):
        from rdesim.chemcheck import ignition_delay

        times = [0.0, 1.0, 2.0, 3.0, 4.0]
        temperatures = [1000.0, 1010.0, 1100.0, 2500.0, 2510.0]
        self.assertAlmostEqual(ignition_delay(times, temperatures), 2.5)
        self.assertIsNone(ignition_delay([0.0, 1.0], [1000.0, 1000.0]))


class EquilibriumReferenceTest(unittest.TestCase):
    def test_cantera_cj_speed_for_stoichiometric_h2_air(self):
        try:
            import cantera  # noqa: F401
        except ImportError:
            self.skipTest("cantera is not installed")
        from rdesim.cjcheck import cj_equilibrium

        config = _detailed()
        config["raw"]["ignition"]["initialTemperature"] = 300.0
        result = cj_equilibrium(config)
        # Published stoichiometric H2-air CJ speeds at 1 atm and 300 K are
        # close to 1970 m/s.
        self.assertTrue(1950.0 < result["speed"] < 2000.0, result["speed"])
        self.assertTrue(2900.0 < result["temperature"] < 3020.0, result["temperature"])


if __name__ == "__main__":
    unittest.main()

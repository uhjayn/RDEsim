"""Multi-species thermodynamics read from an OpenFOAM mechanism ``thermo`` file.

The detailed reaction model carries the species of a published mechanism. This
module reads the same ``thermo`` file that OpenFOAM reads (molWeight and the
NASA-7 janaf coefficients per specie, kg/kmol and J/kg/K) and gives the Python
side what it needs without OpenFOAM: the unburned mass fractions for an
equivalence ratio, frozen mixture properties (gas constant, cp, cv, gamma,
sound speed) and a coverage check of the janaf temperature range. Nothing is
invented here; every number is in the file.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from .config import ConfigError
from .foamdict import parse_file

RR = 8314.462618  # J/(kmol K), the OpenFOAM universal gas constant


@dataclass
class Specie:
    name: str
    weight: float  # kg/kmol
    t_low: float
    t_high: float
    t_common: float
    high: list[float]
    low: list[float]
    elements: dict

    def _coefficients(self, temperature: float) -> list[float]:
        return self.high if temperature >= self.t_common else self.low

    def cp(self, temperature: float) -> float:
        """Heat capacity at constant pressure, J/(kg K)."""
        a = self._coefficients(temperature)
        t = temperature
        return RR / self.weight * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4]))))

    def cv(self, temperature: float) -> float:
        return self.cp(temperature) - RR / self.weight

    def in_range(self, temperature: float) -> bool:
        return self.t_low <= temperature <= self.t_high


class Mechanism:
    """Species of a mechanism, in the order of the ``species`` list of the file."""

    def __init__(self, directory: Path, species: list[Specie]):
        self.directory = directory
        self.species = species
        self.index = {specie.name: i for i, specie in enumerate(species)}

    @property
    def names(self) -> list[str]:
        return [specie.name for specie in self.species]

    def specie(self, name: str) -> Specie:
        if name not in self.index:
            raise ConfigError(
                f"Species '{name}' is not in the mechanism at {self.directory}. "
                f"It has: {', '.join(self.names)}"
            )
        return self.species[self.index[name]]

    def mixture(self, mass_fractions: dict[str, float], temperature: float) -> dict:
        """Frozen properties of a mixture at one temperature (perfect gas)."""
        total = sum(mass_fractions.values())
        if abs(total - 1.0) > 1e-9:
            raise ConfigError(f"Mass fractions sum to {total:.12g}, not 1")
        inverse_weight = 0.0
        cp = 0.0
        for name, fraction in mass_fractions.items():
            specie = self.specie(name)
            inverse_weight += fraction / specie.weight
            cp += fraction * specie.cp(temperature)
        gas_constant = RR * inverse_weight
        cv = cp - gas_constant
        gamma = cp / cv
        return {
            "weight": 1.0 / inverse_weight,
            "gasConstant": gas_constant,
            "cp": cp,
            "cv": cv,
            "gamma": gamma,
            "soundSpeed": math.sqrt(gamma * gas_constant * temperature),
        }


def mechanism_directory(config: dict) -> Path:
    """Directory named by ``detailed.mechanism`` in the fuelDict, relative to it."""
    detailed = config["fuel"].get("detailed")
    if not isinstance(detailed, dict) or "mechanism" not in detailed:
        raise ConfigError("reaction.model is detailed but the fuel has no detailed.mechanism")
    return (Path(config["fuelPath"]).parent / str(detailed["mechanism"])).resolve()


def load_mechanism(directory: str | Path) -> Mechanism:
    directory = Path(directory)
    thermo_path = directory / "thermo"
    if not thermo_path.is_file():
        raise ConfigError(f"Mechanism thermo file not found: {thermo_path}")
    if not (directory / "reactions").is_file():
        raise ConfigError(f"Mechanism reactions file not found: {directory / 'reactions'}")
    data = parse_file(thermo_path)
    listing = data.get("species")
    if not isinstance(listing, list) or len(listing) != 2:
        raise ConfigError(f"{thermo_path}: expected 'species N ( ... );'")
    names = [str(item) for item in listing[1]]
    species = []
    for name in names:
        entry = data.get(name)
        if not isinstance(entry, dict):
            raise ConfigError(f"{thermo_path}: no dictionary for specie {name}")
        try:
            thermodynamics = entry["thermodynamics"]
            species.append(
                Specie(
                    name=name,
                    weight=float(entry["specie"]["molWeight"]),
                    t_low=float(thermodynamics["Tlow"]),
                    t_high=float(thermodynamics["Thigh"]),
                    t_common=float(thermodynamics["Tcommon"]),
                    high=[float(v) for v in thermodynamics["highCpCoeffs"]],
                    low=[float(v) for v in thermodynamics["lowCpCoeffs"]],
                    elements=dict(entry.get("elements", {})),
                )
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ConfigError(f"{thermo_path}: specie {name} is incomplete ({exc})") from exc
    return Mechanism(directory, species)


def unburned_mass_fractions(config: dict, mechanism: Mechanism) -> dict[str, float]:
    """Mass fractions at ``operating.equivalenceRatio`` from the fuelDict stoichiometry.

    Per mole of oxidizer mixture the oxidizer carries its own mole fractions and
    the fuel is ``phi * nuFuel / nuOxygen`` times the O2 present. Molecular
    weights are the mechanism's.
    """
    fuel = config["fuel"]
    phi = float(config["raw"]["operating"]["equivalenceRatio"])
    nu_fuel = float(fuel["stoichiometry"]["nuFuel"])
    nu_oxygen = float(fuel["stoichiometry"]["nuOxygen"])
    fuel_name = str(fuel["fuel"]["name"])
    oxidizer = fuel["oxidizer"]
    moles: dict[str, float] = {}
    for item in oxidizer["species"]:
        name = str(item["name"])
        mechanism.specie(name)
        moles[name] = moles.get(name, 0.0) + float(item["moleFraction"])
    if "O2" not in moles:
        raise ConfigError("The detailed model needs O2 in fuel.oxidizer.species")
    mechanism.specie(fuel_name)
    # Per mole of oxidizer mixture: x_O2 moles of O2, and phi times the
    # stoichiometric fuel for that oxygen.
    moles[fuel_name] = moles.get(fuel_name, 0.0) + phi * nu_fuel / nu_oxygen * moles["O2"]

    mass = {name: n * mechanism.specie(name).weight for name, n in moles.items()}
    total = sum(mass.values())
    fractions = {name: value / total for name, value in mass.items()}
    return {name: fractions.get(name, 0.0) for name in mechanism.names}

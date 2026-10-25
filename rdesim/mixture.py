"""Mixture state from a fuel dictionary and the operating-condition entries.

Molecular weights and the lower heating value come from the fuel file. The
equivalence ratio, pressure, and temperature come from the case. The universal
gas constant comes from ``constants.RUniversal``.
"""

from __future__ import annotations

import math

from .config import ConfigError, as_bool


class Mixture:
    def __init__(self, **values):
        self.__dict__.update(values)


def mixture_from_config(config: dict) -> Mixture:
    fuel = config["fuel"]
    operating = config["raw"]["operating"]
    gas_constant = float(config["raw"]["constants"]["RUniversal"])
    phi = float(operating["equivalenceRatio"])
    if phi <= 0:
        raise ConfigError("'operating.equivalenceRatio' must be positive")

    fuel_weight = float(fuel["fuel"]["molWeight"])
    nu_fuel = float(fuel["stoichiometry"]["nuFuel"])
    nu_oxygen = float(fuel["stoichiometry"]["nuOxygen"])
    oxygen_fraction = float(fuel["oxidizer"]["oxygenMoleFraction"])
    if fuel_weight <= 0 or nu_fuel <= 0 or nu_oxygen <= 0:
        raise ConfigError("Fuel stoichiometric coefficients and molWeight must be positive")
    if oxygen_fraction <= 0 or oxygen_fraction > 1:
        raise ConfigError("'fuel.oxidizer.oxygenMoleFraction' must lie in (0, 1]")

    oxidizer_weight = 0.0
    fraction_sum = 0.0
    for specie in fuel["oxidizer"]["species"]:
        fraction = float(specie["moleFraction"])
        weight = float(specie["molWeight"])
        if weight <= 0:
            raise ConfigError(f"Oxidizer species '{specie['name']}' has a non-positive molWeight")
        oxidizer_weight += fraction * weight
        fraction_sum += fraction
    if oxidizer_weight <= 0:
        raise ConfigError("Oxidizer mixture molWeight is not positive")
    if abs(fraction_sum - 1.0) > 1e-6:
        raise ConfigError(
            "Oxidizer mole fractions must sum to 1. "
            f"The fuel file sums to {fraction_sum:.8g}."
        )

    # Moles of oxidizer mixture per mole of fuel at stoichiometry.
    oxidizer_moles = nu_oxygen / (nu_fuel * oxygen_fraction)
    mass_fuel = fuel_weight
    mass_oxidizer = oxidizer_moles * oxidizer_weight
    y_fuel_stoich = mass_fuel / (mass_fuel + mass_oxidizer)
    y_fuel = (phi * y_fuel_stoich) / (1.0 + (phi - 1.0) * y_fuel_stoich)
    if y_fuel <= 0 or y_fuel >= 1:
        raise ConfigError(
            f"Equivalence ratio {phi} produced fuel mass fraction {y_fuel}, "
            "which is outside (0, 1)"
        )

    thermo = fuel["thermo"]
    composition = thermo["composition"]
    if composition == "frozenReactants":
        mol_weight = 1.0 / (y_fuel / fuel_weight + (1.0 - y_fuel) / oxidizer_weight)
        gamma = float(thermo["gamma"])
    elif composition == "frozenProducts":
        mol_weight = float(thermo["productMolWeight"])
        gamma = float(thermo["productGamma"])
    else:
        raise ConfigError(f"Unknown fuel.thermo.composition '{composition}'")

    if mol_weight <= 0:
        raise ConfigError("Mixture molWeight is not positive")
    if gamma <= 1.0:
        raise ConfigError(
            f"Mixture gamma is {gamma}. A perfect gas with gamma <= 1 has "
            "c^2 <= 0 and is not hyperbolic."
        )

    specific_gas = gas_constant / mol_weight
    cp = gamma / (gamma - 1.0) * specific_gas
    cv = cp / gamma
    lhv = float(fuel["heatRelease"]["LHV"])
    heat_release = y_fuel * lhv
    pressure = float(operating["pInlet"])
    temperature = float(operating["TInlet"])
    if pressure <= 0 or temperature <= 0:
        raise ConfigError("Inlet pressure and temperature must be positive")
    density = pressure / (specific_gas * temperature)
    sound_speed_squared = gamma * pressure / density
    sound_speed = math.sqrt(sound_speed_squared)

    flow_model = config["raw"]["numerics"]["flowModel"]
    mu = 0.0 if flow_model == "euler" else float(fuel["transport"]["mu"])
    prandtl = float(fuel["transport"]["Pr"])

    return Mixture(
        phi=phi,
        y_fuel=y_fuel,
        y_fuel_stoich=y_fuel_stoich,
        mol_weight=mol_weight,
        gamma=gamma,
        specific_gas=specific_gas,
        cp=cp,
        cv=cv,
        heat_release=heat_release,
        lhv=lhv,
        density=density,
        pressure=pressure,
        temperature=temperature,
        sound_speed=sound_speed,
        sound_speed_squared=sound_speed_squared,
        mu=mu,
        prandtl=prandtl,
        equation_of_state=str(thermo["equationOfState"]),
        composition=composition,
        hf=float(thermo["hf"]),
        chemistry=fuel["chemistry"],
    )


def inlet_speed(config: dict, state: Mixture, area: float) -> float:
    """Physical inlet speed from the selected inlet model. ``area`` is m^2."""
    operating = config["raw"]["operating"]
    model = operating["inletModel"]
    if area <= 0:
        raise ConfigError("Injector area is not positive, so the inlet speed is undefined")
    if model == "fixedVelocity":
        velocity = operating["inletVelocity"]
        circumferential = float(velocity["circumferential"])
        axial = float(velocity["axial"])
        radial = float(velocity["radial"])
        return math.sqrt(circumferential**2 + axial**2 + radial**2)
    if model == "massFlow":
        mass_flow = float(operating["massFlowRate"])
        return mass_flow / (state.density * area)
    if model == "isentropicPlenum":
        return _isentropic_speed(config, state)
    raise ConfigError(f"Unknown operating.inletModel '{model}'")


def _isentropic_speed(config: dict, state: Mixture) -> float:
    operating = config["raw"]["operating"]
    p0 = float(operating["plenumPressure"])
    t0 = float(operating["plenumTemperature"])
    pressure = float(operating["pInlet"])
    if p0 <= 0 or t0 <= 0:
        raise ConfigError("Plenum pressure and temperature must be positive")
    if pressure > p0:
        raise ConfigError(
            "operating.pInlet is above operating.plenumPressure, so the "
            "isentropic plenum model has no real injection speed. Lower pInlet "
            "or raise plenumPressure."
        )
    gamma = state.gamma
    # Critical pressure ratio of a perfect gas, from gamma in the fuel file.
    critical = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    if pressure / p0 <= critical:
        throat_temperature = t0 * 2.0 / (gamma + 1.0)
        return math.sqrt(gamma * state.specific_gas * throat_temperature)
    cp = state.cp
    speed_squared = 2.0 * cp * t0 * (1.0 - (pressure / p0) ** ((gamma - 1.0) / gamma))
    if speed_squared < 0:
        raise ConfigError("Isentropic plenum model produced a negative speed squared")
    return math.sqrt(speed_squared)


def velocity_components(config: dict, speed: float) -> dict[str, float]:
    """Direction from ``injectionAngle`` when the model sets a speed.

    ``fixedVelocity`` keeps the three components the user wrote. Other models
    place ``speed`` in the axial-circumferential plane at ``injectionAngle``
    degrees from the axial direction. The radial component stays 0.
    """
    operating = config["raw"]["operating"]
    if operating["inletModel"] == "fixedVelocity":
        velocity = operating["inletVelocity"]
        return {
            "circumferential": float(velocity["circumferential"]),
            "axial": float(velocity["axial"]),
            "radial": float(velocity["radial"]),
        }
    angle = math.radians(float(config["raw"]["geometry"]["injector"]["injectionAngle"]))
    return {
        "circumferential": speed * math.sin(angle),
        "axial": speed * math.cos(angle),
        "radial": 0.0,
    }


def mesh_velocity(components: dict[str, float], topology: str, dimension: int, injection: bool) -> tuple[float, float, float]:
    circumferential = components["circumferential"]
    axial = components["axial"]
    radial = components["radial"]
    if dimension == 3 and abs(circumferential) > 0:
        raise ConfigError(
            "A uniform 3D vector cannot represent circumferential velocity. "
            "Set geometry.injector.injectionAngle and operating.inletVelocity.circumferential "
            "to 0 for a 3D mesh, or use the unwrapped topology."
        )
    if topology == "unwrapped":
        return (circumferential, axial, radial)
    if topology == "meridional":
        return (axial, radial, circumferential)
    return (0.0, 0.0, axial)


def uses_viscosity(config: dict) -> bool:
    return config["raw"]["numerics"]["flowModel"] == "navierStokes" and not as_bool(
        config["raw"]["numerics"].get("flowModel") == "euler"
    )

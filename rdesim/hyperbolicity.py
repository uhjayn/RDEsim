"""Report whether the selected model stays a hyperbolic balance law.

The homogeneous reactive Euler equations with a perfect gas are hyperbolic
when density, pressure, and c^2 are positive. A local chemistry source does
not change that. Viscous terms, a steady ddt scheme, a non-local source, or
an equation of state outside the user list are reported as threats.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import reaction_model
from .mixture import Mixture


IMPLEMENTED_FLUX = ("Kurganov", "Tadmor")
IMPLEMENTED_SOURCE = ("pointImplicit", "explicit")
IMPLEMENTED_EOS = ("perfectGas",)

# The reaction models only add transported scalars with local sources. They
# leave the Euler flux Jacobian, and so the frozen sound speed, unchanged.
REACTION_NOTES = {
    "oneStepProgress": "",
    "twoStep": (
        "twoStep adds xi and lambda as passive scalars with local sources; "
        "the Euler flux Jacobian and frozen c^2 are those of the single perfect gas."
    ),
}


@dataclass
class HyperbolicityReport:
    hyperbolic: bool
    failures: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    gamma: float = 0.0
    density: float = 0.0
    pressure: float = 0.0
    temperature: float = 0.0
    sound_speed_squared: float = 0.0
    flux_scheme: str = ""
    flow_model: str = ""
    reaction_model: str = "oneStepProgress"
    reaction_note: str = ""

    def text(self) -> str:
        state = "hyperbolic" if self.hyperbolic else "NOT hyperbolic"
        lines = [
            f"Hyperbolicity: {state}",
            f"  gamma                 {self.gamma:.8g}",
            f"  density               {self.density:.8g} kg/m^3",
            f"  pressure              {self.pressure:.8g} Pa",
            f"  temperature           {self.temperature:.8g} K",
            f"  c^2                   {self.sound_speed_squared:.8g} m^2/s^2",
            f"  flux scheme           {self.flux_scheme}",
            f"  flow model            {self.flow_model}",
            f"  reaction model        {self.reaction_model}",
        ]
        if self.reaction_note:
            lines.append(f"  {self.reaction_note}")
        if self.failures:
            lines.append("Failures:")
            lines.extend(f"  - {item}" for item in self.failures)
        if self.warnings:
            lines.append("Warnings:")
            lines.extend(f"  - {item}" for item in self.warnings)
        if not self.failures and not self.warnings:
            lines.append("The selected equation of state, source, and scheme satisfy the checks.")
        return "\n".join(lines)


def check_hyperbolicity(config: dict, state: Mixture) -> HyperbolicityReport:
    raw = config["raw"]
    numerics = raw["numerics"]
    hyper = raw["hyperbolicity"]
    boundaries = raw["boundaries"]
    report = HyperbolicityReport(
        hyperbolic=True,
        gamma=state.gamma,
        density=state.density,
        pressure=state.pressure,
        temperature=state.temperature,
        sound_speed_squared=state.sound_speed_squared,
        flux_scheme=str(numerics["fluxScheme"]),
        flow_model=str(numerics["flowModel"]),
        reaction_model=reaction_model(config),
        reaction_note=REACTION_NOTES.get(reaction_model(config), ""),
    )

    def fail(message: str) -> None:
        report.failures.append(message)
        report.hyperbolic = False

    if state.gamma <= float(hyper["gammaMin"]):
        fail(
            f"gamma {state.gamma:.8g} is not above hyperbolicity.gammaMin "
            f"{float(hyper['gammaMin']):.8g}"
        )
    if state.density <= float(hyper["rhoMin"]):
        fail("Inlet density is not above hyperbolicity.rhoMin")
    if state.pressure <= float(hyper["pMin"]):
        fail("Inlet pressure is not above hyperbolicity.pMin")
    if state.temperature <= float(hyper["TMin"]):
        fail("Inlet temperature is not above hyperbolicity.TMin")
    if state.sound_speed_squared <= float(hyper["c2Min"]):
        fail("c^2 is not above hyperbolicity.c2Min, so the flux Jacobian is not real")

    eos = state.equation_of_state
    allowed_eos = [str(item) for item in hyper["hyperbolicEquationsOfState"]]
    if eos not in allowed_eos:
        fail(
            f"Equation of state '{eos}' is not in "
            "hyperbolicity.hyperbolicEquationsOfState"
        )
    if eos not in IMPLEMENTED_EOS:
        fail(
            f"Equation of state '{eos}' is not implemented. "
            f"This version writes {', '.join(IMPLEMENTED_EOS)}."
        )

    scheme = str(numerics["fluxScheme"])
    allowed_flux = [str(item) for item in hyper["hyperbolicFluxSchemes"]]
    if scheme not in allowed_flux:
        fail(f"Flux scheme '{scheme}' is not in hyperbolicity.hyperbolicFluxSchemes")
    if scheme not in IMPLEMENTED_FLUX:
        fail(
            f"Flux scheme '{scheme}' is not implemented by rdeFluid. "
            f"Implemented schemes: {', '.join(IMPLEMENTED_FLUX)}."
        )

    source = str(numerics["sourceIntegration"])
    allowed_source = [str(item) for item in hyper["localSourceIntegrations"]]
    if source not in allowed_source:
        fail(
            f"Source integration '{source}' is not in "
            "hyperbolicity.localSourceIntegrations"
        )
    if source not in IMPLEMENTED_SOURCE:
        fail(
            f"Source integration '{source}' is not implemented. "
            f"Implemented: {', '.join(IMPLEMENTED_SOURCE)}."
        )

    if _yes(hyper["sourceDependsOnGradients"]):
        fail(
            "hyperbolicity.sourceDependsOnGradients is yes. A source that depends "
            "on spatial derivatives is not a hyperbolic balance law."
        )
    coupling = str(hyper["boundaryCoupling"])
    if coupling != "local":
        fail(
            f"hyperbolicity.boundaryCoupling is '{coupling}'. "
            "Only 'local' keeps the boundary treatment from adding an elliptic coupling."
        )

    strict = _yes(hyper["requireStrictlyHyperbolic"])
    if numerics["flowModel"] == "navierStokes":
        if strict:
            fail(
                "numerics.flowModel is navierStokes and "
                "hyperbolicity.requireStrictlyHyperbolic is yes. "
                "Viscous terms make the system incompletely parabolic. "
                "Set flowModel euler, or set requireStrictlyHyperbolic no."
            )
        else:
            report.warnings.append(
                "Navier-Stokes transport is incompletely parabolic. "
                "The inviscid flux is still the reactive Euler flux."
            )
    elif numerics["flowModel"] != "euler":
        fail(f"Unknown numerics.flowModel '{numerics['flowModel']}'")

    if str(numerics["ddtScheme"]) == "steadyState":
        if strict:
            fail("ddtScheme steadyState is not a hyperbolic time evolution")
        else:
            report.warnings.append("ddtScheme steadyState does not advance a time-hyperbolic system")

    cfl = float(numerics["CFL"])
    limit = float(hyper["maxCoStability"])
    if cfl > limit:
        report.warnings.append(
            f"numerics.CFL {cfl} is above hyperbolicity.maxCoStability {limit}. "
            "That is a stability limit for the explicit flux, not a loss of hyperbolicity."
        )

    if boundaries["walls"]["velocityModel"] == "acoustic":
        model = boundaries["acoustics"]["timeDomainModel"]
        if model not in ("slip", "resistive", "helmholtz"):
            fail(f"Acoustic time-domain model '{model}' is not local")

    if reaction_model(config) == "detailed":
        _check_species(config, state, report, fail)
    return report


# Temperatures at which the frozen mixture is checked: the inlet, the hot
# kernel, and a hot burnt-gas value. The burnt composition itself is not
# known before the run, so the check uses the unburned composition.
def _check_species(config: dict, state: Mixture, report: HyperbolicityReport, fail) -> None:
    mechanism = state.mechanism
    fractions = state.species
    if min(fractions.values()) < 0 or abs(sum(fractions.values()) - 1.0) > 1e-9:
        fail("Unburned species mass fractions are not non-negative and summing to 1")
    ignition = config["raw"]["ignition"]
    temperatures = {"inlet": state.temperature}
    if ignition["type"] == "hotKernel":
        temperatures["hot kernel"] = float(ignition["TKernel"])
    for label, temperature in temperatures.items():
        for specie in mechanism.species:
            if fractions[specie.name] > 0 and not specie.in_range(temperature):
                fail(
                    f"{label} temperature {temperature:g} K is outside the janaf range "
                    f"{specie.t_low:g} to {specie.t_high:g} K of {specie.name}"
                )
        properties = mechanism.mixture(fractions, temperature)
        if properties["gamma"] <= float(config["raw"]["hyperbolicity"]["gammaMin"]):
            fail(f"Frozen gamma {properties['gamma']:.6g} at the {label} temperature is not above gammaMin")
        if properties["soundSpeed"] ** 2 <= float(config["raw"]["hyperbolicity"]["c2Min"]):
            fail(f"Frozen c^2 at the {label} temperature is not above c2Min")
    report.warnings.append(
        "detailed: gamma and c^2 were checked for the unburned composition at the inlet "
        "and kernel temperatures. Burnt-gas composition is not known before the run; "
        "the solver checks the initial fields only and does not monitor gamma or c^2 afterwards."
    )


def _yes(value) -> bool:
    return str(value).lower() in ("yes", "true", "on")

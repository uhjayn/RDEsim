"""Normal and oblique reflection from a locally reacting wall.

The microstructure models are lumped acoustic approximations. They are not a
resolved simulation of the cavities. Complex reflection uses the e^{i omega t}
convention, so a closed-pipe impedance is ``-i rho c cot(k L)``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .config import ConfigError


@dataclass
class AcousticResult:
    z_real: float
    z_imag: float
    r_real: float
    r_imag: float
    magnitude: float
    phase_deg: float
    absorption: float
    time_domain_model: str
    inertia: float
    resistance: float
    stiffness: float
    porosity: float
    warning: str


def evaluate_acoustics(config: dict, density: float, sound_speed: float) -> AcousticResult:
    if density <= 0 or sound_speed <= 0:
        raise ConfigError("Acoustic reference density and sound speed must be positive")
    acoustics = config["raw"]["boundaries"]["acoustics"]
    mode = acoustics["mode"]
    angle = math.radians(float(acoustics["incidenceAngleDeg"]))
    cosine = math.cos(angle)
    if cosine <= 0:
        raise ConfigError(
            "boundaries.acoustics.incidenceAngleDeg must be in (-90, 90) "
            "so the wave has a component into the wall"
        )
    if mode == "impedance":
        z_real = float(acoustics["ZReal"])
        z_imag = float(acoustics["ZImag"])
        warning = ""
    elif mode == "reflection":
        reflection = complex(float(acoustics["RReal"]), float(acoustics["RImag"]))
        if abs(1.0 - reflection) < 1e-14:
            raise ConfigError(
                "A reflection coefficient of +1 is a hard wall. Set mode impedance "
                "with a large ZReal instead of inverting that reflection for Z."
            )
        # zeta = (1+R)/(1-R), then Z = zeta * rho c / cos(theta)
        zeta = (1.0 + reflection) / (1.0 - reflection)
        impedance = zeta * density * sound_speed / cosine
        z_real = impedance.real
        z_imag = impedance.imag
        warning = ""
    elif mode == "microstructure":
        z_real, z_imag, warning = _microstructure(config, density, sound_speed)
    else:
        raise ConfigError(f"Unknown boundaries.acoustics.mode '{mode}'")

    reflection = _reflection(z_real, z_imag, density, sound_speed, cosine)
    magnitude = abs(reflection)
    absorption = 1.0 - magnitude ** 2
    phase = math.degrees(math.atan2(reflection.imag, reflection.real))
    coefficients = _time_domain(config, density, sound_speed, z_real)
    if coefficients["warning"]:
        warning = (warning + " " + coefficients["warning"]).strip()
    return AcousticResult(
        z_real=z_real,
        z_imag=z_imag,
        r_real=reflection.real,
        r_imag=reflection.imag,
        magnitude=magnitude,
        phase_deg=phase,
        absorption=absorption,
        time_domain_model=config["raw"]["boundaries"]["acoustics"]["timeDomainModel"],
        inertia=coefficients["inertia"],
        resistance=coefficients["resistance"],
        stiffness=coefficients["stiffness"],
        porosity=coefficients["porosity"],
        warning=warning,
    )


def helmholtz_resonance_frequency(config: dict, sound_speed: float) -> float:
    micro = config["raw"]["boundaries"]["acoustics"]["microstructure"]
    length = _effective_neck(micro)
    area = math.pi * float(micro["neckRadius"]) ** 2
    volume = float(micro["cavityVolume"])
    if length <= 0 or area <= 0 or volume <= 0 or sound_speed <= 0:
        raise ConfigError("Helmholtz geometry and sound speed must be positive")
    return sound_speed / (2.0 * math.pi) * math.sqrt(area / (length * volume))


def _reflection(z_real: float, z_imag: float, density: float, sound_speed: float, cosine: float) -> complex:
    impedance = density * sound_speed
    zeta = complex(z_real, z_imag) / impedance
    return (zeta * cosine - 1.0) / (zeta * cosine + 1.0)


def _microstructure(config: dict, density: float, sound_speed: float):
    acoustics = config["raw"]["boundaries"]["acoustics"]
    micro = acoustics["microstructure"]
    frequency = float(acoustics["frequency"])
    if frequency <= 0:
        raise ConfigError("'boundaries.acoustics.frequency' must be positive")
    wave_number = 2.0 * math.pi * frequency / sound_speed
    kind = micro["type"]
    warning = ""
    if kind == "quarterWave":
        depth = float(micro["cavityDepth"])
        sigma = 1.0
        impedance = _closed_pipe(density, sound_speed, wave_number, depth, float(micro["viscousResistance"]))
        impedance /= sigma
        warning = (
            "A quarter-wave cavity is a distributed duct. The time-domain boundary "
            "uses the selected timeDomainModel and does not reproduce Im(Z) exactly."
        )
    elif kind == "groove":
        width = float(micro["grooveWidth"])
        land = float(micro["grooveLand"])
        if width <= 0 or land < 0 or width + land <= 0:
            raise ConfigError("Groove width must be positive and land must be >= 0")
        sigma = width / (width + land)
        depth = float(micro["cavityDepth"])
        impedance = _closed_pipe(density, sound_speed, wave_number, depth, float(micro["viscousResistance"]))
        impedance /= sigma
        warning = (
            "Groove impedance uses sigma = grooveWidth / (grooveWidth + grooveLand) "
            "and a closed-pipe cavity. The reactive part is reported; the time-domain "
            "boundary follows timeDomainModel."
        )
    elif kind == "helmholtz":
        impedance = _helmholtz(micro, density, sound_speed, wave_number)
        warning = ""
    else:
        raise ConfigError(f"Unknown microstructure type '{kind}'")
    return impedance.real, impedance.imag, warning


def _closed_pipe(density, sound_speed, wave_number, depth, resistance) -> complex:
    if depth <= 0:
        raise ConfigError("Cavity depth must be positive")
    argument = wave_number * depth
    # Guard the cotangent poles. A tiny offset keeps the hard-wall limit finite
    # for reporting; the value is still set by the user's depth and frequency.
    sine = math.sin(argument)
    cosine = math.cos(argument)
    if abs(sine) < 1e-12:
        cotangent = math.copysign(1e12, cosine if cosine != 0 else 1.0)
    else:
        cotangent = cosine / sine
    return complex(resistance, -density * sound_speed * cotangent)


def _helmholtz(micro, density, sound_speed, wave_number) -> complex:
    area = math.pi * float(micro["neckRadius"]) ** 2
    if area <= 0:
        raise ConfigError("'neckRadius' must be positive")
    length = _effective_neck(micro)
    volume = float(micro["cavityVolume"])
    cell = float(micro["cellArea"])
    if volume <= 0 or cell <= 0:
        raise ConfigError("Helmholtz cavityVolume and cellArea must be positive")
    sigma = area / cell
    if sigma <= 0 or sigma > 1:
        raise ConfigError(
            "Helmholtz neck area / cellArea must lie in (0, 1]. "
            f"It is {sigma:.6g}."
        )
    resistance = float(micro["normalizedResistance"])
    specific = density * sound_speed * (
        resistance + 1j * (wave_number * length - area / (wave_number * volume))
    )
    return specific / sigma


def _effective_neck(micro) -> float:
    return float(micro["neckLength"]) + float(micro["endCorrection"]) * float(micro["neckRadius"])


def _time_domain(config, density, sound_speed, z_real) -> dict:
    acoustics = config["raw"]["boundaries"]["acoustics"]
    model = acoustics["timeDomainModel"]
    micro = acoustics["microstructure"]
    warning = ""
    if model == "slip":
        return {"inertia": 0.0, "resistance": 0.0, "stiffness": 0.0, "porosity": 0.0, "warning": warning}
    if model == "resistive":
        if abs(z_real) < float(acoustics["resistiveFloor"]):
            warning = (
                "Resistive |ZReal| is below boundaries.acoustics.resistiveFloor. "
                "A near-zero real impedance is a pressure-release wall, and the "
                "boundary condition refuses to divide by it."
            )
        return {
            "inertia": 0.0,
            "resistance": z_real,
            "stiffness": 0.0,
            "porosity": 1.0,
            "warning": warning,
        }
    if model == "helmholtz":
        area = math.pi * float(micro["neckRadius"]) ** 2
        length = _effective_neck(micro)
        volume = float(micro["cavityVolume"])
        cell = float(micro["cellArea"])
        if length <= 0 or area <= 0 or volume <= 0 or cell <= 0:
            raise ConfigError(
                "timeDomainModel helmholtz needs positive neckLength, neckRadius, "
                "cavityVolume, and cellArea"
            )
        porosity = area / cell
        inertia = density * length
        resistance = density * sound_speed * float(micro["normalizedResistance"])
        stiffness = density * sound_speed ** 2 * area / volume
        return {
            "inertia": inertia,
            "resistance": resistance,
            "stiffness": stiffness,
            "porosity": porosity,
            "warning": warning,
        }
    raise ConfigError(f"Unknown timeDomainModel '{model}'")

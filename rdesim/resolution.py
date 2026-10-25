"""Cell size compared with the one-step chemical length.

The length here is inlet sound speed times the one-step chemical time.
It is not a ZND induction length and not a detonation-cell width.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .config import reaction_model
from .geometry import build_profile
from .mesh import Mesh, build_mesh
from .mixture import mixture_from_config

# Fewer than this many cells across the estimated chemical length, in either
# in-plane direction, is too coarse for detonation cells.
CELLS_ACROSS_CHEMICAL_LENGTH = 10


@dataclass
class ResolutionReport:
    text: str
    coarse: bool
    n_cells: int
    axial_m: float
    circumferential_m: float
    chemical_length_m: float
    cells_axial: float
    cells_circumferential: float


def resolution_report(config: dict) -> ResolutionReport:
    state = mixture_from_config(config)
    profile = build_profile(config)
    mesh = build_mesh(config, profile)
    chemistry = state.chemistry
    edges = _edge_lengths(mesh)
    n_cells = sum(item[3] for item in edges)
    axial = [item[1] for item in edges for _ in range(item[3])]
    circumferential = [item[0] for item in edges for _ in range(item[3])]
    axial_m = _median(axial)
    circumferential_m = _median(circumferential)

    kernel_temperature = float(config["raw"]["ignition"]["TKernel"])
    ignition_temperature = float(chemistry["ignitionTemperature"])
    lines = [
        "Mesh resolution",
        f"  cells                 {n_cells}",
        (
            "  circumferential edge  "
            f"min {_min(circumferential):.6g} m  median {circumferential_m:.6g} m  "
            f"max {_max(circumferential):.6g} m"
        ),
        (
            "  axial edge            "
            f"min {_min(axial):.6g} m  median {axial_m:.6g} m  "
            f"max {_max(axial):.6g} m"
        ),
        (
            "  One-step chemical length uses the inlet sound speed "
            f"{state.sound_speed:.6g} m/s times 1/rate."
        ),
        "  This is not a ZND induction length and not a measured cell width.",
    ]
    if reaction_model(config) != "oneStepProgress":
        lines.append(
            f"  Reaction model is {reaction_model(config)}. The chemical length below "
            "still comes from the one-step chemistry block, not from that model's "
            "induction and release rates, so read it as a mesh-scale reference only."
        )

    kernel_length = _chemical_length(chemistry, state.sound_speed, kernel_temperature)
    ignition_length = _chemical_length(
        chemistry, state.sound_speed, ignition_temperature
    )
    lines.append(
        _length_line("kernel temperature", kernel_temperature, kernel_length, axial_m, circumferential_m)
    )
    lines.append(
        _length_line(
            "ignition temperature",
            ignition_temperature,
            ignition_length,
            axial_m,
            circumferential_m,
        )
    )

    chemical_length = kernel_length if kernel_length and kernel_length > 0 else 0.0
    cells_axial = chemical_length / axial_m if axial_m > 0 and chemical_length > 0 else 0.0
    cells_circ = (
        chemical_length / circumferential_m
        if circumferential_m > 0 and chemical_length > 0
        else 0.0
    )
    coarse = (
        chemical_length <= 0
        or cells_axial < CELLS_ACROSS_CHEMICAL_LENGTH
        or cells_circ < CELLS_ACROSS_CHEMICAL_LENGTH
    )
    if coarse:
        lines.append(
            "  WARNING: the grid is too coarse to resolve detonation cells. "
            f"At the kernel temperature the estimate is {cells_axial:.3g} axial cells "
            f"and {cells_circ:.3g} circumferential cells across the chemical length. "
            f"A cellular structure wants at least {CELLS_ACROSS_CHEMICAL_LENGTH} in both "
            "directions, and a one-step reactive Euler model still tends to make "
            "weak or irregular cells on a fine mesh."
        )
    else:
        lines.append(
            "  The estimate meets "
            f"{CELLS_ACROSS_CHEMICAL_LENGTH} cells across the chemical length "
            "in both in-plane directions. That does not validate a cell size."
        )
    lines.append(
        "  A reactive Euler model with one progress variable does not produce "
        "a validated soot foil. Read pMax as a pressure history, not as a cell width."
    )
    return ResolutionReport(
        text="\n".join(lines),
        coarse=coarse,
        n_cells=n_cells,
        axial_m=axial_m,
        circumferential_m=circumferential_m,
        chemical_length_m=chemical_length,
        cells_axial=cells_axial,
        cells_circumferential=cells_circ,
    )


def _edge_lengths(mesh: Mesh) -> list[tuple[float, float, float, int]]:
    rows = []
    for block in mesh.hexes:
        corners = block.corners
        nx, ny, nz = (max(int(count), 1) for count in block.cells)
        circumferential = _distance(corners[0], corners[1]) / nx
        axial = _distance(corners[0], corners[3]) / ny
        count = nx * ny * nz
        rows.append((circumferential, axial, _distance(corners[0], corners[4]) / nz, count))
    return rows


def _chemical_length(chemistry: dict, sound_speed: float, temperature: float) -> float | None:
    if chemistry.get("model") != "oneStepProgress":
        return None
    if temperature < float(chemistry["ignitionTemperature"]):
        return None
    rate = _rate(chemistry, temperature)
    if rate <= 0 or sound_speed <= 0:
        return None
    return sound_speed / rate


def _rate(chemistry: dict, temperature: float) -> float:
    reference = float(chemistry["TRef"])
    if temperature <= 0 or reference <= 0:
        return 0.0
    theta = temperature / reference
    return (
        float(chemistry["A"])
        * theta ** float(chemistry["beta"])
        * math.exp(-float(chemistry["Ta"]) / temperature)
    )


def _length_line(label, temperature, length, axial, circumferential) -> str:
    if length is None or length <= 0:
        return f"  {label} {temperature:.6g} K: rate is off, so no chemical length is defined."
    return (
        f"  {label} {temperature:.6g} K: chemical length {length:.6g} m, "
        f"{length / axial:.3g} axial cells, "
        f"{length / circumferential:.3g} circumferential cells."
    )


def _distance(a, b) -> float:
    return math.dist(tuple(a), tuple(b))


def _median(values: list[float]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def _min(values: list[float]) -> float:
    return min(values)


def _max(values: list[float]) -> float:
    return max(values)

"""Annulus wall profile: angles, width, steps, nozzle, and injectors."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .config import ConfigError, as_bool


@dataclass
class Segment:
    y0: float
    y1: float
    r_in0: float
    r_out0: float
    r_in1: float
    r_out1: float


@dataclass
class StepFace:
    y: float
    r0: float
    r1: float
    side: str


@dataclass
class Profile:
    segments: list[Segment]
    step_faces: list[StepFace]
    length: float
    recess: float
    unwrap_radius: float
    circumferential_length: float
    slots: list[tuple[float, float]]
    injector_area: float
    radial_steps_resolved: bool
    injector_slots_resolved: bool
    notes: list[str]


def build_profile(config: dict) -> Profile:
    raw = config["raw"]
    domain = raw["domain"]
    geometry = raw["geometry"]
    mesh = raw["mesh"]
    tolerance = float(mesh["stationTolerance"])
    length = float(domain["axialLength"])
    if length <= 0:
        raise ConfigError("'domain.axialLength' must be positive")
    r_inner = float(domain["RInner"])
    r_outer = float(domain["ROuter"])
    if r_inner <= 0 or r_outer <= r_inner:
        raise ConfigError("'domain.ROuter' must be greater than 'domain.RInner', and both positive")

    stations = _stations(domain, geometry, length, tolerance)
    base = {y: _base_radii(y, domain, geometry) for y in stations}
    _apply_nozzle(base, geometry, tolerance)
    segments, step_faces = _apply_steps(base, geometry, tolerance)

    for segment in segments:
        _check_gap(segment)
    recess = float(geometry["injector"]["recessDepth"])
    if recess < 0:
        raise ConfigError("'geometry.injector.recessDepth' must be >= 0")

    unwrap_radius = _unwrap_radius(domain)
    angle = math.radians(float(domain["sectorAngle"]))
    if angle <= 0 or angle > 2.0 * math.pi + 1e-12:
        raise ConfigError("'domain.sectorAngle' must be in (0, 360]")
    circumferential_length = angle * unwrap_radius
    slots = _slots(geometry, circumferential_length)
    area = _injector_area(domain, geometry, slots, circumferential_length, segments)

    dimension = int(domain["dimension"])
    topology = domain["topology"]
    notes = []
    radial_steps = dimension == 3 or topology == "meridional"
    slots_resolved = topology == "unwrapped" and dimension == 2
    if not radial_steps and step_faces:
        notes.append(
            "Radial wall steps are stored in the wall profile. This mesh does "
            "not cut them into the boundary. Set domain.dimension 3, or "
            "domain.dimension 2 with domain.topology meridional, to resolve them."
        )
    if topology == "unwrapped" and dimension == 2:
        notes.append(
            "Unwrapped 2D domain: x is circumferential, y is axial, z is an "
            "empty numerical thickness. Injector slots and the recess are meshed."
        )
    if topology == "meridional":
        notes.append(
            "Meridional 2D domain: x is axial, y is radial, z is empty. "
            "Injector slots are not split; the recess is an annular ring."
        )
    if dimension == 3:
        notes.append(
            "3D domain: z is axial and x-y is the cross-section. Injector slots "
            "are not split circumferentially in this version; the recess is an annular ring."
        )
        if abs(float(geometry["injector"]["injectionAngle"])) > 0:
            raise ConfigError(
                "geometry.injector.injectionAngle must be 0 on a 3D mesh. "
                "A uniform inlet vector cannot carry swirl. Use the unwrapped topology "
                "for a circumferential injection component."
            )

    return Profile(
        segments=segments,
        step_faces=step_faces,
        length=length,
        recess=recess,
        unwrap_radius=unwrap_radius,
        circumferential_length=circumferential_length,
        slots=slots,
        injector_area=area,
        radial_steps_resolved=radial_steps,
        injector_slots_resolved=slots_resolved,
        notes=notes,
    )


def _stations(domain, geometry, length, tolerance) -> list[float]:
    values = {0.0, length}
    nozzle = geometry["nozzle"]
    if as_bool(nozzle["enabled"]):
        start = float(nozzle["start"])
        nozzle_length = float(nozzle["length"])
        end = start + nozzle_length
        if start < 0 or end > length + tolerance:
            raise ConfigError(
                "geometry.nozzle must lie inside the axial length. "
                f"It runs from {start} to {end}, domain length is {length}."
            )
        if nozzle_length <= 0:
            raise ConfigError("'geometry.nozzle.length' must be positive when the nozzle is enabled")
        fraction = float(nozzle["throatFraction"])
        if not 0.0 < fraction < 1.0:
            raise ConfigError("'geometry.nozzle.throatFraction' must lie in (0, 1)")
        values.add(start)
        values.add(start + fraction * nozzle_length)
        values.add(end)
    for entry in geometry["widthProfile"]:
        location = float(entry["location"])
        if location < -tolerance or location > length + tolerance:
            raise ConfigError(f"Width-profile location {location} is outside the axial length")
        values.add(min(max(location, 0.0), length))
    for entry in geometry["steps"]:
        location = float(entry["location"])
        angle = float(entry["angle"])
        if location <= tolerance or location >= length - tolerance:
            raise ConfigError(
                f"Step location {location} must lie strictly inside the axial length"
            )
        if angle <= 0 or angle > 90:
            raise ConfigError(
                f"Step angle {angle} degrees is outside (0, 90]. "
                "90 degrees is a vertical step face."
            )
        values.add(location)
        if angle < 90:
            ramp = abs(float(entry["height"])) / math.tan(math.radians(angle))
            if location + ramp > length + tolerance:
                raise ConfigError(
                    f"Step ramp starting at {location} with angle {angle} "
                    "extends past the axial length. Increase the angle or move the step."
                )
            values.add(location + ramp)
    ordered = sorted(values)
    unique = [ordered[0]]
    for value in ordered[1:]:
        if value - unique[-1] > tolerance:
            unique.append(value)
    return unique


def _base_radii(y: float, domain, geometry) -> tuple[float, float]:
    reference = float(domain["axialReference"])
    dy = y - reference
    r_inner = float(domain["RInner"]) + dy * math.tan(math.radians(float(geometry["innerWallAngle"])))
    r_outer = float(domain["ROuter"]) + dy * math.tan(math.radians(float(geometry["outerWallAngle"])))
    profile = geometry["widthProfile"]
    if profile:
        gap = _interpolate(
            y,
            [float(entry["location"]) for entry in profile],
            [float(entry["gap"]) for entry in profile],
        )
        r_inner, r_outer = _apply_gap(r_inner, r_outer, gap, geometry["widthAnchor"])
    return r_inner, r_outer


def _apply_nozzle(base: dict, geometry, tolerance: float) -> None:
    nozzle = geometry["nozzle"]
    if not as_bool(nozzle["enabled"]):
        return
    start = float(nozzle["start"])
    end = start + float(nozzle["length"])
    if nozzle["shape"] == "piecewise":
        if len(nozzle["profile"]) < 2:
            raise ConfigError(
                "geometry.nozzle.profile needs at least two (location, gap) entries "
                "when geometry.nozzle.shape is piecewise"
            )
        locations = [float(entry["location"]) for entry in nozzle["profile"]]
        gaps = [float(entry["gap"]) for entry in nozzle["profile"]]
        for y in list(base):
            if start - tolerance <= y <= end + tolerance:
                r_inner, r_outer = base[y]
                gap = _interpolate(y, locations, gaps)
                base[y] = _apply_gap(r_inner, r_outer, gap, nozzle["anchor"])
        return
    gap_start = base[_nearest(base, start)][1] - base[_nearest(base, start)][0]
    # Use the pre-nozzle gap exactly at the nozzle start station.
    r_inner, r_outer = base[_nearest(base, start)]
    gap_start = r_outer - r_inner
    throat_y = start + float(nozzle["throatFraction"]) * float(nozzle["length"])
    throat_gap = float(nozzle["throatGap"])
    exit_gap = float(nozzle["exitGap"])
    for y in list(base):
        if y < start - tolerance or y > end + tolerance:
            continue
        if y <= throat_y:
            span = throat_y - start
            fraction = 0.0 if span <= tolerance else (y - start) / span
            gap = gap_start + fraction * (throat_gap - gap_start)
        else:
            span = end - throat_y
            fraction = 1.0 if span <= tolerance else (y - throat_y) / span
            gap = throat_gap + fraction * (exit_gap - throat_gap)
        r_inner, r_outer = base[y]
        base[y] = _apply_gap(r_inner, r_outer, gap, nozzle["anchor"])


def _apply_steps(base: dict, geometry, tolerance: float):
    """Return segments and vertical step faces.

    A 90 degree step shifts the selected wall at every station downstream of
    the step. The station at the step itself keeps the upstream radius; the
    segment that starts there uses the downstream radius. Shallower angles
    spread the same shift over a ramp and do not create a step face.
    """
    stations = sorted(base)
    radii = {y: (base[y][0], base[y][1]) for y in stations}
    for entry in geometry["steps"]:
        if abs(float(entry["height"])) <= tolerance:
            raise ConfigError("Step height must be non-zero. Remove the step or set a height.")
        if float(entry["angle"]) >= 90.0 - 1e-9:
            continue
        location = float(entry["location"])
        height = float(entry["height"])
        ramp = abs(height) / math.tan(math.radians(float(entry["angle"])))
        side = entry["side"]
        for y in list(radii):
            if location + tolerance < y <= location + ramp + tolerance:
                fraction = (y - location) / ramp
                radii[y] = _shift_side(radii[y][0], radii[y][1], side, height * fraction)
            elif y > location + ramp + tolerance:
                radii[y] = _shift_side(radii[y][0], radii[y][1], side, height)

    jump_locations = []
    for entry in geometry["steps"]:
        if float(entry["angle"]) < 90.0 - 1e-9:
            continue
        location = float(entry["location"])
        if location not in jump_locations:
            jump_locations.append(location)
        side = entry["side"]
        height = float(entry["height"])
        for y in list(radii):
            if y > location + tolerance:
                radii[y] = _shift_side(radii[y][0], radii[y][1], side, height)

    ordered = sorted(radii)
    segments = []
    for y0, y1 in zip(ordered, ordered[1:]):
        start = radii[y0]
        end = radii[y1]
        for entry in geometry["steps"]:
            if float(entry["angle"]) < 90.0 - 1e-9:
                continue
            if abs(float(entry["location"]) - y0) <= tolerance:
                start = _shift_side(start[0], start[1], entry["side"], float(entry["height"]))
        segments.append(Segment(y0, y1, start[0], start[1], end[0], end[1]))

    step_faces = []
    for location in jump_locations:
        upstream = radii[location]
        downstream = upstream
        for entry in geometry["steps"]:
            if float(entry["angle"]) < 90.0 - 1e-9:
                continue
            if abs(float(entry["location"]) - location) > tolerance:
                continue
            before = downstream
            downstream = _shift_side(
                downstream[0], downstream[1], entry["side"], float(entry["height"])
            )
            if entry["side"] == "outer":
                r_a, r_b = before[1], downstream[1]
            else:
                r_a, r_b = before[0], downstream[0]
            if abs(r_b - r_a) <= tolerance:
                raise ConfigError(f"Step at {location} does not change the wall radius")
            step_faces.append(StepFace(location, min(r_a, r_b), max(r_a, r_b), entry["side"]))
    return segments, step_faces


def _shift_side(r_inner: float, r_outer: float, side: str, height: float):
    if side == "outer":
        return r_inner, r_outer - height
    return r_inner + height, r_outer


def _check_gap(segment: Segment) -> None:
    for r_inner, r_outer, where in (
        (segment.r_in0, segment.r_out0, segment.y0),
        (segment.r_in1, segment.r_out1, segment.y1),
    ):
        if r_inner <= 0 or r_outer <= r_inner:
            raise ConfigError(
                f"Wall profile at axial coordinate {where} has inner radius "
                f"{r_inner} and outer radius {r_outer}. The gap must stay positive."
            )


def _unwrap_radius(domain) -> float:
    choice = domain["unwrapRadius"]
    if choice == "inner":
        return float(domain["RInner"])
    if choice == "outer":
        return float(domain["ROuter"])
    return 0.5 * (float(domain["RInner"]) + float(domain["ROuter"]))


def _slots(geometry, circumferential_length: float) -> list[tuple[float, float]]:
    injector = geometry["injector"]
    count = int(injector["nSlots"])
    width = float(injector["slotWidth"])
    phase = float(injector["phase"])
    if count < 0:
        raise ConfigError("'geometry.injector.nSlots' must be >= 0")
    if count == 0:
        return [(0.0, circumferential_length)]
    if width <= 0:
        raise ConfigError("'geometry.injector.slotWidth' must be positive when nSlots > 0")
    spacing = circumferential_length / count
    if width >= spacing:
        raise ConfigError(
            "Injector slots overlap. Reduce geometry.injector.slotWidth or "
            "geometry.injector.nSlots, or increase domain.sectorAngle."
        )
    slots = []
    for index in range(count):
        centre = phase + (index + 0.5) * spacing
        start = centre - 0.5 * width
        end = centre + 0.5 * width
        if start < 0 or end > circumferential_length:
            raise ConfigError(
                "An injector slot crosses the periodic boundary. "
                "Change geometry.injector.phase so every slot lies inside the sector."
            )
        slots.append((start, end))
    return slots


def _injector_area(domain, geometry, slots, circumferential_length, segments: list[Segment]) -> float:
    inlet = segments[0]
    r_inner = inlet.r_in0
    r_outer = inlet.r_out0
    sector = math.radians(float(domain["sectorAngle"]))
    full = 0.5 * sector * (r_outer**2 - r_inner**2)
    if int(geometry["injector"]["nSlots"]) == 0:
        return full
    coverage = sum(end - start for start, end in slots) / circumferential_length
    return full * coverage


def _apply_gap(r_inner: float, r_outer: float, gap: float, anchor: str):
    if gap <= 0:
        raise ConfigError(f"A requested gap of {gap} is not positive")
    if anchor == "outer":
        return r_outer - gap, r_outer
    if anchor == "inner":
        return r_inner, r_inner + gap
    mid = 0.5 * (r_inner + r_outer)
    return mid - 0.5 * gap, mid + 0.5 * gap


def _interpolate(y: float, locations: list[float], values: list[float]) -> float:
    pairs = sorted(zip(locations, values))
    if y <= pairs[0][0]:
        return pairs[0][1]
    if y >= pairs[-1][0]:
        return pairs[-1][1]
    for (y0, v0), (y1, v1) in zip(pairs, pairs[1:]):
        if y0 <= y <= y1:
            if y1 == y0:
                return v0
            fraction = (y - y0) / (y1 - y0)
            return v0 + fraction * (v1 - v0)
    return pairs[-1][1]


def _nearest(mapping: dict, y: float) -> float:
    return min(mapping, key=lambda station: abs(station - y))


def _nearest_below(mapping: dict, y: float) -> float:
    below = [station for station in mapping if station <= y + 1e-12]
    if not below:
        return min(mapping)
    return max(below)

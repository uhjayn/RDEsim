"""blockMeshDict generation for the unwrapped, meridional, and 3D annuli."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .config import ConfigError
from .geometry import Profile, Segment, StepFace


@dataclass
class Hex:
    corners: list
    cells: tuple
    grading: tuple
    faces: dict


@dataclass
class Mesh:
    hexes: list
    topology: str
    dimension: int
    notes: list[str] = field(default_factory=list)

    def block_mesh_dict(self) -> str:
        points = []
        index = {}

        def pid(point):
            key = tuple(round(float(coord), 8) for coord in point)
            if key not in index:
                index[key] = len(points)
                points.append(key)
            return index[key]

        indexed = []
        arcs = {}
        for hex_block in self.hexes:
            ids = [pid(corner) for corner in hex_block.corners]
            indexed.append((ids, hex_block))
            self._collect_arcs(hex_block, ids, arcs)

        lines = ["vertices", "("]
        for point in points:
            lines.append(
                "    ("
                + " ".join(_num(coord) for coord in point)
                + ")"
            )
        lines.append(");")
        lines.append("")
        if arcs and self.dimension == 3:
            lines.append("edges")
            lines.append("(")
            for (a, b), point in sorted(arcs.items()):
                lines.append(
                    "    arc "
                    + str(a)
                    + " "
                    + str(b)
                    + " ("
                    + " ".join(_num(coord) for coord in point)
                    + ")"
                )
            lines.append(");")
            lines.append("")
        lines.append("blocks")
        lines.append("(")
        for ids, hex_block in indexed:
            cells = " ".join(str(int(n)) for n in hex_block.cells)
            grading = " ".join(_num(g) for g in hex_block.grading)
            lines.append(
                "    hex ("
                + " ".join(str(i) for i in ids)
                + ") ("
                + cells
                + ") simpleGrading ("
                + grading
                + ")"
            )
        lines.append(");")
        lines.append("")
        lines.append("defaultPatch")
        lines.append("{")
        lines.append("    name unassigned;")
        lines.append("    type patch;")
        lines.append("}")
        lines.append("")
        patches = {}
        for ids, hex_block in indexed:
            for name, local_faces in hex_block.faces.items():
                patches.setdefault(name, [])
                for local in local_faces:
                    patches[name].append(tuple(ids[vertex] for vertex in local))
        lines.append("boundary")
        lines.append("(")
        for name in _patch_order(patches):
            kind = _patch_type(name)
            lines.append(f"    {name}")
            lines.append("    {")
            lines.append(f"        type {kind};")
            if kind == "cyclic":
                lines.append(f"        neighbourPatch {_cyclic_neighbour(name)};")
            lines.append("        faces")
            lines.append("        (")
            for face in patches[name]:
                lines.append("            (" + " ".join(str(i) for i in face) + ")")
            lines.append("        );")
            lines.append("    }")
        lines.append(");")
        lines.append("")
        return "\n".join(lines)

    def _collect_arcs(self, hex_block: Hex, ids: list[int], arcs: dict) -> None:
        if self.dimension != 3:
            return
        # Theta edges in the vertex order used by ``_sweep``.
        for start, end in ((0, 1), (3, 2), (4, 5), (7, 6)):
            a = hex_block.corners[start]
            b = hex_block.corners[end]
            key = tuple(sorted((ids[start], ids[end])))
            midpoint = _arc_point(a, b)
            arcs.setdefault(key, midpoint)


def build_mesh(config: dict, profile: Profile) -> Mesh:
    raw = config["raw"]
    domain = raw["domain"]
    dimension = int(domain["dimension"])
    topology = domain["topology"]
    _check_resolution(raw, dimension, topology)
    if dimension == 2 and topology == "unwrapped":
        hexes = _unwrapped(raw, profile)
    elif dimension == 2 and topology == "meridional":
        hexes = _meridional(raw, profile)
    elif dimension == 3:
        hexes = _annulus_3d(raw, profile)
    else:
        raise ConfigError(
            f"No mesh layout for dimension {dimension} and topology {topology}"
        )
    if not hexes:
        raise ConfigError("The mesh description produced no blocks")
    return Mesh(hexes=hexes, topology=topology, dimension=dimension, notes=list(profile.notes))


def _check_resolution(raw, dimension: int, topology: str) -> None:
    mesh = raw["mesh"]
    n_radial = int(mesh["nRadial"])
    n_circ = int(mesh["nCircumferential"])
    if dimension == 2 and topology == "unwrapped" and n_radial != 1:
        raise ConfigError(
            "mesh.nRadial must be 1 for a 2D unwrapped mesh. "
            "That direction is the empty numerical thickness."
        )
    if dimension == 2 and topology == "meridional" and n_circ != 1:
        raise ConfigError(
            "mesh.nCircumferential must be 1 for a meridional mesh. "
            "That direction is the empty numerical thickness."
        )
    if dimension == 3:
        blocks = int(mesh["nCircumferentialBlocks"])
        sector = float(raw["domain"]["sectorAngle"])
        if blocks < 1:
            raise ConfigError("'mesh.nCircumferentialBlocks' must be at least 1")
        if sector / blocks >= 180.0:
            raise ConfigError(
                "Each circumferential block must span less than 180 degrees. "
                "Increase mesh.nCircumferentialBlocks."
            )


def _unwrapped(raw, profile: Profile) -> list[Hex]:
    mesh = raw["mesh"]
    thickness = float(raw["domain"]["emptyThickness"])
    if thickness <= 0:
        raise ConfigError("'domain.emptyThickness' must be positive")
    length = profile.circumferential_length
    intervals = _x_intervals(profile.slots, length)
    y_blocks = []
    if profile.recess > 0:
        y_blocks.append((-profile.recess, 0.0, "recess"))
    y_blocks.append((0.0, profile.length, "chamber"))
    x_cells = _distribute(
        [end - start for start, end, _kind in intervals],
        int(mesh["nCircumferential"]),
        int(mesh["minCellsPerBlock"]),
        "circumferential",
    )
    y_cells = _distribute(
        [end - start for start, end, _kind in y_blocks],
        int(mesh["nAxial"]),
        int(mesh["minCellsPerBlock"]),
        "axial",
    )
    grading = (
        float(mesh["circumferentialGrading"]),
        float(mesh["axialGrading"]),
        float(mesh["radialGrading"]),
    )
    hexes = []
    for (y0, y1, y_kind), ny in zip(y_blocks, y_cells):
        for (x0, x1, x_kind), nx in zip(intervals, x_cells):
            if y_kind == "recess" and x_kind != "slot":
                continue
            faces = {
                "front": [(0, 3, 2, 1)],
                "back": [(4, 5, 6, 7)],
            }
            if abs(x0) <= 1e-12:
                faces["periodicLeft"] = [(0, 4, 7, 3)]
            if abs(x1 - length) <= 1e-10:
                faces["periodicRight"] = [(1, 5, 6, 2)]
            if y_kind == "recess":
                faces["injector"] = [(0, 1, 5, 4)]
                faces["injectorSide"] = []
                # Side walls of a slot that is not on the periodic boundary.
                if abs(x0) > 1e-12:
                    faces["injectorSide"].append((0, 4, 7, 3))
                if abs(x1 - length) > 1e-10:
                    faces["injectorSide"].append((1, 5, 6, 2))
                if not faces["injectorSide"]:
                    del faces["injectorSide"]
            else:
                if abs(y1 - profile.length) <= 1e-12:
                    faces["outlet"] = [(3, 7, 6, 2)]
                if abs(y0) <= 1e-12 and profile.recess == 0 and x_kind == "slot":
                    faces["injector"] = [(0, 1, 5, 4)]
                if abs(y0) <= 1e-12 and x_kind == "land":
                    faces["inletWall"] = [(0, 1, 5, 4)]
                if abs(y0) <= 1e-12 and profile.recess > 0 and x_kind == "land":
                    faces["inletWall"] = [(0, 1, 5, 4)]
            corners = [
                (x0, y0, 0.0),
                (x1, y0, 0.0),
                (x1, y1, 0.0),
                (x0, y1, 0.0),
                (x0, y0, thickness),
                (x1, y0, thickness),
                (x1, y1, thickness),
                (x0, y1, thickness),
            ]
            hexes.append(Hex(corners, (nx, ny, 1), grading, faces))
    return hexes


def _x_intervals(slots, length):
    edges = {0.0, length}
    for start, end in slots:
        edges.add(start)
        edges.add(end)
    ordered = sorted(edges)
    intervals = []
    slots = list(slots)
    for left, right in zip(ordered, ordered[1:]):
        mid = 0.5 * (left + right)
        kind = "slot" if any(start - 1e-12 <= mid <= end + 1e-12 for start, end in slots) else "land"
        intervals.append((left, right, kind))
    return intervals


def _meridional(raw, profile: Profile) -> list[Hex]:
    annular = _annular_blocks(profile, raw)
    thickness = float(raw["domain"]["emptyThickness"])
    mesh = raw["mesh"]
    grading = (
        float(mesh["axialGrading"]),
        float(mesh["radialGrading"]),
        1.0,
    )
    hexes = []
    for block in annular:
        corners = [
            (block.y0, block.r_in0, 0.0),
            (block.y1, block.r_in1, 0.0),
            (block.y1, block.r_out1, 0.0),
            (block.y0, block.r_out0, 0.0),
            (block.y0, block.r_in0, thickness),
            (block.y1, block.r_in1, thickness),
            (block.y1, block.r_out1, thickness),
            (block.y0, block.r_out0, thickness),
        ]
        # Edges: 0->1 axial, 1->2 radial at the downstream end, 0->4 empty.
        # Radial cell count belongs on edge 1->2, which is not a uniform
        # radial edge if the downstream radius differs from the upstream one.
        # Use the same count the user requested for this radial piece.
        faces = {"front": [(0, 3, 2, 1)], "back": [(4, 5, 6, 7)]}
        _map_annular_faces(block, faces, axial_is_x=True)
        hexes.append(Hex(corners, (block.n_axial, block.n_radial, 1), grading, faces))
    return hexes


def _annulus_3d(raw, profile: Profile) -> list[Hex]:
    annular = _annular_blocks(profile, raw)
    mesh = raw["mesh"]
    domain = raw["domain"]
    n_blocks = int(mesh["nCircumferentialBlocks"])
    theta0 = math.radians(float(domain["theta0"]))
    sector = math.radians(float(domain["sectorAngle"]))
    closed = abs(float(domain["sectorAngle"]) - 360.0) <= 1e-8
    circ_cells = _distribute(
        [1.0] * n_blocks,
        int(mesh["nCircumferential"]),
        int(mesh["minCellsPerBlock"]),
        "circumferential",
    )
    grading = (
        float(mesh["circumferentialGrading"]),
        float(mesh["radialGrading"]),
        float(mesh["axialGrading"]),
    )
    hexes = []
    for block in annular:
        for index in range(n_blocks):
            a0 = theta0 + sector * index / n_blocks
            a1 = theta0 + sector * (index + 1) / n_blocks
            corners = _ring_corners(block, a0, a1)
            faces = {}
            _map_annular_faces(block, faces, axial_is_x=False)
            if not closed:
                if index == 0:
                    faces["periodicLeft"] = [(0, 3, 7, 4)]
                if index == n_blocks - 1:
                    faces["periodicRight"] = [(1, 5, 6, 2)]
            hexes.append(
                Hex(corners, (circ_cells[index], block.n_radial, block.n_axial), grading, faces)
            )
    return hexes


def _ring_corners(block, theta0, theta1):
    def point(axial, radius, theta):
        return (radius * math.cos(theta), radius * math.sin(theta), axial)

    return [
        point(block.y0, block.r_in0, theta0),
        point(block.y0, block.r_in0, theta1),
        point(block.y0, block.r_out0, theta1),
        point(block.y0, block.r_out0, theta0),
        point(block.y1, block.r_in1, theta0),
        point(block.y1, block.r_in1, theta1),
        point(block.y1, block.r_out1, theta1),
        point(block.y1, block.r_out1, theta0),
    ]


def _map_annular_faces(block, faces, axial_is_x: bool) -> None:
    # Meridional vertex order uses axial as edge 0->1.
    # 3D vertex order uses theta as edge 0->1 and axial as edge 0->4.
    if axial_is_x:
        inlet = (0, 4, 7, 3)
        outlet = (1, 2, 6, 5)
        inner = (0, 1, 5, 4)
        outer = (3, 7, 6, 2)
    else:
        inlet = (0, 3, 2, 1)
        outlet = (4, 5, 6, 7)
        inner = (0, 4, 5, 1)
        outer = (3, 2, 6, 7)
    for name, local in (
        (block.patch_y0, inlet),
        (block.patch_y1, outlet),
        (block.inner, inner),
        (block.outer, outer),
    ):
        if name:
            faces.setdefault(name, []).append(local)


def _annular_blocks(profile: Profile, raw) -> list:
    mesh = raw["mesh"]
    minimum = int(mesh["minCellsPerBlock"])
    n_axial = int(mesh["nAxial"])
    n_radial = int(mesh["nRadial"])
    pieces = _split_radial(profile.segments, profile.step_faces, float(mesh["stationTolerance"]))
    if profile.recess > 0:
        inlet = pieces[0]
        pieces = [
            _Piece(
                -profile.recess,
                0.0,
                inlet.r_in0,
                inlet.r_out0,
                inlet.r_in0,
                inlet.r_out0,
                y0p="injector",
                y1p="",
                inner="innerWall",
                outer="outerWall",
            )
        ] + pieces
        pieces[1].y0p = ""
    spans = []
    span_pieces: dict[tuple, list] = {}
    for piece in pieces:
        key = (round(piece.y0, 8), round(piece.y1, 8))
        if key not in span_pieces:
            spans.append(key)
            span_pieces[key] = []
        span_pieces[key].append(piece)
    axial_counts = _distribute(
        [key[1] - key[0] for key in spans],
        n_axial,
        minimum,
        "axial",
    )
    blocks = []
    for key, axial_count in zip(spans, axial_counts):
        siblings = span_pieces[key]
        if len(siblings) == 1:
            radial_counts = [n_radial]
        else:
            radial_counts = _distribute(
                [max(item.r_out0 - item.r_in0, 0.0) for item in siblings],
                n_radial,
                minimum,
                "radial",
            )
        for piece, radial_count in zip(siblings, radial_counts):
            blocks.append(_view(piece, axial_count, radial_count))
    return blocks


@dataclass
class _Piece:
    y0: float
    y1: float
    r_in0: float
    r_out0: float
    r_in1: float
    r_out1: float
    y0p: str = ""
    y1p: str = ""
    inner: str = ""
    outer: str = ""

    # Aliases used by the face mapper. The field names ``y0`` and ``y1`` are
    # coordinates, so patch names live in y0p/y1p and are exposed below.


class _View:
    def __init__(self, piece: _Piece, n_axial: int, n_radial: int):
        self.y0 = piece.y0
        self.y1 = piece.y1
        self.r_in0 = piece.r_in0
        self.r_out0 = piece.r_out0
        self.r_in1 = piece.r_in1
        self.r_out1 = piece.r_out1
        self.n_axial = n_axial
        self.n_radial = n_radial
        self.y0 = piece.y0
        # Patch names cannot share attribute names with coordinates.
        self.patch_y0 = piece.y0p
        self.patch_y1 = piece.y1p
        self.inner = piece.inner
        self.outer = piece.outer


def _view(piece, n_axial, n_radial):
    view = _View(piece, n_axial, n_radial)
    return view


def _split_radial(segments: list[Segment], step_faces: list[StepFace], tolerance: float) -> list[_Piece]:
    pieces = []
    for index, segment in enumerate(segments):
        piece = _Piece(
            segment.y0,
            segment.y1,
            segment.r_in0,
            segment.r_out0,
            segment.r_in1,
            segment.r_out1,
            y0p="injector" if index == 0 else "",
            y1p="outlet" if index == len(segments) - 1 else "",
            inner="innerWall",
            outer="outerWall",
        )
        pieces.extend(_cut_piece(piece, segments, index, step_faces, tolerance))
    return pieces


def _cut_piece(piece, segments, index, step_faces, tolerance) -> list[_Piece]:
    if index + 1 >= len(segments):
        return [piece]
    nxt = segments[index + 1]
    cuts = []
    if abs(piece.r_out1 - nxt.r_out0) > tolerance and piece.r_out1 > nxt.r_out0:
        cuts.append(("outer", nxt.r_out0))
    if abs(piece.r_in1 - nxt.r_in0) > tolerance and piece.r_in1 < nxt.r_in0:
        cuts.append(("inner", nxt.r_in0))
    # Expansion (the downstream radius sticks out) is cut from the downstream piece.
    if abs(piece.r_out1 - nxt.r_out0) > tolerance and piece.r_out1 < nxt.r_out0:
        # Handled when the downstream piece is created: it sees the previous end.
        pass
    previous = segments[index - 1] if index else None
    if previous is not None and abs(previous.r_out1 - piece.r_out0) > tolerance and previous.r_out1 < piece.r_out0:
        cuts.append(("outer-inlet", previous.r_out1))
    if previous is not None and abs(previous.r_in1 - piece.r_in0) > tolerance and previous.r_in1 > piece.r_in0:
        cuts.append(("inner-inlet", previous.r_in1))
    if not cuts:
        return [piece]
    result = [piece]
    for kind, radius in cuts:
        result = _apply_cut(result, kind, radius, tolerance)
    return result


def _apply_cut(pieces: list[_Piece], kind: str, radius: float, tolerance: float) -> list[_Piece]:
    output = []
    for piece in pieces:
        if kind == "outer":
            if not (
                piece.r_in0 + tolerance < radius < piece.r_out0 - tolerance
                and piece.r_in1 + tolerance < radius < piece.r_out1 - tolerance
            ):
                raise ConfigError(
                    "Cannot cut a radial step out of the wall segment ending at "
                    f"{piece.y1}. Increase the gap or reduce the step height so the "
                    "downstream radius stays inside the upstream wall."
                )
            inner = _Piece(
                piece.y0, piece.y1, piece.r_in0, radius, piece.r_in1, radius,
                piece.y0p, piece.y1p, piece.inner, "",
            )
            outer = _Piece(
                piece.y0, piece.y1, radius, piece.r_out0, radius, piece.r_out1,
                "", "outerStep", "", piece.outer,
            )
            output.extend((inner, outer))
        elif kind == "inner":
            if not (
                piece.r_in0 + tolerance < radius < piece.r_out0 - tolerance
                and piece.r_in1 + tolerance < radius < piece.r_out1 - tolerance
            ):
                raise ConfigError(
                    "Cannot cut an inner-wall step out of the segment ending at "
                    f"{piece.y1}."
                )
            inner = _Piece(
                piece.y0, piece.y1, piece.r_in0, radius, piece.r_in1, radius,
                "", "innerStep", piece.inner, "",
            )
            outer = _Piece(
                piece.y0, piece.y1, radius, piece.r_out0, radius, piece.r_out1,
                piece.y0p, piece.y1p, "", piece.outer,
            )
            output.extend((inner, outer))
        elif kind == "outer-inlet":
            if not (
                piece.r_in0 + tolerance < radius < piece.r_out0 - tolerance
                and piece.r_in1 + tolerance < radius < piece.r_out1 - tolerance
            ):
                raise ConfigError(
                    "Cannot cut an expansion step at the inlet of the segment "
                    f"starting at {piece.y0}."
                )
            inner = _Piece(
                piece.y0, piece.y1, piece.r_in0, radius, piece.r_in1, radius,
                piece.y0p, piece.y1p, piece.inner, "",
            )
            outer = _Piece(
                piece.y0, piece.y1, radius, piece.r_out0, radius, piece.r_out1,
                "outerStep", "", "", piece.outer,
            )
            output.extend((inner, outer))
        elif kind == "inner-inlet":
            if not (
                piece.r_in0 + tolerance < radius < piece.r_out0 - tolerance
                and piece.r_in1 + tolerance < radius < piece.r_out1 - tolerance
            ):
                raise ConfigError(
                    "Cannot cut an inner expansion at the inlet of the segment "
                    f"starting at {piece.y0}."
                )
            inner = _Piece(
                piece.y0, piece.y1, piece.r_in0, radius, piece.r_in1, radius,
                "innerStep", "", piece.inner, "",
            )
            outer = _Piece(
                piece.y0, piece.y1, radius, piece.r_out0, radius, piece.r_out1,
                piece.y0p, piece.y1p, "", piece.outer,
            )
            output.extend((inner, outer))
        else:
            output.append(piece)
    return output


def _distribute(lengths, total: int, minimum: int, label: str) -> list[int]:
    count = len(lengths)
    if count == 0:
        return []
    if minimum < 1:
        raise ConfigError("'mesh.minCellsPerBlock' must be at least 1")
    if total < count * minimum:
        raise ConfigError(
            f"mesh resolution in the {label} direction is {total} cells for "
            f"{count} blocks, below the minimum of {minimum} cells per block. "
            "Increase the cell count or lower mesh.minCellsPerBlock."
        )
    extra = total - count * minimum
    span = sum(lengths)
    if span <= 0:
        raise ConfigError(f"A {label} mesh interval has zero length")
    ideal = [extra * length / span for length in lengths]
    base = [int(math.floor(value)) for value in ideal]
    remainder = extra - sum(base)
    order = sorted(range(count), key=lambda i: (ideal[i] - base[i], lengths[i]), reverse=True)
    for index in range(remainder):
        base[order[index]] += 1
    return [minimum + value for value in base]


def _patch_type(name: str) -> str:
    if name in ("periodicLeft", "periodicRight"):
        return "cyclic"
    if name in ("front", "back"):
        return "empty"
    if name in ("innerWall", "outerWall", "inletWall", "injectorSide", "outerStep", "innerStep"):
        return "wall"
    return "patch"


def _cyclic_neighbour(name: str) -> str:
    if name == "periodicLeft":
        return "periodicRight"
    if name == "periodicRight":
        return "periodicLeft"
    raise ConfigError(f"Patch '{name}' is not cyclic")


def _patch_order(patches: dict) -> list[str]:
    preferred = [
        "injector",
        "inletWall",
        "injectorSide",
        "outlet",
        "innerWall",
        "outerWall",
        "innerStep",
        "outerStep",
        "periodicLeft",
        "periodicRight",
        "front",
        "back",
    ]
    names = list(patches)
    return [name for name in preferred if name in patches] + [
        name for name in names if name not in preferred
    ]


def _arc_point(a, b):
    # Mid-angle point at the average radius and axial coordinate.
    radius_a = math.hypot(a[0], a[1])
    radius_b = math.hypot(b[0], b[1])
    theta_a = math.atan2(a[1], a[0])
    theta_b = math.atan2(b[1], b[0])
    if theta_b - theta_a > math.pi:
        theta_b -= 2.0 * math.pi
    if theta_a - theta_b > math.pi:
        theta_a -= 2.0 * math.pi
    theta = 0.5 * (theta_a + theta_b)
    radius = 0.5 * (radius_a + radius_b)
    axial = 0.5 * (a[2] + b[2])
    return (radius * math.cos(theta), radius * math.sin(theta), axial)


def _num(value) -> str:
    text = format(float(value), ".10g")
    if text in ("-0", "-0.0"):
        return "0"
    return text

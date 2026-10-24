"""Legacy ASCII VTK for the block outline.

ParaView opens this file directly. It is one hexahedron per blockMesh block,
not the refined mesh ``blockMesh`` writes.
"""

from __future__ import annotations

from pathlib import Path


def write_blocks(path: str | Path, hexes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    points = []
    cells = []
    for hex_block in hexes:
        start = len(points)
        points.extend(hex_block.corners)
        cells.append(tuple(start + index for index in range(8)))
    lines = [
        "# vtk DataFile Version 3.0",
        "RDEsim block outline",
        "ASCII",
        "DATASET UNSTRUCTURED_GRID",
        f"POINTS {len(points)} float",
    ]
    for point in points:
        lines.append(" ".join(f"{float(coord):.8g}" for coord in point))
    size = len(cells) * 9
    lines.append(f"CELLS {len(cells)} {size}")
    for cell in cells:
        lines.append("8 " + " ".join(str(index) for index in cell))
    lines.append(f"CELL_TYPES {len(cells)}")
    lines.append("\n".join(["12"] * len(cells)))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

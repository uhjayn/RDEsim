#!/usr/bin/env python3
"""Render the running pressure maximum and the mesh edges.

    pvbatch --force-offscreen-rendering paraview/cell_structure.py <case>

pMax is the per-cell running maximum of pressure. A grey image of that field
is only a pressure history on this mesh. It is not a measured detonation-cell
size. The one-step reactive Euler model tends to make weak or irregular cells,
and this annulus is too coarse for them.
"""

import sys
from pathlib import Path

from paraview.simple import (
    ColorBy,
    GetActiveViewOrCreate,
    GetColorTransferFunction,
    OpenFOAMReader,
    Render,
    ResetCamera,
    SaveScreenshot,
    Show,
)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: pvbatch paraview/cell_structure.py <case-directory>")
    case = Path(sys.argv[1]).resolve()
    foam = case / f"{case.name}.foam"
    if not foam.is_file():
        raise SystemExit(f"missing OpenFOAM reader stub: {foam}")

    out = case / "postProcessing" / "cellStructure"
    out.mkdir(parents=True, exist_ok=True)

    reader = OpenFOAMReader(FileName=str(foam))
    reader.CaseType = "Reconstructed Case"
    reader.MeshRegions = ["internalMesh"]
    available = list(reader.CellArrays.Available)
    wanted = [name for name in ("pMax", "tAtPMax", "p", "U") if name in available]
    if not wanted:
        raise SystemExit(
            "The reader has no cell fields yet. Run the case before the cell-structure view."
        )
    reader.CellArrays = wanted
    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [1200, 700]
    display = Show(reader, view)
    display.Representation = "Surface With Edges"

    print("This image is a pressure history, not a validated cell size.")
    print("Cell arrays:", ", ".join(wanted))
    if "pMax" in wanted:
        ColorBy(display, ("CELLS", "pMax"))
        lut = GetColorTransferFunction("pMax")
        lut.ApplyPreset("Grayscale", True)
        display.SetScalarBarVisibility(view, True)
        ResetCamera(view)
        Render()
        image = out / "pMax.png"
        SaveScreenshot(str(image), view)
        print(f"Wrote {image}")
        print("Grey pMax is high pressure in white. Do not read a cell width from it.")
    else:
        print("pMax was not written. Run the case so the peakPressure function object can store it.")
        ColorBy(display, ("CELLS", wanted[0]))

    display.Representation = "Wireframe"
    ResetCamera(view)
    Render()
    edges = out / "mesh-edges.png"
    SaveScreenshot(str(edges), view)
    print(f"Wrote {edges}")
    print("The wireframe is the mesh. Compare its spacing with constant/meshResolution.")


if __name__ == "__main__":
    main()

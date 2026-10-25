#!/usr/bin/env python3
"""Open an RDEsim case with ParaView's OpenFOAM reader.

From the repository, after ParaView is on PATH (the openfoam14 package
installs one when it is missing):

    pvbatch paraview/view_rde.py cases/annulus2D

The script colours pressure. In the ParaView UI, change the colour array to
T, rho, U, or lambda. Internal energy e is there after the solver has written
it. There is no separate detonation field: the wave is the structure in
pressure, temperature, density, and lambda.

VTK/geometryBlocks.vtk is only the block outline, not the refined mesh.
Use the .foam file for the mesh and the fields.

This script was not executed in the environment that wrote it. A .pvsm state
file is not shipped, because those files store absolute paths.
"""

import sys
from pathlib import Path

from paraview.simple import (
    ColorBy,
    GetActiveViewOrCreate,
    OpenFOAMReader,
    Render,
    Show,
)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: pvbatch paraview/view_rde.py <case-directory>")
    case = Path(sys.argv[1]).resolve()
    foam = case / f"{case.name}.foam"
    if not foam.is_file():
        raise SystemExit(f"missing OpenFOAM reader stub: {foam}")

    reader = OpenFOAMReader(FileName=str(foam))
    reader.CaseType = "Reconstructed Case"
    reader.MeshRegions = ["internalMesh"]
    wanted = ["p", "T", "U", "rho", "lambda", "e"]
    available = list(reader.CellArrays.Available)
    reader.CellArrays = [name for name in wanted if name in available]
    view = GetActiveViewOrCreate("RenderView")
    display = Show(reader, view)
    ColorBy(display, ("CELLS", "p"))
    view.ResetCamera()
    Render()
    print("Coloured by p. Also available:", ", ".join(reader.CellArrays))
    print("Switch the colour array to T, rho, U, or lambda to see the wave.")


if __name__ == "__main__":
    main()

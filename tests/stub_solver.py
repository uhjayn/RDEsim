"""Stand-in solver for the optimization smoke test.

Reads the step height written into constant/wallProfile and writes the
objective file rdeFluid would write after a real run. The value grows with
step height so the driver has a geometry dependence to optimise.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdesim.foamdict import parse_file


def main() -> None:
    case = Path.cwd()
    profile = parse_file(case / "constant" / "wallProfile")
    height = 0.0
    for face in profile.get("stepFaces", []):
        # stepFaces are dictionaries with r0 and r1. The configured height is
        # the change in radius.
        height = abs(float(face["r1"]) - float(face["r0"]))
    destination = case / "postProcessing" / "rdeObjectives"
    destination.mkdir(parents=True, exist_ok=True)
    text = (
        "# stub objectives\n"
        "pMeanInlet 100000\n"
        f"pMeanOutlet {100000.0 + height * 1.0e7}\n"
        "TMax 300\n"
        "lambdaMax 0\n"
    )
    (destination / "objectives.dat").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()

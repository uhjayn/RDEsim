#!/usr/bin/env python3
"""Rebuild the OpenFOAM mechanism files from the published Burke 2012 data.

Input:  source/Burke2012.yaml (Cantera YAML converted from the Burke et al.
        mechanism, Caltech Explosion Dynamics Laboratory SDToolbox copy).
Output: thermo and reactions (OpenFOAM format, via chemkinToFoam). The
        intermediate CHEMKIN-II files are kept under build/.

Only the H2-air subset is kept: species H H2 O OH H2O O2 HO2 H2O2 N2. The
argon, helium, CO and CO2 species are dropped together with the four
reactions that exist only for the argon and helium baths (Burke reactions 7,
8, 10, 11). Every kept rate constant, third-body efficiency, falloff
parameter and NASA-7 coefficient is copied from the source file unchanged.

Needs PyYAML and an OpenFOAM 14 shell (chemkinToFoam). This is a developer
step: the generated `thermo` and `reactions` are committed, and nothing in
rdesim imports this file.
"""

import re
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
SPECIES = ["H", "H2", "O", "OH", "H2O", "O2", "HO2", "H2O2", "N2"]
DROPPED = {"AR", "HE", "CO", "CO2"}

# N2 is published for 300-5000 K, the other eight species for 200 K upward.
# The OpenFOAM janaf evaluation warns on every call below the lower limit, and
# a cold-gas cell at the 300 K ambient temperature (the default) rounds just
# under it, which floods the log. N2's lower limit is therefore set to 200 K
# like the others. Its coefficients are unchanged, so between 200 and 300 K it
# is an extrapolation of the published low-temperature polynomial: cp is about
# 1.1 % below tabulated N2 at 200 K and 0.1 % at 290 K. No gas this cold is
# part of a detonation; it only keeps rarefied fresh gas from tripping the check.
LOW_LIMIT_OVERRIDE = {"N2": 200.0}

# Standard CHEMKIN-II NASA-7 record: 18 name, 6 date, 20 formula, 1 phase,
# 10 Tlow, 10 Thigh, 8 Tcommon.
ELEMENT_ORDER = ["H", "O", "N"]


def thermo_record(specie: dict) -> str:
    name = specie["name"]
    composition = specie["composition"]
    formula = ""
    for element in ELEMENT_ORDER:
        if element in composition:
            formula += f"{element:<2s}{int(composition[element]):3d}"
    formula = formula.ljust(20)
    thermo = specie["thermo"]
    low, common, high = thermo["temperature-ranges"]
    low = LOW_LIMIT_OVERRIDE.get(name, low)
    low_coeffs, high_coeffs = thermo["data"]
    head = f"{name:<18s}{'':<6s}{formula}G{low:10.3f}{high:10.3f}{common:8.2f}{'':>6s}1"
    coeffs = list(high_coeffs) + list(low_coeffs)
    lines = [head]
    for row, label in zip((coeffs[0:5], coeffs[5:10], coeffs[10:14]), ("2", "3", "4")):
        text = "".join(f"{value:15.8E}" for value in row)
        lines.append(f"{text:<75s}{label:>5s}")
    return "\n".join(lines)


def arrhenius(rate: dict) -> str:
    return f"{rate['A']:.10g} {rate['b']:.10g} {rate['Ea']:.10g}"


def efficiencies_line(effs: dict) -> str:
    keep = {k: v for k, v in effs.items() if k not in DROPPED}
    return " ".join(f"{k}/{v:.4f}/" for k, v in keep.items())


def reaction_lines(reaction: dict) -> list[str]:
    equation = reaction["equation"]
    kind = reaction.get("type", "elementary")
    lines = []
    if kind == "falloff":
        lines.append(f"{equation}  {arrhenius(reaction['high-P-rate-constant'])}")
        low = reaction["low-P-rate-constant"]
        lines.append(f"LOW /{low['A']:.10g} {low['b']:.10g} {low['Ea']:.10g}/")
        troe = reaction["Troe"]
        lines.append(f"TROE /{troe['A']:.10g} {troe['T3']:.10g} {troe['T1']:.10g}/")
    else:
        lines.append(f"{equation}  {arrhenius(reaction['rate-constant'])}")
    if "efficiencies" in reaction:
        lines.append(efficiencies_line(reaction["efficiencies"]))
    if reaction.get("duplicate"):
        lines.append("DUPLICATE")
    return lines


def main() -> int:
    source = yaml.safe_load((HERE / "source" / "Burke2012.yaml").read_text())
    species = [s for s in source["species"] if s["name"] in SPECIES]
    species.sort(key=lambda s: SPECIES.index(s["name"]))
    reactions = []
    for reaction in source["reactions"]:
        tokens = set(re.split(r"[\s+()<=>]+", reaction["equation"]))
        if tokens & DROPPED:
            continue
        reactions.append(reaction)
    if len(reactions) != 23:
        print(f"expected 23 reactions, found {len(reactions)}", file=sys.stderr)
        return 1

    build = HERE / "build"
    build.mkdir(exist_ok=True)
    inp = ["ELEM", "H O N", "END", "", "SPECIES", "  ".join(SPECIES), "END", "",
           "REACTIONS CAL/MOLE MOLES"]
    for reaction in reactions:
        inp.extend(reaction_lines(reaction))
    inp.append("END")
    (build / "chem.inp").write_text("\n".join(inp) + "\n")

    dat = ["THERMO ALL", "   200.000  1000.000  5000.000"]
    dat.extend(thermo_record(s) for s in species)
    dat.append("END")
    (build / "thermo.dat").write_text("\n".join(dat) + "\n")
    # chemkinToFoam wants a transport dictionary. The solver is Euler only, so
    # the Sutherland coefficient As is 0, which makes mu exactly 0.
    (build / "transport.dict").write_text(
        '".*"\n{\n    transport\n    {\n        As 0;\n        Ts 1;\n    }\n}\n'
    )

    subprocess.run(
        [
            "chemkinToFoam",
            str(build / "chem.inp"),
            str(build / "thermo.dat"),
            str(build / "transport.dict"),
            str(HERE / "reactions"),
            str(HERE / "thermo"),
        ],
        check=True,
    )
    print("Wrote reactions and thermo.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Detonation-tube validation against reference H2-air data.

    python3 -m rdesim.validation analytic [--config PATH]
    python3 -m rdesim.validation calibrate [--write] [--config PATH]
    python3 -m rdesim.validation run [--all | --runs NAME ...] [--measure-only] [--config PATH]
    python3 -m rdesim.validation compare [--config PATH]

``analytic`` and ``calibrate`` need only Python. ``run`` needs OpenFOAM 14 on
PATH (``blockMesh`` and ``foamRun``) and the built rdeFluid library. Every
number the tool uses comes from ``cases/detonationTubeH2/validationConfig``,
``referenceData`` or the fuel files. See README.md, section "Validation".
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

from .detonation import (
    Chemistry,
    Gas,
    calibrate,
    cj_state,
    reduced_activation_energy,
    von_neumann,
    znd_structure,
)
from .foamdict import FoamDict, foam_header, format_data, format_number, parse_file, render_entries
from .mixture import mixture_from_config

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "cases" / "detonationTubeH2" / "validationConfig"


class ValidationError(Exception):
    pass


class Setup:
    def __init__(self, path: Path):
        self.path = path.resolve()
        self.raw = parse_file(self.path)
        for key in (
            "reference", "referenceMechanism", "fuelLibrary", "constants", "operating",
            "fuels", "calibration", "tube", "initiation", "numerics", "hyperbolicity",
            "measurement", "runs", "defaultRuns", "tolerances",
        ):
            if key not in self.raw:
                raise ValidationError(f"{self.path}: missing required entry '{key}'")
        ref_path = self.path.parent / str(self.raw["reference"])
        if not ref_path.is_file():
            raise ValidationError(
                f"Reference data {ref_path} is missing. It is written by "
                "tools/validation_reference.py (needs cantera)."
            )
        reference = parse_file(ref_path)
        mech = str(self.raw["referenceMechanism"])
        if mech not in reference["computed"]:
            raise ValidationError(f"referenceData.computed has no block '{mech}'")
        self.reference = reference
        self.target = reference["computed"][mech]
        self.case_dir = self.path.parent
        self.run_dir = self.case_dir / "run"
        self.library = self._find_library(str(self.raw["fuelLibrary"]))
        op = self.raw["operating"]
        computed = reference["computed"]
        if abs(float(op["T0"]) - float(computed["initialTemperature"])) > 1e-9 or abs(
            float(op["p0"]) - float(computed["initialPressure"])
        ) > 1e-6:
            raise ValidationError(
                "operating.p0/T0 in validationConfig differ from referenceData.computed. "
                "Regenerate the reference with tools/validation_reference.py or fix the config."
            )

    def _find_library(self, name: str) -> Path:
        for parent in [self.path.parent, *self.path.parents]:
            found = parent / name
            if found.is_dir():
                return found
        raise ValidationError(f"Fuel library '{name}' not found above {self.path}")

    def fuel_file(self, name: str) -> Path:
        path = self.library / name / "fuelDict"
        if not path.is_file():
            raise ValidationError(f"Fuel '{name}' not found at {path}")
        return path

    def fuel_name(self, role: str) -> str:
        fuels = self.raw["fuels"]
        if role not in fuels:
            raise ValidationError(f"validationConfig fuels has no entry '{role}'")
        return str(fuels[role])

    def load_fuel(self, name: str) -> FoamDict:
        data = parse_file(self.fuel_file(name))
        return FoamDict((k, v) for k, v in data.items() if k != "FoamFile")


class Model:
    """One fuel as the solver sees it, plus its analytic detonation numbers."""

    def __init__(self, setup: Setup, fuel_name: str):
        self.setup = setup
        self.fuel_name = fuel_name
        self.fuel = setup.load_fuel(fuel_name)
        op = setup.raw["operating"]
        cfg = {
            "fuel": self.fuel,
            "raw": {
                "operating": {
                    "equivalenceRatio": op["equivalenceRatio"],
                    "pInlet": op["p0"],
                    "TInlet": op["T0"],
                },
                "constants": {"RUniversal": setup.raw["constants"]["RUniversal"]},
                "numerics": {"flowModel": setup.raw["numerics"]["flowModel"]},
            },
        }
        self.mixture = mixture_from_config(cfg)
        m = self.mixture
        self.gas = Gas(m.gamma, m.specific_gas, float(op["T0"]), float(op["p0"]), m.heat_release)
        c = self.fuel["chemistry"]
        self.chem = Chemistry(
            float(c["A"]), float(c["beta"]), float(c["Ta"]), float(c["TRef"]),
            float(c["TFloor"]), float(c["ignitionTemperature"]), float(c["lambdaMax"]),
        )
        self.cj = cj_state(self.gas)
        self.vn = von_neumann(self.gas, self.cj["D"])
        points = int(setup.raw["calibration"]["zndPoints"])
        self.znd = znd_structure(self.gas, self.chem, self.cj["D"], points=points)
        delta = float(setup.raw["calibration"]["thetaShockSpeedDelta"])
        self.theta = (
            reduced_activation_energy(self.gas, self.chem, self.cj["D"], delta, "t50", points // 2)
            if math.isfinite(self.znd["x50"]) else math.nan
        )

    @property
    def x_ind(self) -> float:
        return self.znd["x_ind"]

    @property
    def x_half(self) -> float:
        """Length to half the heat release in the model's own steady ZND profile."""
        return self.znd["x50"]


def _ref(setup: Setup) -> dict:
    t = setup.target
    return {
        "D": float(t["DCJ"]), "p": float(t["pCJ"]), "T": float(t["TCJ"]),
        "pvN": float(t["pvN"]), "Tvn": float(t["TvN"]),
        "x_ind": float(t["xIndTmax"]), "theta": float(t["theta"]),
        "x50": float(t["xHeatRelease50"]), "theta50": float(t["thetaHalf"]),
        "x90": float(t["xHeatRelease90"]),
    }


# --------------------------------------------------------------------------
# Calibration


def target_theta(setup: Setup) -> float:
    """Reduced activation energy the calibration aims for, from the config."""
    value = setup.raw["calibration"]["reducedActivationEnergy"]
    if isinstance(value, str):
        if value != "reference":
            raise ValidationError("calibration.reducedActivationEnergy must be a number or the word 'reference'")
        return _ref(setup)["theta50"]
    return float(value)


def calibrate_fuel(setup: Setup) -> dict:
    ref = _ref(setup)
    base = Model(setup, setup.fuel_name("original"))
    delta = float(setup.raw["calibration"]["thetaShockSpeedDelta"])
    result = calibrate(base.gas, base.chem, ref["D"], ref["x50"], target_theta(setup), delta)
    result["LHV"] = result["q"] / base.mixture.y_fuel
    result["yFuel"] = base.mixture.y_fuel
    return result


def _six(value: float) -> str:
    return f"{value:.6g}"


def render_calibrated_fuel(setup: Setup, result: dict) -> str:
    original_name = setup.fuel_name("original")
    text = setup.fuel_file(original_name).read_text(encoding="utf-8")
    base = Model(setup, original_name)
    ref = _ref(setup)
    c = base.fuel["chemistry"]

    def swap(block: str, key: str, value: str, text: str) -> str:
        pattern = re.compile(r"(^\s*" + block + r"\s*\{.*?^\s*" + key + r"\s+)([^;]+)(;)", re.S | re.M)
        new, count = pattern.subn(lambda m: m.group(1) + value + m.group(3), text, count=1)
        if count != 1:
            raise ValidationError(f"Could not find {block}.{key} in {setup.fuel_file(original_name)}")
        return new

    text = swap("heatRelease", "LHV", _six(result["LHV"]), text)
    text = swap("chemistry", "A", _six(result["A"]), text)
    text = swap("chemistry", "Ta", _six(result["Ta"]), text)
    text = re.sub(
        r"// Illustrative hydrogen-air library entry\..*?\n\n",
        "// Calibrated copy of fuels/%s/fuelDict for the detonation-tube validation.\n"
        "// Only heatRelease.LHV, chemistry.A and chemistry.Ta differ from the original.\n"
        "// Written by `python3 -m rdesim.validation calibrate --write`; see the\n"
        "// calibration block at the end of this file. Valid for stoichiometric H2-air\n"
        "// at 300 K and 1 atm. Not a predictive detonation chemistry.\n\n" % original_name,
        text, count=1, flags=re.S,
    )
    block = f"""
// What was calibrated. Not read by the solver. The original values are kept
// here and in fuels/{original_name}/fuelDict.
calibration
{{
    reference       "Cantera {setup.reference['computed']['gri30'].get('mechanism', 'gri30.yaml')} CJ state and steady ZND profile, tools/validation_reference.py";
    mixture         "stoichiometric H2-air, p0 {_six(setup.raw['operating']['p0'])} Pa, T0 {_six(setup.raw['operating']['T0'])} K";
    targets
    {{
        waveSpeed           {_six(ref['D'])};   // m/s, CJ speed
        halfHeatReleaseLength {_six(ref['x50'])};   // m, distance to half the heat release in the ZND profile
        reducedActivationEnergy {_six(target_theta(setup))};   // Ea/(R T_vN) of the time to half heat release
        referenceReducedActivationEnergy {_six(ref['theta50'])};   // the Cantera value, for comparison
    }}
    original
    {{
        LHV     {_six(float(base.fuel['heatRelease']['LHV']))};
        A       {_six(float(c['A']))};
        Ta      {_six(float(c['Ta']))};
    }}
    calibrated
    {{
        LHV     {_six(result['LHV'])};
        A       {_six(result['A'])};
        Ta      {_six(result['Ta'])};
    }}
    note            "LHV is q/Yf with q set by the CJ speed for one gamma. It is an effective number, about half the real H2 lower heating value, that absorbs what one gamma of 1.4 cannot represent. It is not a heat of reaction.";
}}
"""
    marker = "// ************************************************************************* //"
    head, sep, tail = text.rpartition(marker)
    if not sep:
        raise ValidationError("Original fuelDict has no closing banner line")
    return head.rstrip() + "\n" + block + "\n" + sep + tail


def cmd_calibrate(setup: Setup, write: bool) -> int:
    start = time.time()
    result = calibrate_fuel(setup)
    base = setup.fuel_name("original")
    base_chem = parse_file(setup.fuel_file(base))
    print(f"Calibration of '{base}' against {setup.reference['computed'][str(setup.raw['referenceMechanism'])]['mechanism']} "
          f"({time.time() - start:.1f} s)")
    print(f"  LHV  {float(base_chem['heatRelease']['LHV']):.6g} -> {result['LHV']:.6g} J/kg fuel  (q = {result['q']:.6g} J/kg mixture, Yf = {result['yFuel']:.6g})")
    print(f"  A    {float(base_chem['chemistry']['A']):.6g} -> {result['A']:.6g} 1/s")
    print(f"  Ta   {float(base_chem['chemistry']['Ta']):.6g} -> {result['Ta']:.6g} K")
    new_name = setup.fuel_name("calibrated")
    target = setup.library / new_name / "fuelDict"
    text = render_calibrated_fuel(setup, result)
    if write:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
        print(f"Wrote {target}")
    else:
        committed = target.read_text(encoding="utf-8") if target.is_file() else None
        if committed is None:
            print(f"{target} does not exist. Run with --write to create it.")
        elif committed == text:
            print(f"{target} is identical to this calibration.")
        else:
            print(f"{target} differs from this calibration. Run with --write to overwrite it.")
    return 0


# --------------------------------------------------------------------------
# Analytic table


def _pct(value: float, reference: float) -> float:
    return 100.0 * (value - reference) / reference


def analytic_rows(setup: Setup, models: list[tuple[str, Model]]) -> list[list[str]]:
    ref = _ref(setup)
    rows = [["quantity", "reference"] + [f"{n} (error)" for n, _ in models]]

    def add(label, ref_value, fmt, getter):
        row = [label, fmt.format(ref_value) if ref_value is not None else "-"]
        for _, model in models:
            value = getter(model)
            if value is None or (isinstance(value, float) and not math.isfinite(value)):
                row.append("none")
            elif ref_value is None:
                row.append(fmt.format(value))
            else:
                row.append(f"{fmt.format(value)} ({_pct(value, ref_value):+.1f}%)")
        rows.append(row)

    add("CJ speed D (m/s)", ref["D"], "{:.1f}", lambda m: m.cj["D"])
    add("CJ pressure (atm)", ref["p"] / 101325.0, "{:.2f}", lambda m: m.cj["p"] / 101325.0)
    add("CJ temperature (K)", ref["T"], "{:.0f}", lambda m: m.cj["T"])
    add("vN pressure (atm)", ref["pvN"] / 101325.0, "{:.2f}", lambda m: m.vn["p"] / 101325.0)
    add("vN temperature (K)", ref["Tvn"], "{:.0f}", lambda m: m.vn["T"])
    add("length to 50% heat release (mm)", ref["x50"] * 1e3, "{:.4f}", lambda m: m.znd["x50"] * 1e3)
    add("induction length, max dT/dx (mm)", ref["x_ind"] * 1e3, "{:.4f}", lambda m: m.x_ind * 1e3)
    add("reduced activation energy Ea/(R T_vN)", ref["theta50"], "{:.2f}", lambda m: m.theta)
    add("length to 90% heat release (mm)", ref["x90"] * 1e3, "{:.4f}", lambda m: m.znd["x90"] * 1e3)
    return rows


def print_table(rows: list[list[str]]) -> None:
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    for index, row in enumerate(rows):
        print("  ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)))
        if index == 0:
            print("  ".join("-" * w for w in widths))


def print_reference(setup: Setup) -> None:
    ref = setup.reference
    t = setup.target
    cea = ref["computed"].get("cea")
    print("Reference data: stoichiometric H2-air, 101325 Pa, 300 K")
    print(f"  Cantera {t['mechanism']}: D_CJ {float(t['DCJ']):.1f} m/s, p_CJ {float(t['pCJ'])/101325:.2f} atm, "
          f"T_CJ {float(t['TCJ']):.0f} K, vN {float(t['pvN'])/101325:.2f} atm / {float(t['TvN']):.0f} K")
    if cea:
        print(f"  NASA CEA {cea['version']} (independent CJ check): D_CJ {float(cea['DCJ']):.1f} m/s, "
              f"p_CJ {float(cea['pCJ'])/101325:.2f} atm, T_CJ {float(cea['TCJ']):.0f} K")
    other = [k for k in ref["computed"] if k not in ("cea", str(setup.raw["referenceMechanism"]))
             and isinstance(ref["computed"][k], dict)]
    for key in other:
        o = ref["computed"][key]
        print(f"  Cantera {o['mechanism']} (mechanism spread): D_CJ {float(o['DCJ']):.1f} m/s, "
              f"x_ind {float(o['xIndTmax'])*1e3:.4f} mm")
    print(f"  ZND induction length {float(t['xIndTmax'])*1e3:.4f} mm (max dT/dx), "
          f"{float(t['xIndThermicity'])*1e3:.4f} mm (max thermicity); "
          f"const-P reactor check of the induction time: {float(t['tIndConstPReactor'])*1e6:.4f} us vs ZND {float(t['tIndTmax'])*1e6:.4f} us")
    widths = [float(e["width"]) for e in ref["literature"]["cellWidth"]]
    print(f"  Cell width (Shepherd detonation database): {min(widths):g} to {max(widths):g} mm from "
          f"{len(widths)} data points, 100-101.3 kPa, 293-300 K, phi 1.0-1.02")


def cmd_analytic(setup: Setup) -> int:
    print_reference(setup)
    print()
    models = [
        ("before: " + setup.fuel_name("original"), Model(setup, setup.fuel_name("original"))),
    ]
    cal_file = setup.library / setup.fuel_name("calibrated") / "fuelDict"
    if cal_file.is_file():
        models.append(("after: " + setup.fuel_name("calibrated"), Model(setup, setup.fuel_name("calibrated"))))
    else:
        print(f"{cal_file} is missing. Run `python3 -m rdesim.validation calibrate --write`.")
    print("Steady ZND/CJ theory of the one-step model, no CFD")
    print_table(analytic_rows(setup, models))
    return 0


# --------------------------------------------------------------------------
# Case writer


def _run_entry(setup: Setup, name: str) -> dict:
    for entry in setup.raw["runs"]:
        if str(entry["name"]) == name:
            return entry
    raise ValidationError(f"No run named '{name}' in validationConfig")


def _sizes(setup: Setup, run: dict, model: Model) -> dict:
    if not math.isfinite(model.x_half):
        raise ValidationError(f"Fuel {model.fuel_name} does not react at its CJ speed (T_vN below ignitionTemperature)")
    n_per = int(run["cellsPerHalfLength"])
    tube = setup.raw["tube"]
    init = setup.raw["initiation"]
    cells = int(round(float(tube["lengthInHalfLengths"]) * n_per))
    dx = model.x_half / n_per
    length = cells * dx
    kernel_cells = max(1, int(round(float(init["lengthInHalfLengths"]) * n_per)))
    return {
        "cells": cells, "dx": dx, "length": length, "kernelCells": kernel_cells,
        "kernelLength": kernel_cells * dx, "cellsPerHalfLength": n_per,
    }


def _kernel_state(setup: Setup, model: Model) -> dict:
    init = setup.raw["initiation"]
    p = float(init["pressureOverCJ"]) * model.cj["p"]
    t = float(init["temperatureOverCJ"]) * model.cj["T"]
    return {"p": p, "T": t}


def _dict_file(location: str, name: str, body: str, class_name: str = "dictionary") -> str:
    return foam_header(class_name, name, location) + body.rstrip() + "\n\n// ************************************************************************* //\n"


def write_case(setup: Setup, run: dict, model: Model, out: Path) -> dict:
    sizes = _sizes(setup, run, model)
    num = setup.raw["numerics"]
    m = model.mixture
    cj_speed = model.cj["D"]
    wave_time = (sizes["length"] - sizes["kernelLength"]) / (cj_speed * float(num["endSpeedOverCJ"]))
    end_time = float(num["endFraction"]) * wave_time
    write_interval = float(num["writeEveryHalfLength"]) * model.x_half / cj_speed
    cfl = float(num["CFL"])
    dt0 = float(num["deltaTFraction"]) * cfl * sizes["dx"] / cj_speed
    kernel = _kernel_state(setup, model)

    for sub in ("system", "constant", "0"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    for path in out.iterdir():
        if path.is_dir():
            try:
                value = float(path.name)
            except ValueError:
                continue
            if value > 0:
                shutil.rmtree(path)

    h = float(setup.raw["tube"]["crossSection"])
    length = sizes["length"]
    block = f"""convertToMeters 1;

vertices
(
    (0 0 0)
    ({format_number(length)} 0 0)
    ({format_number(length)} {format_number(h)} 0)
    (0 {format_number(h)} 0)
    (0 0 {format_number(h)})
    ({format_number(length)} 0 {format_number(h)})
    ({format_number(length)} {format_number(h)} {format_number(h)})
    (0 {format_number(h)} {format_number(h)})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({sizes['cells']} 1 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    closedEnd
    {{
        type patch;
        faces ((0 4 7 3));
    }}
    outlet
    {{
        type patch;
        faces ((1 2 6 5));
    }}
    sides
    {{
        type empty;
        faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7));
    }}
);

mergePatchPairs
(
);
"""
    (out / "system" / "blockMeshDict").write_text(_dict_file("system", "blockMeshDict", block), encoding="utf-8")

    control = f"""solver          rdeFluid;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {format_number(end_time)};
deltaT          {format_number(dt0)};
adjustTimeStep  yes;
maxCo           {format_number(cfl)};
maxDeltaT       {format_number(write_interval)};
deltaTFactor    {format_number(float(num['deltaTFactor']))};

writeControl    adjustableRunTime;
writeInterval   {format_number(write_interval)};
purgeWrite      0;
writeFormat     ascii;
writePrecision  {int(num['writePrecision'])};
writeCompression off;
timeFormat      general;
timePrecision   {int(num['timePrecision'])};
runTimeModifiable no;
"""
    (out / "system" / "controlDict").write_text(_dict_file("system", "controlDict", control), encoding="utf-8")

    schemes = f"""fluxScheme      {format_data(num['fluxScheme'])};

ddtSchemes
{{
    default         {format_data(num['ddtScheme'])};
}}

gradSchemes
{{
    default         {format_data(num['gradScheme'])};
}}

divSchemes
{{
    default         none;
}}

laplacianSchemes
{{
    default         {format_data(num['laplacianScheme'])};
}}

interpolationSchemes
{{
    default         none;
    reconstruct(rho) {format_data(num['reconstructScalar'])};
    reconstruct(U)  {format_data(num['reconstructVector'])};
    reconstruct(T)  {format_data(num['reconstructScalar'])};
    reconstruct(e)  {format_data(num['reconstructScalar'])};
    reconstruct(lambda) {format_data(num['reconstructScalar'])};
    flux(rhoU)      {format_data(num['fluxInterpolation'])};
}}

snGradSchemes
{{
    default         orthogonal;
}}
"""
    (out / "system" / "fvSchemes").write_text(_dict_file("system", "fvSchemes", schemes), encoding="utf-8")
    (out / "system" / "fvSolution").write_text(
        _dict_file("system", "fvSolution", render_entries(num["fvSolution"])), encoding="utf-8")

    thermo = f"""thermoType
{{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState {m.equation_of_state};
    specie          specie;
    energy          sensibleInternalEnergy;
}}

mixture
{{
    specie
    {{
        molWeight   {format_number(8314.462618 / m.specific_gas)};
    }}
    thermodynamics
    {{
        Cp          {format_number(m.cp)};
        hf          {format_number(m.hf)};
    }}
    transport
    {{
        mu          0;
        Pr          {format_number(m.prandtl)};
    }}
}}
"""
    (out / "constant" / "physicalProperties").write_text(_dict_file("constant", "physicalProperties", thermo), encoding="utf-8")

    chem = m.chemistry
    hyper = setup.raw["hyperbolicity"]
    props = f"""flowModel          {num['flowModel']};
equationOfState    {m.equation_of_state};
ddtScheme          {format_data(num['ddtScheme'])};

chemistry
{{
    model                  {chem['model']};
    A                      {format_number(float(chem['A']))};
    beta                   {format_number(float(chem['beta']))};
    Ta                     {format_number(float(chem['Ta']))};
    TRef                   {format_number(float(chem['TRef']))};
    TFloor                 {format_number(float(chem['TFloor']))};
    ignitionTemperature    {format_number(float(chem['ignitionTemperature']))};
    lambdaMin              {format_number(float(chem['lambdaMin']))};
    lambdaMax              {format_number(float(chem['lambdaMax']))};
    q                      {format_number(m.heat_release)};
    sourceIntegration      {num['sourceIntegration']};
    sourceDependsOnGradients {hyper['sourceDependsOnGradients']};
}}

hyperbolicity
{{
{render_entries(hyper, indent=1)}
}}

objectives
(
);
"""
    (out / "constant" / "rdeProperties").write_text(_dict_file("constant", "rdeProperties", props), encoding="utf-8")
    (out / "constant" / "momentumTransport").write_text(
        _dict_file("constant", "momentumTransport", "simulationType    laminar;\n"), encoding="utf-8")
    (out / "constant" / "thermophysicalTransport").write_text(
        _dict_file("constant", "thermophysicalTransport",
                   f"laminar\n{{\n    model    Fourier;\n    Pr    {format_number(m.prandtl)};\n}}\n"),
        encoding="utf-8")
    (out / "constant" / "fvConstraints").write_text(
        _dict_file("constant", "fvConstraints", "// No fvConstraints in this case.\n"), encoding="utf-8")
    (out / "constant" / "fvModels").write_text(
        _dict_file("constant", "fvModels", "// No fvModels in this case.\n"), encoding="utf-8")

    lo = format_number(-h)
    hi = format_number(2 * h)
    kernel_end = format_number(sizes["kernelLength"])

    def zonal(default: float, value: float) -> str:
        return f"""{{
    type            zonal;
    defaultValue    {format_number(default)};
    zones
    {{
        kernel
        {{
            type            box;
            box             ({lo} {lo} {lo}) ({kernel_end} {hi} {hi});
            value           {format_number(value)};
        }}
    }}
}}"""

    def field(cls: str, name: str, dims: str, internal: str, left: str, right: str) -> str:
        body = f"""dimensions      {dims};

internalField   {internal};

boundaryField
{{
    sides
    {{
        type            empty;
    }}
    closedEnd
    {{
{left}
    }}
    outlet
    {{
{right}
    }}
}}
"""
        return foam_header(cls, name, "0") + body + "\n// ************************************************************************* //\n"

    zero_grad = "        type            zeroGradient;"
    (out / "0" / "p").write_text(field("volScalarField", "p", "[pressure]",
        zonal(model.gas.p0, kernel["p"]), zero_grad, zero_grad), encoding="utf-8")
    (out / "0" / "T").write_text(field("volScalarField", "T", "[temperature]",
        zonal(model.gas.T0, kernel["T"]), zero_grad, zero_grad), encoding="utf-8")
    (out / "0" / "U").write_text(field("volVectorField", "U", "[velocity]",
        "uniform (0 0 0)", "        type            fixedValue;\n        value           uniform (0 0 0);", zero_grad),
        encoding="utf-8")
    (out / "0" / "lambda").write_text(field("volScalarField", "lambda", "[0 0 0 0 0 0 0]",
        zonal(0.0, float(chem["lambdaMax"])), zero_grad, zero_grad), encoding="utf-8")
    (out / f"{out.name}.foam").write_text("", encoding="utf-8")
    meta = {
        "name": run["name"], "fuel": model.fuel_name, **sizes,
        "endTime": end_time, "writeInterval": write_interval,
        "kernelPressure": kernel["p"], "kernelTemperature": kernel["T"],
    }
    (out / "validationCase.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return meta


# --------------------------------------------------------------------------
# Running OpenFOAM


def _need_openfoam() -> None:
    missing = [c for c in ("blockMesh", "foamRun") if shutil.which(c) is None]
    if missing:
        raise ValidationError(
            f"{' and '.join(missing)} not on PATH. Source OpenFOAM 14 first "
            "(. /opt/openfoam14/etc/bashrc) or use ./bin/rdesim-wsl validate."
        )


def _run_logged(command: list[str], cwd: Path, log: Path) -> float:
    start = time.time()
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=False)
    elapsed = time.time() - start
    if result.returncode != 0:
        tail = "".join(log.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)[-25:])
        raise ValidationError(f"{' '.join(command)} failed in {cwd} (exit {result.returncode}). Last lines of {log}:\n{tail}")
    return elapsed


# --------------------------------------------------------------------------
# Field reading and measurement


_LIST = re.compile(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", re.S)
_UNIFORM = re.compile(r"internalField\s+uniform\s+([^;]+);")
_VECTOR = re.compile(r"\(\s*(\S+)\s+\S+\s+\S+\s*\)")


def read_field(path: Path, cells: int) -> list[float]:
    """Scalar field, or the x component of a vector field, one value per cell."""
    text = path.read_text(encoding="utf-8")
    match = _LIST.search(text)
    if match:
        kind, count = match.group(1), int(match.group(2))
        if count != cells:
            raise ValidationError(f"{path}: {count} values, expected {cells}")
        if kind == "scalar":
            end = text.index(")", match.end())
            return [float(v) for v in text[match.end():end].split()]
        if kind == "vector":
            values = []
            for item in _VECTOR.finditer(text, match.end()):
                values.append(float(item.group(1)))
                if len(values) == count:
                    break
            return values
        raise ValidationError(f"{path}: unsupported List<{kind}>")
    uniform = _UNIFORM.search(text)
    if uniform:
        token = uniform.group(1).strip()
        if token.startswith("("):
            token = token.strip("()").split()[0]
        return [float(token)] * cells
    raise ValidationError(f"{path}: internalField not understood")


def result_times(case: Path) -> list[tuple[float, str]]:
    times = []
    for path in case.iterdir():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value > 0 and (path / "p").is_file():
            times.append((value, path.name))
    return sorted(times)


def _crossing(values: list[float], threshold: float) -> int | None:
    for i in range(len(values) - 1, -1, -1):
        if values[i] > threshold:
            return i
    return None


def measure(setup: Setup, run: dict, model: Model, out: Path, meta: dict) -> dict:
    meas = setup.raw["measurement"]
    cells = int(meta["cells"])
    dx = float(meta["dx"])
    length = float(meta["length"])
    p0 = model.gas.p0
    threshold = float(meas["shockPressureRatio"]) * p0
    gamma = model.gas.gamma
    gas_r = model.gas.R
    skip = int(meas["excludeShockCells"])
    search = int(round(float(meas["searchHalfLengths"]) * meta["cellsPerHalfLength"]))

    frames = []
    for t, name in result_times(out):
        p = read_field(out / name / "p", cells)
        i = _crossing(p, threshold)
        if i is None or i >= cells - 1:
            continue
        xs = (i + 0.5) * dx + (p[i] - threshold) / (p[i] - p[i + 1]) * dx
        frames.append({"t": t, "name": name, "i": i, "xs": xs, "p": p})
    if len(frames) < 4:
        raise ValidationError(f"{out}: fewer than four result times with a shock; the wave did not start")

    lo = float(meas["fitFromFraction"]) * length
    hi = float(meas["fitToFraction"]) * length
    window = [f for f in frames if lo <= f["xs"] <= hi]
    if len(window) < 4:
        raise ValidationError(
            f"{out}: only {len(window)} frames with the shock between {lo:.4g} and {hi:.4g} m. "
            "The wave did not travel through the fit window."
        )
    ts = [f["t"] for f in window]
    xs = [f["xs"] for f in window]
    tm, xm = statistics.fmean(ts), statistics.fmean(xs)
    speed = sum((t - tm) * (x - xm) for t, x in zip(ts, xs)) / sum((t - tm) ** 2 for t in ts)
    local = [(xs[k + 1] - xs[k]) / (ts[k + 1] - ts[k]) for k in range(len(ts) - 1)]

    sonic_p, sonic_T, vn_peak, mach_max, half_len = [], [], [], [], []
    for f in window:
        name = f["name"]
        p = f["p"]
        temp = read_field(out / name / "T", cells)
        u = read_field(out / name / "U", cells)
        i0 = f["i"] - skip
        peak = max(p[max(0, f["i"] - search): f["i"] + 1])
        vn_peak.append(peak)
        best = None
        prev_mach = None
        top = 0.0
        for i in range(i0, -1, -1):
            c = math.sqrt(gamma * gas_r * temp[i])
            mach = (speed - u[i]) / c
            top = max(top, mach)
            if prev_mach is not None and prev_mach < 1.0 <= mach:
                w = (1.0 - prev_mach) / (mach - prev_mach)
                best = (i, i + 1, w)
                break
            prev_mach = mach
        mach_max.append(top)
        if best is not None:
            ia, ib, w = best
            sonic_p.append(p[ib] + w * (p[ia] - p[ib]))
            sonic_T.append(temp[ib] + w * (temp[ia] - temp[ib]))
        lam = read_field(out / name / "lambda", cells)
        for i in range(f["i"], max(f["i"] - search, 0), -1):
            if lam[i] >= 0.5 > lam[i + 1]:
                frac = (0.5 - lam[i + 1]) / (lam[i] - lam[i + 1])
                crossing = ((i + 1) + 0.5) * dx - frac * dx
                half_len.append(f["xs"] - crossing)
                break

    out_dict = {
        "run": run["name"], "fuel": model.fuel_name, "cells": cells, "dx": dx,
        "cellsPerHalfLength": meta["cellsPerHalfLength"], "length": length,
        "frames": len(frames), "fitFrames": len(window),
        "D": speed,
        "Dlocal": {"min": min(local), "max": max(local), "mean": statistics.fmean(local)},
        "pSonic": statistics.fmean(sonic_p) if sonic_p else None,
        "pSonicSpread": (max(sonic_p) - min(sonic_p)) / 2 if sonic_p else None,
        "TSonic": statistics.fmean(sonic_T) if sonic_T else None,
        "sonicFound": len(sonic_p),
        "machRelMax": statistics.fmean(mach_max),
        "pVnPeak": statistics.fmean(vn_peak),
        "xHalf": statistics.median(half_len) if half_len else None,
        "xHalfMin": min(half_len) if half_len else None,
        "xHalfMax": max(half_len) if half_len else None,
    }
    return out_dict


def _clock_seconds(log: Path) -> float:
    last = None
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ExecutionTime"):
            last = line
    if last is None:
        return math.nan
    return float(last.split("ClockTime =")[1].split()[0])


def cmd_run(setup: Setup, names: list[str], measure_only: bool = False) -> int:
    if not measure_only:
        _need_openfoam()
    for name in names:
        run = _run_entry(setup, name)
        model = Model(setup, setup.fuel_name(str(run["fuel"])))
        out = setup.run_dir / name
        if measure_only:
            meta_path = out / "validationCase.json"
            if not meta_path.is_file():
                raise ValidationError(f"{meta_path} is missing; run '{name}' first")
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            t_run = _clock_seconds(out / "log.foamRun")
        else:
            meta = write_case(setup, run, model, out)
            print(f"[{name}] fuel {model.fuel_name}: {meta['cells']} cells, dx {meta['dx']*1e6:.2f} um "
                  f"({meta['cellsPerHalfLength']} per half-release length {model.x_half*1e3:.4f} mm), "
                  f"tube {meta['length']*1e3:.2f} mm, end time {meta['endTime']*1e6:.2f} us")
            shutil.rmtree(out / "constant" / "polyMesh", ignore_errors=True)
            t_mesh = _run_logged(["blockMesh"], out, out / "log.blockMesh")
            print(f"[{name}] blockMesh done ({t_mesh:.1f} s). Running foamRun, log {out / 'log.foamRun'}")
            t_run = _run_logged(["foamRun"], out, out / "log.foamRun")
            print(f"[{name}] foamRun done ({t_run:.1f} s)")
        result = measure(setup, run, model, out, meta)
        result["runtime"] = t_run
        (out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0


# --------------------------------------------------------------------------
# Comparison table


def _verdict(error_pct: float | None, band: float) -> str:
    if error_pct is None:
        return "n/a"
    return "PASS" if abs(error_pct) <= band else "FAIL"


def cmd_compare(setup: Setup, names: list[str] | None = None) -> int:
    ref = _ref(setup)
    tol = setup.raw["tolerances"]
    print_reference(setup)
    print()
    models = [("before: " + setup.fuel_name("original"), Model(setup, setup.fuel_name("original")))]
    cal_file = setup.library / setup.fuel_name("calibrated") / "fuelDict"
    if cal_file.is_file():
        models.append(("after: " + setup.fuel_name("calibrated"), Model(setup, setup.fuel_name("calibrated"))))
    print("Steady ZND/CJ theory of the one-step model (no CFD)")
    print_table(analytic_rows(setup, models))
    print()

    results = []
    for run in setup.raw["runs"]:
        name = str(run["name"])
        if names is not None and name not in names:
            continue
        path = setup.run_dir / name / "result.json"
        if path.is_file():
            data = json.loads(path.read_text(encoding="utf-8"))
            data["gate"] = str(run.get("gate", "no")) == "yes"
            results.append(data)
    if not results:
        print("No CFD results yet. Run: python3 -m rdesim.validation run")
        return 0

    failed = False
    for r in results:
        gate = " (gated)" if r["gate"] else " (reported, not gated)"
        print(f"CFD run {r['run']}{gate}: fuel {r['fuel']}, {r['cells']} cells, dx {r['dx']*1e6:.2f} um, "
              f"{r['cellsPerHalfLength']} cells per model half-release length, runtime {r['runtime']:.0f} s")
        rows = [["quantity", "reference", "CFD", "error", "band", "result"]]
        entries = [
            ("wave speed (m/s)", ref["D"], r["D"], "{:.1f}", float(tol["waveSpeedPercent"])),
            ("pressure at sonic plane (atm)", ref["p"] / 101325.0, None if r["pSonic"] is None else r["pSonic"] / 101325.0, "{:.2f}", float(tol["pressurePercent"])),
            ("temperature at sonic plane (K)", ref["T"], r["TSonic"], "{:.0f}", float(tol["temperaturePercent"])),
            ("length to 50% heat release (mm)", ref["x50"] * 1e3, None if r["xHalf"] is None else r["xHalf"] * 1e3, "{:.4f}", float(tol["halfReleaseLengthPercent"])),
        ]
        for label, rv, cv, fmt, band in entries:
            if cv is None:
                rows.append([label, fmt.format(rv), "not found", "-", f"{band:g}%", "FAIL"])
                if r["gate"]:
                    failed = True
                continue
            err = _pct(cv, rv)
            verdict = _verdict(err, band)
            if verdict == "FAIL" and r["gate"]:
                failed = True
            rows.append([label, fmt.format(rv), fmt.format(cv), f"{err:+.1f}%", f"{band:g}%", verdict])
        print_table(rows)
        print(f"  speed between frames: {r['Dlocal']['min']:.0f} to {r['Dlocal']['max']:.0f} m/s (mean {r['Dlocal']['mean']:.0f}); "
              f"max relative Mach number behind the shock {r['machRelMax']:.3f}; "
              f"sonic plane found in {r['sonicFound']}/{r['fitFrames']} frames; "
              f"p spike behind shock {r['pVnPeak']/101325:.1f} atm (reference vN {ref['pvN']/101325:.1f} atm, not gated)")
        own = Model(setup, r["fuel"])
        own_rows = [
            ("D", r["D"], own.cj["D"], "{:.1f} m/s"),
            ("p", None if r["pSonic"] is None else r["pSonic"], own.cj["p"], "{:.2f} atm"),
            ("T", r["TSonic"], own.cj["T"], "{:.0f} K"),
        ]
        parts = []
        for label, cfd, theory, fmt in own_rows:
            if cfd is None:
                continue
            scale = 1.0 / 101325.0 if "atm" in fmt else 1.0
            parts.append(f"{label} {fmt.format(cfd * scale)} vs {fmt.format(theory * scale)} ({_pct(cfd, theory):+.1f}%)")
        print("  solver check against the model's own CJ theory (not the reference): " + "; ".join(parts))
        if r["xHalfMin"] is not None:
            print(f"  half-release length over frames: {r['xHalfMin']*1e3:.4f} to {r['xHalfMax']*1e3:.4f} mm")
        print()
    print("Not compared, because the 1D one-step model cannot or does not predict it: detonation cell width "
          "(reference 8 to 15 mm), reaction-zone length beyond the induction zone, multi-step chemistry, "
          "radial or transverse structure.")
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python3 -m rdesim.validation")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("analytic", help="Theory for the original and calibrated fuels, no OpenFOAM")
    cal = sub.add_parser("calibrate", help="Fit q, Ta and A to the reference")
    cal.add_argument("--write", action="store_true", help="Write the calibrated fuelDict")
    run = sub.add_parser("run", help="Mesh, run foamRun, measure and print the table")
    run.add_argument("--all", action="store_true", help="Every run in validationConfig")
    run.add_argument("--runs", nargs="+", help="Run names")
    run.add_argument("--measure-only", action="store_true", help="Re-measure results already on disk, no OpenFOAM")
    sub.add_parser("compare", help="Print the table from saved results")
    args = parser.parse_args(argv)
    try:
        setup = Setup(Path(args.config))
        if args.command == "analytic":
            return cmd_analytic(setup)
        if args.command == "calibrate":
            return cmd_calibrate(setup, args.write)
        if args.command == "compare":
            return cmd_compare(setup)
        if args.all:
            names = [str(r["name"]) for r in setup.raw["runs"]]
        elif args.runs:
            names = args.runs
        else:
            names = [str(n) for n in setup.raw["defaultRuns"]]
        cmd_run(setup, names, args.measure_only)
        return cmd_compare(setup, names)
    except ValidationError as exc:
        print(f"validation: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

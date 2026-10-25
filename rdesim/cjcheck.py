"""One-dimensional Chapman-Jouguet sanity check for the reaction models.

``python3 -m rdesim.cjcheck --config cases/annulus2D/rdeConfig --model twoStep``
writes a straight 1D tube next to the usual case files (the same thermo,
rdeProperties, schemes and solution files that ``prepare`` writes, so the
selected reaction model is the one that runs), ignites the closed end with a hot
kernel, runs ``foamRun``, tracks the reaction front through the written times,
and compares its speed with the Chapman-Jouguet speed.

Reference used: for a perfect gas that releases q per unit mass, with
H = (gamma^2 - 1) q / (2 c0^2), the Chapman-Jouguet Mach number is
M = sqrt(1 + H) + sqrt(H) and p_CJ / p0 = (1 + gamma M^2) / (gamma + 1)
(Zeldovich-von Neumann-Doering, see Fickett and Davis, Detonation, 1979,
chapter 2). It depends only on gamma and q, not on the rates, so the two-step
and one-step models share it and a mismatch points at the heat release or at the
resolution.

What this is: a check that a 1D wave settles near the CJ speed for the chosen
resolution. What it is not: a validation against measured detonation data, and
it says nothing about cell size or about the multi-species model's equilibrium.
Run from an OpenFOAM 14 shell (``. /opt/openfoam14/etc/bashrc``) or through
``./bin/rdesim-wsl cjcheck``.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

from .casewriter import prepare_case
from .config import ConfigError, load_config
from .foamdict import FoamDict, write_file
from .mixture import mixture_from_config


def cj_perfect_gas(gamma: float, q: float, c0: float, p0: float) -> dict:
    """CJ Mach number, speed and pressure for one perfect gas with heat release q."""
    if gamma <= 1 or q < 0 or c0 <= 0 or p0 <= 0:
        raise ValueError("gamma must exceed 1 and q, c0, p0 must be positive")
    h = (gamma * gamma - 1.0) * q / (2.0 * c0 * c0)
    mach = math.sqrt(1.0 + h) + math.sqrt(h)
    return {
        "mach": mach,
        "speed": mach * c0,
        "pressureRatio": (1.0 + gamma * mach * mach) / (gamma + 1.0),
        "pressure": p0 * (1.0 + gamma * mach * mach) / (gamma + 1.0),
    }


def tube_block_mesh(length: float, cells: int, width: float = 1.0e-3) -> str:
    """blockMeshDict for a 1D tube along x. y and z faces are empty patches."""
    return f"""scale   1;
vertices
(
    (0 0 0) ({length:.10g} 0 0) ({length:.10g} {width:.10g} 0) (0 {width:.10g} 0)
    (0 0 {width:.10g}) ({length:.10g} 0 {width:.10g}) ({length:.10g} {width:.10g} {width:.10g}) (0 {width:.10g} {width:.10g})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) ({cells} 1 1) simpleGrading (1 1 1)
);
edges ();
boundary
(
    inletWall {{ type wall; faces ((0 4 7 3)); }}
    outlet    {{ type patch; faces ((1 2 6 5)); }}
    front     {{ type empty; faces ((0 1 5 4) (3 7 6 2)); }}
    back      {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
mergePatchPairs ();
"""


_LIST = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(([^)]*)\)", re.S
)


def read_scalar_field(path: Path) -> list[float] | float:
    """Internal values of a volScalarField: a list, or one number if uniform."""
    text = path.read_text(encoding="utf-8")
    match = _LIST.search(text)
    if match:
        values = [float(item) for item in match.group(2).split()]
        if len(values) != int(match.group(1)):
            raise ValueError(f"{path}: list length does not match its count")
        return values
    uniform = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", text)
    if uniform:
        return float(uniform.group(1))
    raise ValueError(f"{path}: internalField is neither uniform nor a scalar list")


def front_position(values: list[float], dx: float, level: float) -> float | None:
    """Right-most x where ``values`` falls through ``level``, linearly interpolated.

    Used on pressure for the lead shock and on lambda for the reaction front.
    Cell i is centred at (i + 0.5) dx, so the crossing between cells i-1 and i
    lies at (i - 0.5 + fraction) dx.
    """
    for i in range(len(values) - 1, 0, -1):
        if values[i - 1] >= level > values[i]:
            fraction = (values[i - 1] - level) / (values[i - 1] - values[i])
            return (i - 0.5 + fraction) * dx
    return None


def fitted_speed(times: list[float], positions: list[float]) -> float:
    """Least-squares slope of position against time."""
    if len(times) < 2:
        raise ValueError("need at least two points to fit a speed")
    mean_t = sum(times) / len(times)
    mean_x = sum(positions) / len(positions)
    denom = sum((t - mean_t) ** 2 for t in times)
    if denom <= 0:
        raise ValueError("times are all equal")
    return sum((t - mean_t) * (x - mean_x) for t, x in zip(times, positions)) / denom


def write_tube_case(
    config: dict,
    case_dir: Path,
    length: float,
    cells: int,
    kernel_length: float,
    kernel_temperature: float,
    kernel_pressure: float,
    end_time: float,
    write_interval: float,
) -> dict:
    raw = config["raw"]
    raw["ignition"]["type"] = "hotKernel"
    raw["ignition"]["boxMin"] = [-1.0, -1.0, -1.0]
    raw["ignition"]["boxMax"] = [kernel_length, 1.0, 1.0]
    raw["ignition"]["TKernel"] = kernel_temperature
    raw["boundaries"]["walls"]["velocityModel"] = "slip"
    raw["numerics"]["endTime"] = end_time
    raw["output"]["writeInterval"] = write_interval
    raw["output"]["writeControl"] = "adjustableRunTime"
    raw["output"].pop("checkpoint", None)
    raw["output"]["objectives"] = [
        FoamDict([("name", "TMax"), ("field", "T"), ("region", "volume"), ("operation", "max")]),
        FoamDict([("name", "pMax"), ("field", "p"), ("region", "volume"), ("operation", "max")]),
        FoamDict([("name", "lambdaMax"), ("field", "lambda"), ("region", "volume"), ("operation", "max")]),
    ]
    prepared = prepare_case(config, case_dir)
    write_file(
        case_dir / "system" / "blockMeshDict",
        "dictionary",
        "blockMeshDict",
        tube_block_mesh(length, cells),
    )
    pressure_file = case_dir / "0" / "p"
    text = pressure_file.read_text(encoding="utf-8")
    kernel = (
        "internalField   {\n"
        "    type            zonal;\n"
        f"    defaultValue    {float(raw['ignition']['initialPressure']):.10g};\n"
        "    zones\n    {\n        kernel\n        {\n"
        "            type            box;\n"
        f"            box             (-1 -1 -1) ({kernel_length:.10g} 1 1);\n"
        f"            value           {kernel_pressure:.10g};\n"
        "        }\n    }\n};"
    )
    text = re.sub(r"internalField\s+[^;]*;", lambda _m: kernel, text, count=1)
    pressure_file.write_text(text, encoding="utf-8")
    return prepared


def run_openfoam(case_dir: Path) -> None:
    for command in (["blockMesh"], ["foamRun"]):
        log = case_dir / f"log.{command[0]}"
        with log.open("w", encoding="utf-8") as handle:
            result = subprocess.run(
                command, cwd=case_dir, stdout=handle, stderr=subprocess.STDOUT, check=False
            )
        if result.returncode != 0:
            raise RuntimeError(f"{command[0]} failed with exit code {result.returncode}; see {log}")


def analyse(case_dir: Path, length: float, cells: int, expected: dict, p0: float) -> dict:
    dx = length / cells
    times = []
    for entry in case_dir.iterdir():
        if entry.is_dir() and (entry / "lambda").is_file():
            try:
                value = float(entry.name)
            except ValueError:
                continue
            if value > 0:
                times.append((value, entry))
    times.sort()
    points = []
    zone = []
    peak_pressure = 0.0
    for value, entry in times:
        lam = read_scalar_field(entry / "lambda")
        pressure = read_scalar_field(entry / "p")
        if isinstance(lam, float) or isinstance(pressure, float):
            continue
        shock = front_position(pressure, dx, 1.5 * p0)
        if shock is not None:
            points.append((value, shock))
            reaction = front_position(lam, dx, 0.5)
            if reaction is not None and abs(reaction - shock) < 0.25 * length:
                zone.append(abs(reaction - shock))
            peak_pressure = max(peak_pressure, max(pressure))
    if len(points) < 4:
        raise RuntimeError(
            f"Only {len(points)} written times have a lead shock; "
            "lengthen the run or the kernel. No wave speed was measured."
        )
    second_half = points[len(points) // 2:]
    speed = fitted_speed([t for t, _ in second_half], [x for _, x in second_half])
    return {
        "points": points,
        "fittedSpeed": speed,
        "ratio": speed / expected["speed"],
        "peakPressure": peak_pressure,
        "dx": dx,
        "shockToReactionFront": sum(zone[len(zone) // 2:]) / max(len(zone[len(zone) // 2:]), 1)
        if zone
        else None,
    }


def _zone_text(result: dict) -> str:
    zone = result["shockToReactionFront"]
    if zone is None:
        return "not found"
    return f"{zone * 1e6:.1f} micrometres = {zone / result['dx']:.1f} cells (later-half mean)"


def report(model: str, expected: dict, result: dict, gamma: float, q: float, c0: float) -> str:
    lines = [
        f"1D Chapman-Jouguet sanity check, reaction model {model}",
        f"  gamma {gamma:.6g}  q {q:.6g} J/kg  c0 {c0:.6g} m/s",
        f"  CJ Mach number        {expected['mach']:.5f}",
        f"  CJ speed              {expected['speed']:.2f} m/s",
        f"  CJ pressure           {expected['pressure']:.4g} Pa",
        f"  cell size             {result['dx'] * 1e6:.2f} micrometres",
        f"  lead shock positions  {len(result['points'])} written times, speed fitted on the later half",
        f"  measured shock speed  {result['fittedSpeed']:.2f} m/s",
        f"  measured / CJ         {result['ratio']:.4f}",
        f"  shock to lambda=0.5   {_zone_text(result)}",
        f"  largest cell pressure {result['peakPressure']:.4g} Pa (includes the von Neumann spike and any numerical overshoot)",
        "This checks the heat release and the flux scheme at one resolution. It is not a cell-size or data validation.",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python3 -m rdesim.cjcheck")
    parser.add_argument("--config", required=True)
    parser.add_argument("--model", help="Reaction model; default is the one in rdeConfig")
    parser.add_argument("--case", required=True, help="Directory to write the tube case into")
    parser.add_argument("--length", type=float, default=0.12, help="Tube length, m")
    parser.add_argument("--cells", type=int, default=6000)
    parser.add_argument("--kernel-length", type=float, default=0.008, help="Hot kernel length, m")
    parser.add_argument("--kernel-temperature", type=float, default=2500.0)
    parser.add_argument("--kernel-pressure", type=float, default=4.0e6)
    parser.add_argument("--run-fraction", type=float, default=0.85, help="Fraction of the tube the wave crosses")
    parser.add_argument("--no-run", action="store_true", help="Only write the case")
    parser.add_argument("--tolerance", type=float, default=0.05, help="Allowed |measured/CJ - 1|")
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config, args.model)
        state = mixture_from_config(config)
        expected = cj_perfect_gas(state.gamma, state.heat_release, state.sound_speed, state.pressure)
        travel = args.run_fraction * args.length
        end_time = (args.kernel_length + travel) / expected["speed"]
        case_dir = Path(args.case)
        if case_dir.exists():
            shutil.rmtree(case_dir)
        write_tube_case(
            config,
            case_dir,
            args.length,
            args.cells,
            args.kernel_length,
            args.kernel_temperature,
            args.kernel_pressure,
            end_time,
            end_time / 30.0,
        )
        print(f"Wrote tube case {case_dir}")
        if args.no_run:
            return 0
        if shutil.which("foamRun") is None:
            print("foamRun is not on PATH. Run '. /opt/openfoam14/etc/bashrc' first.", file=sys.stderr)
            return 2
        run_openfoam(case_dir)
        result = analyse(case_dir, args.length, args.cells, expected, state.pressure)
    except (ConfigError, RuntimeError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    from .config import reaction_model

    print(report(reaction_model(config), expected, result, state.gamma, state.heat_release, state.sound_speed))
    ok = abs(result["ratio"] - 1.0) <= args.tolerance
    print(f"Within {args.tolerance:.0%} of CJ: {'yes' if ok else 'NO'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

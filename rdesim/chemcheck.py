"""0D check of the detailed chemistry step against a constant-volume reactor.

A short closed tube is filled with the unburned mixture at one temperature and
pressure, so nothing moves and the solver reduces to its chemistry step. The
cell temperature history is compared with Cantera integrating the same
mechanism (the Burke 2012 YAML the OpenFOAM files were built from) in a
constant-volume reactor. Cantera is only needed for this check; it is not a
dependency of the solver or of ``rdesim`` itself.

What this tests: the stiff ODE integration of the OpenFOAM reaction list, the
constant-volume energy bookkeeping and the thermo (janaf) evaluation, with no
advection. What it does not test: shocks, transport, or any 2D behaviour.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from .cjcheck import read_scalar_field, write_tube_case
from .config import ConfigError, load_config
from .species import load_mechanism, mechanism_directory, unburned_mass_fractions

CLOSED_TUBE = """scale   1;
vertices
(
    (0 0 0) ({length:.10g} 0 0) ({length:.10g} 1e-3 0) (0 1e-3 0)
    (0 0 1e-3) ({length:.10g} 0 1e-3) ({length:.10g} 1e-3 1e-3) (0 1e-3 1e-3)
);
blocks
(
    hex (0 1 2 3 4 5 6 7) ({cells} 1 1) simpleGrading (1 1 1)
);
edges ();
boundary
(
    inletWall {{ type wall; faces ((0 4 7 3) (1 2 6 5)); }}
    outlet    {{ type patch; faces (); }}
    front     {{ type empty; faces ((0 1 5 4) (3 7 6 2)); }}
    back      {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
mergePatchPairs ();
"""


def cantera_reference(yaml_path: Path, fractions: dict, temperature: float, pressure: float, times: list[float]):
    """Constant-volume reactor (T, p) at the given times."""
    import cantera as ct

    gas = ct.Solution(str(yaml_path))
    gas.TPY = temperature, pressure, {name: value for name, value in fractions.items() if value > 0}
    reactor = ct.IdealGasReactor(gas, clone=False)
    network = ct.ReactorNet([reactor])
    network.rtol = 1e-10
    network.atol = 1e-14
    series = []
    for t in times:
        network.advance(t)
        series.append((reactor.T, reactor.phase.P))
    return series


def ignition_delay(times: list[float], temperatures: list[float]) -> float | None:
    best = 0.0
    delay = None
    for i in range(len(times) - 1):
        dt = times[i + 1] - times[i]
        if dt > 0:
            slope = (temperatures[i + 1] - temperatures[i]) / dt
            if slope > best:
                best = slope
                delay = 0.5 * (times[i] + times[i + 1])
    return delay


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python3 -m rdesim.chemcheck")
    parser.add_argument("--config", required=True)
    parser.add_argument("--case", required=True, help="Directory to write the 0D case into")
    parser.add_argument("--temperature", type=float, default=1200.0, help="Initial temperature, K")
    parser.add_argument("--pressure", type=float, default=101325.0, help="Initial pressure, Pa")
    parser.add_argument("--end-time", type=float, default=1.0e-4)
    parser.add_argument("--writes", type=int, default=100)
    parser.add_argument("--cells", type=int, default=4)
    parser.add_argument("--delta-t", type=float, default=1.0e-7, help="Solver time step (fixed), s")
    parser.add_argument("--no-run", action="store_true")
    parser.add_argument("--tolerance", type=float, default=0.02, help="Allowed relative error in ignition delay and final T")
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config, "detailed")
        mechanism = load_mechanism(mechanism_directory(config))
        fractions = unburned_mass_fractions(config, mechanism)
        yaml_path = mechanism_directory(config) / "source" / "Burke2012.yaml"
        case_dir = Path(args.case)
        if case_dir.exists():
            shutil.rmtree(case_dir)
        raw = config["raw"]
        raw["numerics"]["adjustTimeStep"] = "no"
        raw["numerics"]["deltaT"] = args.delta_t
        raw["ignition"]["initialTemperature"] = args.temperature
        length = 1.0e-3 * args.cells
        write_tube_case(
            config,
            case_dir,
            length,
            args.cells,
            2.0 * length,
            args.temperature,
            args.pressure,
            args.end_time,
            args.end_time / args.writes,
        )
        (case_dir / "system" / "blockMeshDict").write_text(
            "FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }\n"
            + CLOSED_TUBE.format(length=length, cells=args.cells),
            encoding="utf-8",
        )
        print(f"Wrote 0D case {case_dir}")
        if args.no_run:
            return 0
        if shutil.which("foamRun") is None:
            print("foamRun is not on PATH. Run '. /opt/openfoam14/etc/bashrc' first.", file=sys.stderr)
            return 2
        for command in (["blockMesh"], ["foamRun"]):
            with (case_dir / f"log.{command[0]}").open("w", encoding="utf-8") as handle:
                code = subprocess.run(
                    command, cwd=case_dir, stdout=handle, stderr=subprocess.STDOUT, check=False
                ).returncode
            if code != 0:
                raise RuntimeError(f"{command[0]} exited {code}; see {case_dir}/log.{command[0]}")
        found = []
        for entry in case_dir.iterdir():
            if entry.is_dir() and (entry / "T").is_file():
                try:
                    value = float(entry.name)
                except ValueError:
                    continue
                if value > 0:
                    found.append((value, entry))
        found.sort()
        times, temps, pressures = [], [], []
        for value, folder in found:
            times.append(value)
            for field, store in (("T", temps), ("p", pressures)):
                data = read_scalar_field(folder / field)
                store.append(sum(data) / len(data) if isinstance(data, list) else data)
    except (ConfigError, RuntimeError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    try:
        reference = cantera_reference(yaml_path, fractions, args.temperature, args.pressure, times)
    except ImportError:
        print("cantera is not installed, so there is no reference. Measured final T "
              f"{temps[-1]:.1f} K.", file=sys.stderr)
        return 2
    cantera_t = [item[0] for item in reference]
    cantera_p = [item[1] for item in reference]
    delay_foam = ignition_delay([0.0] + times, [args.temperature] + temps)
    delay_ct = ignition_delay([0.0] + times, [args.temperature] + cantera_t)
    interval = args.end_time / args.writes
    away = [abs(a - b) for t, a, b in zip(times, temps, cantera_t) if abs(t - delay_ct) > 2.0 * interval]
    worst = max(away) if away else 0.0
    delay_error = abs(delay_foam - delay_ct) / delay_ct
    final_error = abs(temps[-1] - cantera_t[-1]) / cantera_t[-1]
    pressure_error = abs(pressures[-1] - cantera_p[-1]) / cantera_p[-1]
    print(
        "\n".join(
            [
                "0D constant-volume check of the detailed chemistry step",
                f"  mixture        phi {config['raw']['operating']['equivalenceRatio']}  T0 {args.temperature:g} K  p0 {args.pressure:g} Pa",
                f"  solver step    {args.delta_t:g} s, {len(times)} written times to {times[-1]:g} s",
                f"  ignition delay OpenFOAM {delay_foam * 1e6:.2f} us   Cantera {delay_ct * 1e6:.2f} us   (maximum dT/dt, resolved to the write interval {args.end_time / args.writes * 1e6:.2f} us)",
                f"  final T        OpenFOAM {temps[-1]:.1f} K   Cantera {cantera_t[-1]:.1f} K   relative error {final_error:.2%}",
                f"  final p        OpenFOAM {pressures[-1]:.5g} Pa   Cantera {cantera_p[-1]:.5g} Pa   relative error {pressure_error:.2%}",
                f"  largest |T difference| more than two write intervals from the ignition  {worst:.1f} K",
                "This checks the chemistry step alone (no shocks, no transport).",
            ]
        )
    )
    ok = delay_error <= args.tolerance and final_error <= args.tolerance
    print(f"Within {args.tolerance:.0%}: {'yes' if ok else 'NO'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

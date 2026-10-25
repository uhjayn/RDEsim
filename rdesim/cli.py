"""Command line for check, prepare, and optimize."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .acoustics import evaluate_acoustics
from .casewriter import prepare_case
from .config import ConfigError, load_config
from .geometry import build_profile
from .hyperbolicity import check_hyperbolicity
from .mixture import mixture_from_config
from .optimize import optimize
from .resolution import resolution_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="rdesim")
    sub = parser.add_subparsers(dest="command", required=True)

    check = sub.add_parser("check", help="Parse a case and report hyperbolicity")
    check.add_argument("--config", required=True)

    prepare = sub.add_parser("prepare", help="Write an OpenFOAM case from rdeConfig")
    prepare.add_argument("--config", required=True)
    prepare.add_argument("--case", required=True)

    opt = sub.add_parser("optimize", help="Vary design variables and run the solver")
    opt.add_argument("--config", required=True)

    resolution = sub.add_parser("resolution", help="Compare cell size with the one-step chemical length")
    resolution.add_argument("--config", required=True)

    args = parser.parse_args(argv)
    try:
        if args.command == "check":
            return _check(Path(args.config))
        if args.command == "prepare":
            return _prepare(Path(args.config), Path(args.case))
        if args.command == "resolution":
            print(resolution_report(load_config(Path(args.config))).text)
            return 0
        return _optimize(Path(args.config))
    except ConfigError as exc:
        print(exc, file=sys.stderr)
        return 1


def _check(path: Path) -> int:
    config = load_config(path)
    state = mixture_from_config(config)
    profile = build_profile(config)
    acoustics = evaluate_acoustics(config, state.density, state.sound_speed)
    report = check_hyperbolicity(config, state)
    print(f"Fuel file: {config['fuelPath']}")
    print(
        f"Mixture: phi={state.phi:.6g} Yf={state.y_fuel:.6g} "
        f"W={state.mol_weight:.6g} kg/mol gamma={state.gamma:.6g} "
        f"q={state.heat_release:.6g} J/kg"
    )
    print(
        f"Geometry: steps={len(profile.step_faces)} "
        f"radialStepsResolved={'yes' if profile.radial_steps_resolved else 'no'} "
        f"injectorSlotsResolved={'yes' if profile.injector_slots_resolved else 'no'}"
    )
    for note in profile.notes:
        print(f"  {note}")
    print(
        f"Acoustics: |R|={acoustics.magnitude:.6g} phase={acoustics.phase_deg:.4g} deg "
        f"absorption={acoustics.absorption:.6g} model={acoustics.time_domain_model}"
    )
    if acoustics.warning:
        print(f"  {acoustics.warning}")
    print(report.text())
    print(resolution_report(config).text)
    return 0 if report.hyperbolic else 1


def _prepare(config_path: Path, case_dir: Path) -> int:
    config = load_config(config_path)
    prepared = prepare_case(config, case_dir)
    print(f"Wrote {case_dir}")
    print(prepared["hyperbolicity"].text())
    return 0 if prepared["hyperbolicity"].hyperbolic else 1


def _optimize(config_path: Path) -> int:
    summary = optimize(config_path)
    print(summary["text"])
    if summary["best"] is None:
        print("No feasible evaluation.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

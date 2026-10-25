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
from .runs import latest_time, run_directory, write_run_record


MODEL_HELP = "Reaction model, overriding reaction.model in rdeConfig"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="rdesim")
    sub = parser.add_subparsers(dest="command", required=True)

    check = sub.add_parser("check", help="Parse a case and report hyperbolicity")
    check.add_argument("--config", required=True)
    check.add_argument("--model", help=MODEL_HELP)

    prepare = sub.add_parser("prepare", help="Write an OpenFOAM case from rdeConfig")
    prepare.add_argument("--config", required=True)
    prepare.add_argument("--model", help=MODEL_HELP)
    prepare.add_argument("--case")
    prepare.add_argument(
        "--run-folder",
        action="store_true",
        help="Write into runs/<case>-<model>-<fingerprint> instead of --case",
    )
    prepare.add_argument("--runs-root", help="Parent of run folders (default: <repo>/runs)")

    run_dir = sub.add_parser("run-dir", help="Print the run folder for a configuration")
    run_dir.add_argument("--config", required=True)
    run_dir.add_argument("--model", help=MODEL_HELP)
    run_dir.add_argument("--runs-root")

    opt = sub.add_parser("optimize", help="Vary design variables and run the solver")
    opt.add_argument("--config", required=True)
    opt.add_argument("--model", help=MODEL_HELP)
    opt.add_argument(
        "--budget",
        type=int,
        help="Run at most this many new solver evaluations; finished ones are reused",
    )
    opt.add_argument(
        "--fresh",
        action="store_true",
        help="Delete earlier evaluations in the work directory first",
    )

    resolution = sub.add_parser("resolution", help="Compare cell size with the one-step chemical length")
    resolution.add_argument("--config", required=True)
    resolution.add_argument("--model", help=MODEL_HELP)

    args = parser.parse_args(argv)
    try:
        if args.command == "check":
            return _check(Path(args.config), args.model)
        if args.command == "prepare":
            if bool(args.case) == bool(args.run_folder):
                parser.error("prepare needs exactly one of --case and --run-folder")
            return _prepare(Path(args.config), args.case, args.run_folder, args.runs_root, args.model)
        if args.command == "resolution":
            print(resolution_report(load_config(Path(args.config), args.model)).text)
            return 0
        if args.command == "run-dir":
            return _run_dir(Path(args.config), args.runs_root, args.model)
        return _optimize(Path(args.config), args.budget, args.fresh, args.model)
    except ConfigError as exc:
        print(exc, file=sys.stderr)
        return 1


def _check(path: Path, model=None) -> int:
    config = load_config(path, model)
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


def _prepare(config_path: Path, case, run_folder: bool, runs_root, model=None) -> int:
    config = load_config(config_path, model)
    if run_folder:
        case_dir = run_directory(config, runs_root)
        known = latest_time(case_dir)
        # The folder name is the physics fingerprint, so resuming is safe.
        config["raw"]["output"].setdefault("startFrom", "latestTime")
    else:
        case_dir = Path(case)
        known = None
    prepared = prepare_case(config, case_dir)
    if run_folder:
        write_run_record(config, case_dir)
    print(f"Wrote {case_dir}")
    if known is not None:
        print(f"Earlier results found up to time {known}. A run with startFrom latestTime resumes there.")
    print(prepared["hyperbolicity"].text())
    return 0 if prepared["hyperbolicity"].hyperbolic else 1


def _run_dir(config_path: Path, runs_root, model=None) -> int:
    print(run_directory(load_config(config_path, model), runs_root))
    return 0


def _optimize(config_path: Path, budget, fresh: bool, model=None) -> int:
    summary = optimize(config_path, budget=budget, fresh=fresh, model=model)
    print(summary["text"])
    if summary["best"] is None:
        print("No feasible evaluation.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

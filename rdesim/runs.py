"""Run folders that keep results per case.

A run folder is named ``<case>-<model>-<fingerprint>`` under the runs root. The
fingerprint covers everything that changes the physics or the mesh: the whole
``rdeConfig`` and the merged fuel dictionary. It leaves out the entries that
only say how long to run or what to write: ``output``, ``optimization``,
``numerics.endTime`` and ``numerics.runTimeModifiable``. The same physics
therefore always lands in the same folder, so a second run resumes from the
latest time written there, and a run with a longer ``endTime`` continues it.
Changing the geometry, the fuel, the reaction model or any numerics entry opens
a new folder and leaves the earlier results alone.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .config import reaction_model

EXCLUDED_TOP = ("output", "optimization")
EXCLUDED_NUMERICS = ("endTime", "runTimeModifiable", "startFrom")


def _plain(value):
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value)
    return str(value)


def physics_snapshot(config: dict) -> dict:
    raw = _plain(config["raw"])
    for key in EXCLUDED_TOP:
        raw.pop(key, None)
    numerics = raw.get("numerics", {})
    for key in EXCLUDED_NUMERICS:
        numerics.pop(key, None)
    return {"raw": raw, "fuel": _plain(config["fuel"])}


def fingerprint(config: dict) -> str:
    text = json.dumps(physics_snapshot(config), sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:10]


def run_name(config: dict) -> str:
    case = Path(config["path"]).resolve().parent.name
    return f"{case}-{reaction_model(config)}-{fingerprint(config)}"


def runs_root(config: dict, override: str | Path | None = None) -> Path:
    if override is not None:
        return Path(override).resolve()
    here = Path(config["path"]).resolve().parent
    for parent in (here, *here.parents):
        if (parent / "rdesim").is_dir() and (parent / "fuels").is_dir():
            return parent / "runs"
    return here / "runs"


def run_directory(config: dict, override: str | Path | None = None) -> Path:
    return runs_root(config, override) / run_name(config)


def write_run_record(config: dict, run_dir: Path) -> None:
    """Keep a copy of the inputs next to the results they produced."""
    run_dir.mkdir(parents=True, exist_ok=True)
    record = {
        "name": run_dir.name,
        "fingerprint": fingerprint(config),
        "reactionModel": reaction_model(config),
        "config": str(config["path"]),
        "fuelFile": str(config["fuelPath"]),
        "preparedAt": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (run_dir / "rdeRun.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    shutil.copyfile(config["path"], run_dir / "rdeConfig.used")
    shutil.copyfile(config["fuelPath"], run_dir / "fuelDict.used")
    pointer = run_dir.parent / "latest"
    pointer.write_text(run_dir.name + "\n", encoding="utf-8")


def latest_time(run_dir: Path) -> str | None:
    best = None
    for path in run_dir.iterdir() if run_dir.is_dir() else ():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value > 0 and (best is None or value > best[0]):
            best = (value, path.name)
    return best[1] if best else None

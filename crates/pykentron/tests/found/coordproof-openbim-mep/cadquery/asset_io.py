"""Shared CadQuery export helpers."""

from __future__ import annotations

import sys
from pathlib import Path

import cadquery as cq

ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = ROOT / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

from reproducibility import source_timestamp  # noqa: E402


def normalize_step_header(path: Path) -> None:
    """Remove host paths and wall-clock time from a CadQuery STEP header."""

    timestamp = source_timestamp()
    lines = path.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        if not line.startswith("FILE_NAME("):
            continue
        fields = line.split("'")
        if len(fields) < 5:
            raise RuntimeError(f"Unexpected STEP FILE_NAME header in {path}")
        fields[1] = path.name
        fields[3] = timestamp
        lines[index] = "'".join(fields)
        break
    else:
        raise RuntimeError(f"STEP FILE_NAME header not found in {path}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_shape(shape: cq.Workplane, asset_id: str) -> tuple[Path, Path]:
    step_path = ROOT / "exports" / "step" / f"{asset_id}.step"
    stl_path = ROOT / "exports" / "stl" / f"{asset_id}.stl"
    step_path.parent.mkdir(parents=True, exist_ok=True)
    stl_path.parent.mkdir(parents=True, exist_ok=True)
    cq.exporters.export(shape, str(step_path))
    normalize_step_header(step_path)
    cq.exporters.export(shape, str(stl_path), tolerance=0.25, angularTolerance=0.2)
    return step_path, stl_path


def merged(parameters: dict[str, object] | None, defaults: dict[str, object]) -> dict[str, object]:
    values = dict(defaults)
    if parameters:
        values.update(parameters)
    return values

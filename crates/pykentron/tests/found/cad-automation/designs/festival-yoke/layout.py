"""Which modules are on the collar and at what rail angle, for the active kit.

Kits are listed in config.toml; angles tuned for balance by `balance.py --tune` are stored in kits.toml (per size
and kit). Pick a kit with `kit = "full"` in config.toml or FESTIVAL_YOKE_KIT=glow in the environment.
"""
import os
import tomllib
from pathlib import Path

import sizing

ROOT = Path(__file__).parent
cfg = tomllib.loads((ROOT / "config.toml").read_text())
SIZE = sizing.resolve(cfg)[0]
KIT = os.environ.get("FESTIVAL_YOKE_KIT") or cfg["kit"]
KITS = cfg["kits"]

# fallback angles (deg on the rail, 0 = right side, 90 = front, 180 = left, 270 = back)
DEFAULTS = {"cup": 40, "phone": 140, "dock": 225, "fan": 345, "pod": 300, "id": 30, "tube": 15, "blank": 200}


def _tuned() -> dict:
    path = ROOT / "kits.toml"
    return tomllib.loads(path.read_text()).get(f"{SIZE}.{KIT}", {}) if path.exists() else {}


def present() -> list[str]:
    return list(KITS[KIT]["modules"])


def angle(name: str) -> float:
    return _tuned().get(name, DEFAULTS[name])


def angles() -> dict[str, float]:
    return {name: angle(name) for name in present()}


def has_leds() -> bool:
    """The pod powers the strip, so a kit with a pod gets the LED strip and cover."""
    return "pod" in present()

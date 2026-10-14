"""Finding font files, for text: by path, or by family name and style among the system's
fonts (fontconfig's `fc-match` where there is one, else the usual font directories), as
OCCT's font manager finds them for CadQuery."""

from __future__ import annotations

import functools
import os
import shutil
import subprocess
import sys
from typing import List, Optional

# The font files kentron can read: TrueType and OpenType (and collections of them).
_EXTENSIONS = (".ttf", ".otf", ".ttc", ".otc")

_STYLES = {"regular": "Regular", "bold": "Bold", "italic": "Italic"}


def font_file(font: str = "Arial", fontPath: Optional[str] = None, kind: str = "regular") -> str:
    """The font file to set text in: `fontPath` if it names a font file, else the system's
    font best matching the family `font` in the style `kind` (`regular`, `bold`, `italic`),
    else any font the system has. Raises `ValueError` if it has none."""
    if kind not in _STYLES:
        raise ValueError(f"unknown font kind {kind!r}: regular, bold or italic")
    if fontPath and os.path.isfile(fontPath) and fontPath.lower().endswith(_EXTENSIONS):
        return fontPath
    found = _by_name(font, kind)
    if found is None:
        raise ValueError(f"no font found for {font!r} ({kind}), nor any other")
    return found


@functools.lru_cache(maxsize=None)
def _by_name(font: str, kind: str) -> Optional[str]:
    style = _STYLES[kind]
    # fontconfig's best match, the first of its fallbacks kentron can read.
    fc = shutil.which("fc-match")
    if fc:
        try:
            out = subprocess.run(
                [fc, "-s", "--format=%{file}\\n", f"{font}:style={style}"],
                capture_output=True, text=True, timeout=10,
            ).stdout
            for line in out.splitlines():
                if line.lower().endswith(_EXTENSIONS) and os.path.isfile(line):
                    return line
        except (OSError, subprocess.SubprocessError):
            pass
    # Else the files in the usual directories, by name: the family's in the style, the
    # family's, a sans serif's, any.
    files = _font_files()
    if not files:
        return None
    family = font.lower().replace(" ", "")
    wanted = {"regular": ("regular", "book", ""), "bold": ("bold",), "italic": ("italic", "oblique")}[kind]

    def stem(path: str) -> str:
        return os.path.splitext(os.path.basename(path))[0].lower().replace(" ", "").replace("-", "").replace("_", "")

    def styled(path: str) -> bool:
        s = stem(path)
        if kind == "regular":
            return not any(w in s for w in ("bold", "italic", "oblique", "light", "thin", "black", "condensed"))
        return any(w in s for w in wanted)

    for choose in (
        lambda p: stem(p).startswith(family) and styled(p),
        lambda p: stem(p).startswith(family),
        lambda p: ("dejavusans" in stem(p) or "liberationsans" in stem(p) or "arial" in stem(p)) and styled(p),
        styled,
        lambda p: True,
    ):
        matches = sorted((p for p in files if choose(p)), key=lambda p: (len(stem(p)), p))
        if matches:
            return matches[0]
    return None


@functools.lru_cache(maxsize=None)
def _font_files() -> List[str]:
    if sys.platform == "win32":
        dirs = [os.path.join(os.environ.get("WINDIR", "C:\\Windows"), "Fonts")]
    elif sys.platform == "darwin":
        dirs = ["/System/Library/Fonts", "/Library/Fonts", os.path.expanduser("~/Library/Fonts")]
    else:
        dirs = ["/usr/share/fonts", "/usr/local/share/fonts", os.path.expanduser("~/.local/share/fonts"), os.path.expanduser("~/.fonts")]
    found = []
    for d in dirs:
        for root, _, names in os.walk(d):
            found.extend(os.path.join(root, n) for n in names if n.lower().endswith(_EXTENSIONS))
    return found

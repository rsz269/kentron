"""Finding font files, for text: by path, or by family name and style among the system's
fonts (fontconfig's `fc-match` where there is one, else the usual font directories), as
OCCT's font manager finds them for CadQuery: a family not installed set in a generic one,
and characters a font lacks in a fallback font for their script."""

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

# The kinds as fontconfig patterns (italic one way or the other, a family's "Oblique" too).
_PATTERNS = {"regular": "style=Regular", "bold": "bold", "italic": "italic"}


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


# A family not installed is set, as OCCT sets it, in a generic one: serif for these names,
# monospaced for these, sans serif for any other; each the first of these families
# installed.
_SERIF_NAMES = ("serif", "times")
_MONO_NAMES = ("mono", "monospace", "courier")
_GENERIC = {
    "sans": ("FreeSans", "DejaVu Sans", "Liberation Sans", "Noto Sans"),
    "serif": ("FreeSerif", "DejaVu Serif", "Liberation Serif", "Noto Serif"),
    "mono": ("FreeMono", "DejaVu Sans Mono", "Liberation Mono", "Noto Sans Mono"),
}


def _fc_match(fc: str, font: str, kind: str, indices: bool = False) -> List[tuple]:
    """fontconfig's matches for a family in a kind of style, best first: each one's family
    names (lower case) and file (and the face's index in the file, with `indices`)."""
    try:
        out = subprocess.run(
            [fc, "-s", "--format=%{family}|%{index}|%{file}\\n", f"{font}:{_PATTERNS[kind]}"],
            capture_output=True, text=True, timeout=10,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return []
    found = []
    for line in out.splitlines():
        families, _, rest = line.partition("|")
        index, _, path = rest.partition("|")
        if path.lower().endswith(_EXTENSIONS) and os.path.isfile(path):
            names = {f.strip().lower() for f in families.split(",")}
            found.append((names, path, int(index or 0)) if indices else (names, path))
    return found


# The fallback fonts for the characters of each script a font lacks, as OCCT's font manager
# has them on this platform: the first of these families installed (the Western one a
# generic sans serif).
_WINDOWS, _MAC = sys.platform == "win32", sys.platform == "darwin"
_FALLBACKS = {
    "cjk": ("simsun", "droid sans fallback", "noto sans sc", "noto sans cjk jp"),
    "korean": (("malgun gothic", "gulim") if _WINDOWS else ("applegothic", "stfangsong") if _MAC else ())
    + ("nanumgothic", "noto sans kr", "nanummyeongjo", "noto serif cjk jp", "noto sans cjk jp", "droid sans fallback"),
    "arabic": ("times new roman",) if _WINDOWS else ("decotype naskh",) if _MAC else (),
}


def _script(c: str) -> str:
    """The script OCCT sets a character in (as kentron's `text::Script::of`)."""
    u = ord(c)
    if 0x3400 <= u <= 0x9FFF or 0xF900 <= u <= 0xFAFF or 0x20000 <= u <= 0x2A6DF or 0x2F800 <= u <= 0x2FA1F or 0x3040 <= u <= 0x30FF:
        return "cjk"
    if 0x1100 <= u <= 0x11FF or 0x3130 <= u <= 0x318F or 0xAC00 <= u <= 0xD7A3:
        return "korean"
    if 0x600 <= u <= 0x6FF:
        return "arabic"
    return "western"


def fallback_fonts(text: str, primary: str, kind: str = "regular") -> List[tuple]:
    """The fallback fonts for the scripts of `text`'s characters, for glyphs the font file
    `primary` lacks: (script, font file, face index) each, as `native.text` takes them."""
    found = []
    for script in sorted({_script(c) for c in text if not c.isspace()}):
        fallback = _fallback(script, kind)
        if fallback and fallback[0] != primary:
            found.append((script, *fallback))
    return found


@functools.lru_cache(maxsize=None)
def _fallback(script: str, kind: str) -> Optional[tuple]:
    if script == "western":
        path = _by_name("Arial", kind)
        return (path, 0) if path else None
    fc = shutil.which("fc-match")
    if not fc:
        return None
    for family in _FALLBACKS[script]:
        matches = _fc_match(fc, family, kind, indices=True)
        if matches and family in matches[0][0]:
            return matches[0][1:]
    return None


@functools.lru_cache(maxsize=None)
def _by_name(font: str, kind: str) -> Optional[str]:
    fc = shutil.which("fc-match")
    if fc:
        # The family itself where it is installed, else a generic one's font, else
        # fontconfig's best match (the first of its fallbacks kentron can read).
        def installed(family: str) -> Optional[str]:
            matches = _fc_match(fc, family, kind)
            return matches[0][1] if matches and family.lower() in matches[0][0] else None

        name = font.strip().lower()
        generic = "serif" if name in _SERIF_NAMES else "mono" if name in _MONO_NAMES else "sans"
        for family in (font, *_GENERIC[generic]):
            path = installed(family)
            if path:
                return path
        matches = _fc_match(fc, font, kind)
        if matches:
            return matches[0][1]
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

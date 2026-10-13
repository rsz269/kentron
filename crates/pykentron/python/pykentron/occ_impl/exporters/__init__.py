"""CadQuery's `occ_impl.exporters` on pykentron."""

from ...exporters import *  # noqa: F401,F403


def __getattr__(name):
    # (The SVG and DXF exporters on first use: DXF needs ezdxf.)
    if name in ("getSVG", "exportSVG"):
        from . import svg

        return getattr(svg, name)
    if name in ("DxfDocument", "exportDXF"):
        from . import dxf

        return getattr(dxf, name)
    raise AttributeError(name)

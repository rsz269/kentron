"""Rendering is not provided: `show` does nothing (CadQuery's needs VTK)."""

import warnings


def show(*objs, **kwargs) -> None:
    warnings.warn("pykentron.vis.show: rendering is not supported; nothing is shown", stacklevel=2)


show_object = show

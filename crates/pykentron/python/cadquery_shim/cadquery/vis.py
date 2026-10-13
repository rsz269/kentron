from pykentron.vis import show, show_object  # noqa: F401


def __getattr__(name):
    """Other names of cadquery.vis: placeholders (see the `cadquery` shim)."""
    if name.startswith("__"):
        raise AttributeError(name)
    import cadquery

    return cadquery._placeholder("cadquery.vis", name)

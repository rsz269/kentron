"""OCP.BRepAdaptor: a face's surface (`BRepAdaptor_Surface`), standing for it as the face
itself (what pykentron's intersections take); the other names are placeholders."""

from ._placeholder import fallback as _fallback


class BRepAdaptor_Surface:
    """A face's surface: `Surface().Surface()` is the face, for `GeomAPI_IntSS`."""

    def __init__(self, face, restriction: bool = True):
        self._face = face

    def Surface(self):
        return self

    def Face(self):
        return self._face

    def __getattr__(self, name):
        raise NotImplementedError(f"OCP.BRepAdaptor.BRepAdaptor_Surface.{name} is not provided by pykentron")


__getattr__ = _fallback("BRepAdaptor")

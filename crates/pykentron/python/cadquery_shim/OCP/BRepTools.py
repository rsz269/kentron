class BRepTools:
    @staticmethod
    def UVBounds_s(face):
        """The face's parameter bounds (u0, u1, v0, v1)."""
        from pykentron import _kentron

        return _kentron.face_uv_bounds(face)

    @staticmethod
    def Clean_s(shape, force: bool = True):
        """Drops meshes stored on the shape. pykentron bounds a shape by its meshes only once that
        shape object itself was tessellated, so there is nothing to undo here."""

from ._placeholder import fallback as _fallback  # noqa: E402

__getattr__ = _fallback("BRepTools")

"""Analyses of shapes beyond CadQuery's API, done natively (for scripts that check what they
build: occupancy grids, ...). CadQuery has none of this: scripts that also run on CadQuery
use it when `pykentron` is there and keep their own way otherwise."""

from __future__ import annotations

from . import _kentron as native


class Mesh:
    """A shape's triangle mesh, kept natively (no Python object per node)."""

    def __init__(self, native_mesh):
        self._native = native_mesh

    @property
    def counts(self) -> tuple[int, int]:
        """How many nodes and triangles."""
        return self._native.counts()

    def arrays(self):
        """The nodes (an n x 3 float array) and the triangles (an m x 3 uint32 array of node
        indices, counter-clockwise seen from outside), as `tessellate` gives them."""
        import numpy as np

        nodes, triangles = self._native.arrays()
        return np.frombuffer(nodes, dtype=float).reshape(-1, 3), np.frombuffer(triangles, dtype=np.uint32).reshape(-1, 3)


def mesh(shape, tolerance: float, angularTolerance: float = 0.1) -> Mesh:
    """The shape's mesh, as `shape.tessellate(tolerance, angularTolerance)` makes it; like
    that, it leaves the shape meshed (its `BoundingBox()` then holds the mesh too, as
    CadQuery's does)."""
    out = Mesh(native.mesh(shape.wrapped, tolerance, angularTolerance))
    shape._meshed = True
    return out


def voxelize(mesh: Mesh, xs, ys, zs):
    """Which points of the grid `xs` x `ys` x `zs` (increasing coordinates) lie inside the
    closed mesh, as a numpy bool array of shape (len(xs), len(ys), len(zs)).

    Each column (x, y) is followed along +Z: its crossings are where the triangles over it
    (their edges included) meet it, and between successive pairs of them, sorted, the
    points from the first (included) to the second (excluded) are inside. A column through
    an edge or vertex of the mesh counts a crossing for each triangle there, so the grid is
    best offset from the mesh's features."""
    import numpy as np

    xs, ys, zs = (np.asarray(a, dtype=float) for a in (xs, ys, zs))
    raw = native.voxelize(mesh._native, xs.tolist(), ys.tolist(), zs.tolist())
    return np.frombuffer(raw, dtype=np.uint8).astype(bool).reshape(len(xs), len(ys), len(zs))


def _morphology(operation: str, mask, structure):
    import numpy as np

    mask = np.ascontiguousarray(mask, dtype=bool)
    structure = np.ascontiguousarray(structure, dtype=bool)
    if mask.ndim != 3 or structure.ndim != 3:
        raise ValueError("3D grids only")
    raw = native.morphology(operation, mask.tobytes(), list(mask.shape), structure.tobytes(), list(structure.shape))
    # (A bytearray: the array is writable, without a copy.)
    return np.frombuffer(raw, dtype=bool).reshape(mask.shape)


def erosion(mask, structure):
    """The cells of the 3D bool grid `mask` whose `structure` (a bool grid of odd sizes,
    centred), placed on them, lies wholly in it; cells outside the grid count as empty (as
    `scipy.ndimage.binary_erosion` by default). A ball-shaped structure (every offset within
    a radius) is applied in time linear in the grid, whatever its radius."""
    return _morphology("erode", mask, structure)


def dilation(mask, structure):
    """The cells whose `structure`, placed on them, meets `mask` (as
    `scipy.ndimage.binary_dilation`)."""
    return _morphology("dilate", mask, structure)


def opening(mask, structure):
    """Erosion, then dilation: `mask` without the parts `structure` does not fit in (as
    `scipy.ndimage.binary_opening`)."""
    return _morphology("open", mask, structure)


def closing(mask, structure):
    """Dilation, then erosion: `mask` with the gaps `structure` does not fit in filled (as
    `scipy.ndimage.binary_closing`)."""
    return _morphology("close", mask, structure)


def components(mask):
    """The connected parts of the 3D bool grid `mask` (cells sharing a face, as
    `scipy.ndimage.label` joins them by default), numbered as it numbers them (by their
    first cells, x-major): their sizes in cells (an int array) and their centres (an n x 3
    float array of mean cell indices, as `np.argwhere(labels == k).mean(0)`)."""
    import numpy as np

    mask = np.ascontiguousarray(mask, dtype=bool)
    if mask.ndim != 3:
        raise ValueError("3D grids only")
    sizes, centres = native.components(mask.tobytes(), list(mask.shape))
    return np.array(sizes, dtype=np.int64), np.array(centres, dtype=float).reshape(len(sizes), 3)

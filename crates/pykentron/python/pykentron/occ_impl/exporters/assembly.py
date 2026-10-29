"""CadQuery's `occ_impl.exporters.assembly` on pykentron: assemblies written to STEP and glTF files."""

from __future__ import annotations

from typing import Literal

from ... import _kentron as native
from ...shapes import Compound, Shape

class ExportModes:
    DEFAULT = "default"
    FUSED = "fused"


STEPExportModeLiterals = Literal["default", "fused"]


def _shapes_of(obj) -> list:
    if obj is None:
        return []
    if isinstance(obj, Shape):
        return [obj]
    if hasattr(obj, "vals"):
        return [o for o in obj.vals() if isinstance(o, Shape)]
    return []


def _subshapes(assy) -> list:
    """The node's faces, edges and wires with names, colours or layers of their own, as the
    native writer takes them: `(shape, name, (r, g, b, a), layer)`."""
    marks = {}
    for attribute, slot in (("_subshape_names", 1), ("_subshape_colors", 2), ("_subshape_layers", 3)):
        for s, value in getattr(assy, attribute, {}).items():
            entry = marks.setdefault(id(s), [s, None, None, None])
            entry[slot] = tuple(value.toTuple()) if slot == 2 else value
    return [(s.wrapped, name, color, layer) for s, name, color, layer in marks.values()]


def _material(assy):
    m = getattr(assy, "material", None)
    return None if m is None else (m.name, m.description, float(m.density), m.densityUnit)


def _tree(assy, factor: float = 1.0):
    """The assembly as the native writer's nested `(name, shape, placement rows, rgba,
    children, extras)` tuples, lengths scaled by `factor` (shapes and placements'
    translations); the extras its subshapes' names, colours and layers and its material."""
    shapes = _shapes_of(assy.obj)
    shape = None
    if shapes:
        whole = shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes)
        shape = (whole.scale(factor) if factor != 1.0 else whole).wrapped
    rgba = tuple(assy.color.toTuple()) if assy.color is not None else None
    rows = [list(r[:3]) + [r[3] * factor] for r in assy.loc.rows]
    # (Subshapes are found by identity: on a scaled copy, none.)
    extras = {"subshapes": _subshapes(assy) if factor == 1.0 else [], "material": _material(assy)}
    return (assy.name, shape, rows, rgba, [_tree(ch, factor) for ch in assy.children], extras)


def _fused_tree(assy, glue: bool = False, tol=None):
    """The assembly's parts, placed, fused into one shape (as CadQuery's `toFusedCAF`), as one
    product of the assembly's name: each face of the fused shape coloured as the part it
    lies on (the part with a face nearest it)."""
    from ...shapes import Vertex

    parts = [(shape.moved(loc), color) for shape, _, loc, color in assy]
    if not parts:
        raise ValueError(f"Error: Assembly {assy.name} has no shapes.")
    if len(parts) == 1:
        fused = parts[0][0]
        if fused.ShapeType() == "Compound":
            fused = fused.fuse(glue=glue, tol=tol) if len(fused.Solids()) > 1 else fused
    else:
        fused = parts[0][0].fuse(*[p for p, _ in parts[1:]], glue=glue, tol=tol)
    subshapes = []
    coloured = [(p, c) for p, c in parts if c is not None]
    if coloured:
        for face in fused.Faces():
            at = face.Center()
            point = Vertex.makeVertex(at.x, at.y, at.z)
            _, color = min(coloured, key=lambda pc: min(f.distance(point) for f in pc[0].Faces()))
            subshapes.append((face.wrapped, None, tuple(color.toTuple()), None))
    identity = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
    return (assy.name, fused.wrapped, identity, None, [], {"subshapes": subshapes, "material": None})


def exportAssembly(assy, path: str, mode: str = "default", unit: str = "MM", outputUnit=None, **kwargs) -> bool:
    """Writes the assembly to a STEP file, its tree as products and occurrences; fused, as
    one product of all its parts fused, its faces coloured as theirs."""
    import os

    from ...exporters import _mm_per

    # (As OCCT's writer: nothing written to a directory that does not exist, and False.)
    folder = os.path.dirname(os.path.abspath(str(path)))
    if not os.path.isdir(folder):
        return False
    out = (outputUnit or unit).upper()
    if mode == "fused":
        if _mm_per(unit) != _mm_per(out):
            raise NotImplementedError("fused assemblies written in another unit")
        native.write_step_assembly(str(path), [_fused_tree(assy, glue=bool(kwargs.get("glue", False)), tol=kwargs.get("fuzzy_tol"))], unit=out)
        return True
    tree = _tree(assy, _mm_per(unit) / _mm_per(out))
    # (STEP places no root: a top level placed otherwise than at the origin is written as
    # the one child, so placed, of a root of its name, which reading takes back off.)
    if any(abs(v - (1.0 if i == j else 0.0)) > 1e-12 for i, row in enumerate(assy.loc.rows) for j, v in enumerate(row)):
        identity = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
        tree = (assy.name, None, identity, None, [tree], {"subshapes": [], "material": None})
    native.write_step_assembly(str(path), [tree], unit=out)
    return True


def exportStepMeta(assy, path: str, unit: str = "MM", outputUnit=None, **kwargs) -> bool:
    return exportAssembly(assy, path, unit=unit, outputUnit=outputUnit)


def exportCAF(assy, path: str, binary: bool = True) -> bool:
    raise NotImplementedError(f"exporting assemblies to {'XBF' if binary else 'XML'} files")


def exportVRML(assy, path: str, tolerance: float = 1e-3, angularTolerance: float = 0.1) -> None:
    raise NotImplementedError("exporting assemblies to VRML files")


def _quaternion(rows) -> list:
    """The rotation of a placement's rows (3 x 4) as a unit quaternion (x, y, z, w)."""
    import math

    (a, b, c), (d, e, f), (g, h, i) = (r[:3] for r in rows)
    trace = a + e + i
    if trace > 0.0:
        k = 0.5 / math.sqrt(trace + 1.0)
        q = [(h - f) * k, (c - g) * k, (d - b) * k, 0.25 / k]
    elif a > e and a > i:
        k = 2.0 * math.sqrt(1.0 + a - e - i)
        q = [0.25 * k, (b + d) / k, (c + g) / k, (h - f) / k]
    elif e > i:
        k = 2.0 * math.sqrt(1.0 + e - a - i)
        q = [(b + d) / k, 0.25 * k, (f + h) / k, (c - g) / k]
    else:
        k = 2.0 * math.sqrt(1.0 + i - a - e)
        q = [(c + g) / k, (f + h) / k, 0.25 * k, (d - b) / k]
    n = math.sqrt(sum(x * x for x in q))
    return [x / n for x in q]


def exportGLTF(assy, path: str, binary=None, tolerance: float = 1e-3, angularTolerance: float = 0.1) -> bool:
    """Writes the assembly to a glTF file (binary, .glb, unless the path ends in .gltf; a
    .gltf's data in a .bin beside it), as CadQuery's (OCCT's glTF writer): a node for each
    of the assembly's, placed as it is, with a mesh of its shape, a primitive for each face
    (positions and normals), in millimetres; its root turned so that +Z up becomes glTF's
    +Y up; a material for each colour (children without one taking their parent's)."""
    import json
    import math
    import os
    import struct
    from array import array

    if binary is None:
        binary = not str(path).lower().endswith(".gltf")
    data = bytearray()
    accessors, views, meshes, materials, nodes = [], [], [], [], []
    material_of = {}

    def view(blob: bytes, target: int, stride=None) -> int:
        while len(data) % 4:
            data.append(0)
        v = {"buffer": 0, "byteOffset": len(data), "byteLength": len(blob), "target": target}
        if stride:
            v["byteStride"] = stride
        data.extend(blob)
        views.append(v)
        return len(views) - 1

    def material(rgba) -> int:
        key = tuple(round(x, 6) for x in rgba)
        if key not in material_of:
            m = {"name": f"mat_{len(materials)}", "pbrMetallicRoughness": {"baseColorFactor": list(key)}, "doubleSided": True}
            if key[3] < 1.0:
                m["pbrMetallicRoughness"]["metallicFactor"] = 0.0
                m["alphaMode"] = "BLEND"
            materials.append(m)
            material_of[key] = len(materials) - 1
        return material_of[key]

    def mesh(name, shape, rgba):
        primitives = []
        for face in Shape.cast(shape).Faces():
            points, tris = native.triangulate(face.wrapped, tolerance, angularTolerance)
            if not tris:
                continue
            # (Each point's normal: the triangles' round it, weighted by their areas.)
            normals = [[0.0, 0.0, 0.0] for _ in points]
            for t in tris:
                p, q, r = (points[k] for k in t)
                u, v = [q[j] - p[j] for j in range(3)], [r[j] - p[j] for j in range(3)]
                n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
                for k in t:
                    normals[k] = [normals[k][j] + n[j] for j in range(3)]
            for n in normals:
                length = math.sqrt(sum(x * x for x in n)) or 1.0
                n[:] = [x / length for x in n]
            flat = array("f", [x for p in points for x in p])
            lo = [min(p[j] for p in points) for j in range(3)]
            hi = [max(p[j] for p in points) for j in range(3)]
            at = len(accessors)
            accessors.append({"bufferView": view(flat.tobytes(), 34962, 12), "componentType": 5126, "count": len(points), "type": "VEC3", "min": lo, "max": hi})
            accessors.append({"bufferView": view(array("f", [x for n in normals for x in n]).tobytes(), 34962, 12), "componentType": 5126, "count": len(points), "type": "VEC3"})
            accessors.append({"bufferView": view(array("I", [k for t in tris for k in t]).tobytes(), 34963), "componentType": 5125, "count": 3 * len(tris), "type": "SCALAR"})
            primitive = {"attributes": {"POSITION": at, "NORMAL": at + 1}, "indices": at + 2, "mode": 4}
            if rgba is not None:
                primitive["material"] = material(rgba)
            primitives.append(primitive)
        if not primitives:
            return None
        meshes.append({"name": name, "primitives": primitives})
        return len(meshes) - 1

    def node(tree, rgba, rows=None) -> int:
        name, shape, own_rows, own_rgba, children, _ = tree
        rgba = own_rgba if own_rgba is not None else rgba
        rows = rows if rows is not None else own_rows
        n = {"name": name}
        translation = [r[3] for r in rows]
        if any(abs(x) > 1e-12 for x in translation):
            n["translation"] = translation
        if any(abs(v - (1.0 if i == j else 0.0)) > 1e-12 for i, r in enumerate(rows) for j, v in enumerate(r[:3])):
            n["rotation"] = _quaternion(rows)
        nodes.append(n)
        index = len(nodes) - 1
        if shape is not None:
            m = mesh(name, shape, rgba)
            if m is not None:
                n["mesh"] = m
        kids = [node(child, rgba) for child in children]
        if kids:
            n["children"] = kids
        return index

    tree = _tree(assy)
    # (The root's placement followed by a turn about X by -90°, Z up to Y up, as CadQuery's
    # `assy.loc *= Location((0, 0, 0), (1, 0, 0), -90)`.)
    turn = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]
    own = tree[2]
    root = [[sum(own[i][k] * turn[k][j] for k in range(3)) for j in range(3)] + [own[i][3]] for i in range(3)]
    node(tree, None, root)
    while len(data) % 4:
        data.append(0)
    doc = {"asset": {"generator": "pykentron", "version": "2.0"}, "scene": 0, "scenes": [{"nodes": [0]}], "nodes": nodes, "meshes": meshes, "accessors": accessors, "bufferViews": views}
    if materials:
        doc["materials"] = materials
    if binary:
        doc["buffers"] = [{"byteLength": len(data)}]
        text = json.dumps(doc, separators=(",", ":")).encode()
        text += b" " * (-len(text) % 4)
        total = 12 + 8 + len(text) + 8 + len(data)
        with open(path, "wb") as f:
            f.write(struct.pack("<4sII", b"glTF", 2, total))
            f.write(struct.pack("<I4s", len(text), b"JSON") + text)
            f.write(struct.pack("<I4s", len(data), b"BIN\x00") + bytes(data))
    else:
        bin_path = os.path.splitext(str(path))[0] + ".bin"
        doc["buffers"] = [{"byteLength": len(data), "uri": os.path.basename(bin_path)}]
        with open(bin_path, "wb") as f:
            f.write(bytes(data))
        with open(path, "w") as f:
            json.dump(doc, f, indent=1)
    return True


def exportVTKJS(assy, path: str) -> None:
    raise NotImplementedError("exporting assemblies to VTK.js files")

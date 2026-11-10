#!/usr/bin/env python3
"""The drawings of docs/known-gaps.md: docs/known-gaps/*.svg.

    gap_figures.py [NAME ...] [--rebuild] [--list]

Each figure shows one entry's problem: the shape kentron builds beside the one expected
(CadQuery's or OCCT's, or one built from an exact construction where OCCT fails too), or a
diagram where the problem is not a shape (tolerance and copies are sketches, and say so).
The shapes are built by the entry's own repro: a CadQuery script run on both kernels
(kentron's through pykentron's shim, from this source tree: build it first with `maturin
develop`), or a DRAW test, run on kentron by crates/occt-tests' draw_commands and on OCCT
through OCP. They are drawn by the small renderer below as plain vector SVG: shaded faces,
the faces' edges and the silhouettes, with the wrong region marked. Nothing is embedded as
a bitmap and the colours follow the reader's light or dark theme, so the files show as they
are on GitHub.

Run with a Python that has CadQuery, numpy, scipy, contourpy and shapely, and with cargo on
the path. The wire figure reads offset_wire_058.brep of OCCT's test data, looked for under
$OCCT_TEST_DATA and ~/occt-data. Meshes are cached in target/gap-figures; --rebuild drops
them.
"""

import math
import os
import pickle
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "docs" / "known-gaps"
CACHE = REPO / "target" / "gap-figures"
PYKENTRON = REPO / "crates" / "pykentron" / "python"

# Materials of faces: the fill classes of the stylesheet below.
BODY, WRONG, FIRST, SECOND, THIRD, FOURTH = range(6)
MATERIAL_CLASSES = ["body", "wrong", "first", "second", "third", "fourth"]


# ---------------------------------------------------------------------------------------
# Scenes: shapes built through the CadQuery API, on either kernel (in a worker process),
# or by DRAW commands on kentron and through OCP on OCCT.
# A scene returns a dict: name -> list of (shape, material), where material is one of the
# constants above or a function of a face giving one.
# ---------------------------------------------------------------------------------------

def scene_saddle(cq, kernel):
    arms = cq.Workplane().box(15, 6, 6).translate((7.5, 0, 0))
    arms = arms.union(cq.Workplane().box(6, 15, 6).translate((0, 7.5, 0)))
    arms = arms.union(cq.Workplane().box(6, 6, 15).translate((0, 0, 7.5)))
    return {"arms": [(arms.val(), BODY)], "filleted": [(arms.edges().fillet(1).val(), BODY)]}


def scene_sliver(cq, kernel):
    box = cq.Workplane().box(20, 20, 20).edges().chamfer(2)
    return {"shelled": [(box.faces(">Z").shell(-3).val(), BODY)]}


def scene_upto(cq, kernel):
    """After example24 of CadQuery's documentation, as far as its cutBlind."""
    sphere = cq.Workplane().sphere(5)
    base = cq.Workplane(origin=(0, 0, -2)).box(12, 12, 10).cut(sphere).edges("|Z").fillet(2)
    sphere_face = base.faces(">>X[2] and (not |Z) and (not |Y)").val()
    base = base.faces("<Z").workplane().circle(2).extrude(10)
    shaft = cq.Workplane().sphere(4.5).circle(1.5).extrude(20)
    joint = base.union(shaft)
    cut = joint.faces(">X").workplane(centerOption="CenterOfMass").move(0, 4).slot2D(10, 2, 90).cutBlind(sphere_face)
    return {"joint": [(joint.val(), BODY)], "cut": [(cut.val(), BODY)]}


def scene_rounded(cq, kernel):
    """The box of usd-brep's TestSmRemoveFilletNurbsPlanes: 10 a side, every edge rounded by 1."""
    box = cq.Workplane().box(10, 10, 10)
    return {"rounded": [(box.edges().fillet(1).val(), lambda face: BODY if face.geomType() == "PLANE" else FIRST)], "box": [(box.val(), BODY)]}


def scene_twisted(cq, kernel):
    """After knurl() of mryndzionek's cadquery_plaground (util.py): height 10, radius 26.75,
    cut angle 120, depth 0.5, twist 40, 12 prisms a set."""
    depth, radius = 0.5, 26.75
    x1, y1 = 2 * depth * math.cos(math.radians(60)), 2 * depth * math.sin(math.radians(60))

    def prism(start, twist):
        at = cq.Location(cq.Vector(radius * math.cos(math.radians(start)), radius * math.sin(math.radians(start)), 0), cq.Vector(0, 0, 1), start)
        sketch = cq.Sketch().push([(-depth, 0)]).segment((0, 0), (x1, y1)).segment((x1, -y1)).close().assemble().reset().moved(at)
        return cq.Workplane().placeSketch(sketch).twistExtrude(10, twist).val()

    a, b = prism(-60, 40), prism(0, -40)
    out = {"a": [(a, FIRST)], "b": [(b, SECOND)],
           "sets": [(prism(30 * i, 40), FIRST) for i in range(12)] + [(prism(30 * i, -40), SECOND) for i in range(12)]}
    for name, operation in (("union", a.fuse), ("common", a.intersect)):
        try:
            out[name] = [(operation(b), BODY)]
        except Exception as error:
            out[name] = str(error)
    return out


def scene_gyroid(cq, kernel):
    """After LatticeQuery's gyroid cell (lq/topologies/gyroid.py): four of its eight plates."""
    T, U, h, d = 0.1, 10.0, 5.0, 1e-8

    def plate(thickness):
        ends = [((0.5, 0.5), (0.25, 0.0), (0.5, -0.5)), ((-0.5, -0.5), (0.0, -0.25), (0.5, -0.5)), ((-0.5, -0.5), (0.0, -0.25), (0.5, -0.5)),
                ((-0.5, -0.5), (-0.25, 0.0), (-0.5, 0.5)), ((0.5, 0.5), (0.0, 0.25), (-0.5, 0.5)), ((0.5, 0.5), (0.0, 0.25), (-0.5, 0.5))]
        wire = None
        for points, plane, offset in zip(ends, ["XZ", "XY", "YZ", "XZ", "YZ", "XY"], [-1, 1, 1, 1, -1, -1]):
            edge = cq.Workplane(plane).workplane(offset=-offset * h).spline([(x * U, y * U) for x, y in points])
            wire = edge if wire is None else wire.add(edge)
        face = cq.Workplane("XY").interpPlate(wire, [(0, 0, 0)], 0).val()
        both = cq.Workplane("XY").add(face.thicken(0.5 * thickness)).union(cq.Workplane("XY").add(face.thicken(-0.5 * thickness)))
        return cq.Workplane("XY").add(both.val().translate(cq.Vector(h, h, h)))

    g000 = plate(T)
    g100 = g000.mirror(mirrorPlane="ZY", basePointVector=(U, 0, 0)).mirror(mirrorPlane="XZ", basePointVector=(0, h, 0))
    g110 = plate(-T).mirror(mirrorPlane="XZ", basePointVector=(0, U, 0)).translate((U, 0, 0))
    g001 = g110.translate((-U, -U, U))
    g011 = g100.translate((-(1 + d) * U, (1 + d) * U, (1 + d) * U))
    plates = [g000, g110, g001, g011]
    out = {"plates": [(g.val(), material) for g, material in zip(plates, (FIRST, SECOND, THIRD, FOURTH))]}
    union = g000
    for g in plates[1:]:
        union = union.union(g)
    out["union"] = [(union.val(), BODY)]
    return out


def scene_mirrored(cq, kernel):
    """After LatticeQuery's Schwarz-P cell (lq/topologies/schwartz.py): an eighth of a cell
    of size 3.98, thickness 0.9, and its mirror image."""
    size, thickness = 1.99, 0.9
    r = 0.5 - 0.5 / math.sqrt(2)
    convex, concave = [(0.0, 0.5), (r, r), (0.5, 0.0)], [(-0.5, 0.0), (-r, -r), (0.0, -0.5)]
    wire = None
    for i, (plane, offset) in enumerate(zip(["XZ", "XY", "YZ", "XZ", "XY", "YZ"], [-1, -1, 1, 1, 1, -1])):
        a, m, b = [(x * size, y * size) for x, y in (convex if i % 2 == 0 else concave)]
        edge = cq.Workplane(plane).workplane(offset=-offset * size * 0.5).moveTo(*a).threePointArc(m, b)
        wire = edge if wire is None else wire.add(edge)
    face = cq.Workplane("XY").interpPlate(wire, [(0, 0, 0)], 0).val()
    plate = cq.Workplane("XY").add(face.thicken(0.5 * thickness)).union(cq.Workplane("XY").add(face.thicken(-0.5 * thickness)))
    half = cq.Workplane("XY").add(plate.val().translate(cq.Vector(0.995, 0.995, 0.995)))
    image = half.mirror("ZY", (1.99, 0, 0))
    return {"half": [(half.val(), FIRST)], "image": [(image.val(), SECOND)], "union": [(half.union(image).val(), BODY)]}


def scene_lamp(cq, kernel):
    """After jhnphm's spiral lamp (lamp.py), as far as its first union."""
    import numpy as np
    from scipy.special import erf

    width, depth, diffuser, lip, lip_depth = 26, 8, 1.6, 1, 1
    least_radius, radius_scale, z_scale, tolerance, loops = 30, 75, 75, 0.01, 3

    def spiral(t):
        e = -erf(t / (2 * loops * np.pi / (1.2 * 2)) - 1.2)
        radius = -e * radius_scale + least_radius + radius_scale
        return radius * np.cos(t), radius * np.sin(t), e * z_scale

    xs, ys, zs = spiral(0)
    section_width = lip * 2 + diffuser * 2 + width
    section_height = depth + lip_depth + tolerance
    edge = section_height / 2

    def channel():
        return (cq.Workplane("XZ", origin=(xs, ys, zs)).moveTo(-section_width / 2, section_height / 2).line(section_width, 0).line(0, -edge)
                .line(-diffuser, 0).line(0, -(section_height - edge)).line(-lip, 0).line(0, lip_depth).line(-width, 0).line(0, -lip_depth).line(-lip, 0))

    long = xs + lip + width / 2 + width
    arm = channel().line(0, section_height - edge).line(-diffuser, 0).line(0, edge).close().sweep(
        cq.Workplane("XY").parametricCurve(spiral, start=0, stop=loops * np.pi * 2).val())
    negative = (channel().line(0, lip_depth).line(-xs, 0).line(0, -depth).line(long + width, 0).line(0, depth * 3).line(-long - width, 0).line(0, -depth).close()
                .sweep(cq.Workplane("XY").parametricCurve(spiral, start=-np.pi / 2, stop=np.pi / 2)))
    support = cq.Workplane("YZ", origin=(width / 2, 0, 0)).center(0, zs).rect(width + lip * 2, depth * 1.8).extrude(xs + width * .5).val()
    hub = support.cut(negative.val())
    return {"arm": [(arm.val(), FIRST)], "hub": [(hub, SECOND)], "union": [(arm.union(hub).val(), BODY)], "start": (float(xs), float(ys), float(zs))}


def draw_commands(commands):
    """kentron's DRAW interpreter run on commands (crates/occt-tests' draw_commands): what
    they print, or a RuntimeError saying which command failed and why."""
    done = subprocess.run(["cargo", "run", "--quiet", "--release", "-p", "occt-tests", "--bin", "draw_commands"], cwd=REPO, input="\n".join(commands) + "\n", capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError(done.stderr.strip().splitlines()[-1])
    return done.stdout


def data_file(name):
    """A file of OCCT's test data, found by name as the DRAW tests find it."""
    for root in [*os.environ.get("OCCT_TEST_DATA", "").split(os.pathsep), Path.home() / "occt-data"]:
        for found in (Path(root).rglob(name) if root else ()):
            return found
    raise RuntimeError("OCCT test data file not found: %s (set OCCT_TEST_DATA)" % name)


def scene_evolved(cq, kernel):
    """OCCT's tests/evolved/voluved/AGT003: a profile in the plane x + y = 1 swept along a
    polygon, as a solid. kentron's through its DRAW commands, OCCT's through OCP."""
    spine = [(0, 0, 0), (0, 10, 0), (20, 10, 0), (20, 0, 0), (13, 0, 0), (13, 5, 0), (6, 5, 0), (6, 0, 0), (0, 0, 0)]
    profile = [(-1, 2, 0), (0, 1, 5), (0, 1, 10), (0.5, 0.5, 10), (0.5, 0.5, 0)]
    if kernel == "kentron":
        CACHE.mkdir(parents=True, exist_ok=True)
        draw_commands(["polyline p " + " ".join("%g" % x for p in spine for x in p), "polyline t " + " ".join("%g" % x for p in profile for x in p),
                       "evolved result -s p -p t -a -solid -v", "save result %s" % (CACHE / "evolved.brep")])
        result = cq.Shape.importBrep(str(CACHE / "evolved.brep"))
    else:
        from OCP.BRepOffsetAPI import BRepOffsetAPI_MakeEvolved
        from OCP.GeomAbs import GeomAbs_Arc

        maker = BRepOffsetAPI_MakeEvolved(cq.Wire.makePolygon(spine).wrapped, cq.Wire.makePolygon(profile).wrapped, GeomAbs_Arc, True, True, False, 1e-7, True, False)
        maker.Build()
        result = cq.Shape.cast(maker.Shape())
    return {"result": [(result, BODY)], "spine": {"points": spine}, "profile": {"points": profile}}


def scene_throat(cq, kernel):
    """OCCT's tests/bugs/modalg_7/bug28828_2: a cylinder fused to a box tilted by 20 degrees,
    the edge where they meet chamfered with a throat of 2."""
    if kernel == "kentron":
        CACHE.mkdir(parents=True, exist_ok=True)
        draw_commands(["pcylinder a1 10 50", "ttranslate a1 0 -5 0", "box a2 -20 -20 -20 40 40 40", "trotate a2 0 0 0 1 0 0 20", "bfuse b a1 a2", "savehistory hh",
                       "explode a1 f", "generated ee hh a1_1", "chamf_throat result b ee 2.", "save b %s" % (CACHE / "throat-fused.brep"), "save result %s" % (CACHE / "throat.brep")])
        fused, result = (cq.Shape.importBrep(str(CACHE / name)) for name in ("throat-fused.brep", "throat.brep"))
    else:
        from OCP.BRepFilletAPI import BRepFilletAPI_MakeChamfer
        from OCP.ChFiDS import ChFiDS_ConstThroatChamfer

        fused = cq.Solid.makeCylinder(10, 50).translate(cq.Vector(0, -5, 0)).fuse(cq.Solid.makeBox(40, 40, 40, cq.Vector(-20, -20, -20)).rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), 20))
        maker = BRepFilletAPI_MakeChamfer(fused.wrapped)
        maker.SetMode(ChFiDS_ConstThroatChamfer)
        for edge in fused.Edges():
            if edge.geomType() == "ELLIPSE":
                maker.Add(2.0, edge.wrapped)
        maker.Build()
        result = cq.Shape.cast(maker.Shape())
    return {"fused": [(fused, BODY)], "result": [(result, lambda face: FIRST if face.geomType() not in ("PLANE", "CYLINDER") else BODY)]}


def scene_revolved(cq, kernel):
    """OCCT's tests/offset/faces_type_i/F3: a quarter of an ellipse and a polyline revolved
    by 270 degrees, hollowed by 1 inwards, open at the face the revolution starts from.
    kentron's through its DRAW commands. The hollow expected is built on OCCT from the
    profile's offset, revolved whole and cut by the planes of the side faces (the open one
    where it is, the other moved in by 1): OCCT's own offset gives the source back."""
    out = {}
    if kernel == "kentron":
        CACHE.mkdir(parents=True, exist_ok=True)
        try:
            draw_commands(["ellipse w1 0 0 0 15 10", "mkedge w1 w1 0 %r" % (math.pi / 2), "trotate w1 0 0 0 1 0 0 90", "polyline w2 0 0 10 0 0 -10 25 0 -10 15 0 0", "wire w w1 w2", "mkplane w w",
                           "revol s w 0 0 0 0 0 1 270", "save s %s" % (CACHE / "revolved.brep"), "nearest open s F 8.78408 0 -1.83548 317.81", "offsetparameter 1e-7 c i", "offsetload s -1 open", "offsetperform result",
                           "save result %s" % (CACHE / "revolved-hollow.brep")])
            out["hollow"] = [(cq.Shape.importBrep(str(CACHE / "revolved-hollow.brep")), BODY)]
        except RuntimeError as failure:
            out["error"] = str(failure)
        source = cq.Shape.importBrep(str(CACHE / "revolved.brep"))
    else:
        arc = cq.Edge.makeEllipse(15, 10, cq.Vector(0, 0, 0), cq.Vector(0, -1, 0), cq.Vector(1, 0, 0), 0, 90)
        lines = [(0, 10), (0, -10), (25, -10), (15, 0)]
        profile = cq.Wire.assembleEdges([arc] + [cq.Edge.makeLine(cq.Vector(a, 0, b), cq.Vector(c, 0, d)) for (a, b), (c, d) in zip(lines, lines[1:])])
        source = cq.Solid.revolve(cq.Face.makeFromWires(profile), 270, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
        # Its offset by 1 inwards. The ellipse's is a spline through its points moved along
        # their normals, run on past the ellipse's end to where it meets the cone's (the
        # line x + z = 15 - sqrt 2); the axis stays.
        def moved(t):
            nx, nz = 10 * math.cos(t), 15 * math.sin(t)
            return 15 * math.cos(t) - nx / math.hypot(nx, nz), 10 * math.sin(t) - nz / math.hypot(nx, nz)

        low, high = -0.3, 0.0
        for _ in range(60):
            middle = (low + high) / 2
            low, high = (middle, high) if sum(moved(middle)) < 15 - math.sqrt(2) else (low, middle)
        points = [moved(low + (math.pi / 2 - low) * i / 400) for i in range(400)] + [(0.0, 9.0)]
        outline = cq.Workplane("XZ").moveTo(0, -9).lineTo(24 - math.sqrt(2), -9).lineTo(*points[0]).spline(points[1:], includeCurrent=True).close()
        inner = outline.revolve(360, (0, 0, 0), (0, 1, 0)).val()
        inner = inner.cut(cq.Solid.makeBox(60, 60, 60, cq.Vector(-1, -60, -30)))
        out["inner"] = [(inner, BODY)]
        out["hollow"] = [(source.cut(inner), BODY)]
    # The face left open: the plane y = 0.
    out["source"] = [(source, lambda face: FIRST if face.geomType() == "PLANE" and abs(face.Center().y) < 1e-6 and face.Center().x > 1 else BODY)]
    return out


def scene_wire(cq, kernel):
    """OCCT's tests/offset/wire_unclosed_outside_0_025/B6: an open wire offset by 2.5% of
    its length. kentron's contour through its DRAW commands, and the boundary of the points
    within the distance (tools/brute.py), both as points in the wire's plane."""
    import numpy as np
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps

    sys.path.insert(0, str(REPO / "tools"))
    import brute

    path = data_file("offset_wire_058.brep")
    wire = brute.load(path)
    props = GProp_GProps()
    BRepGProp.LinearProperties_s(wire, props)
    distance = 0.025 * props.Mass()
    CACHE.mkdir(parents=True, exist_ok=True)
    draw_commands(["restore %s s" % path, "mkoffset result s 1 %r" % -distance, "save result_1 %s" % (CACHE / "wire.brep")])
    points, _, (origin, axes) = brute.wire_points(wire, 20000)
    contour = cq.Shape.importBrep(str(CACHE / "wire.brep"))
    pieces = [(np.array([p.toTuple() for p in edge.positions(np.linspace(0, 1, 400))]) - origin) @ axes.T for edge in contour.Edges()]
    loops = brute.wire_offset_loops(wire, distance, 200000)
    # The turn at the wire's tip takes a hundredth of its curve's parameters: sampled closely.
    at_tip = points[:, 0].argmax() / (len(points) - 1)
    tip = (np.array([p.toTuple() for p in cq.Shape.cast(wire).Edges()[0].positions(np.linspace(at_tip - 0.02, at_tip + 0.02, 40000), mode="parameter")]) - origin) @ axes.T
    return {"wire": points, "tip": tip, "contour": tuple(pieces), "boundary": tuple(loops), "distance": distance, "length": props.Mass(),
            "contour_length": sum(edge.Length() for edge in contour.Edges())}


SCENES = {name[6:]: f for name, f in list(globals().items()) if name.startswith("scene_")}


def dump_scene(name, kernel):
    """Worker: builds a scene on one kernel and pickles its meshes."""
    if kernel == "kentron":
        sys.path[:0] = [str(PYKENTRON / "pykentron" / "cadquery_shim"), str(PYKENTRON)]
    import cadquery as cq
    import numpy as np

    out = {}
    for part, shapes in SCENES[name](cq, kernel).items():
        if not isinstance(shapes, list):
            out[part] = shapes
            continue
        faces, edges, info = [], [], []
        for shape, material in shapes:
            mesh_faces = shape.Faces()
            if not mesh_faces:
                info.append({"faces": 0, "solids": 0, "volume": 0.0, "area": 0.0, "valid": False})
                continue
            box = shape.BoundingBox()
            size = max(box.xlen, box.ylen, box.zlen)
            info.append({"faces": len(mesh_faces), "solids": len(shape.Solids()), "volume": shape.Volume(), "area": shape.Area(), "valid": shape.isValid()})
            first = len(faces)
            seen = {}
            for index, face in enumerate(mesh_faces):
                vertices, triangles = face.tessellate(size * 2e-4, 0.05)
                faces.append({
                    "v": np.array([(v.x, v.y, v.z) for v in vertices], dtype=float).reshape(-1, 3),
                    "t": np.array(triangles, dtype=np.int32).reshape(-1, 3),
                    "plane": face.geomType() == "PLANE",
                    "area": face.Area(),
                    "centre": face.Center().toTuple(),
                    "material": material(face) if callable(material) else material,
                })
                for edge in face.Edges():
                    if edge.Length() < size * 1e-7:
                        continue
                    a, m, b = edge.positions([0.31, 0.5, 0.77])
                    key = tuple(round(x, 5) for p in (a, m, b) for x in (p.x, p.y, p.z))
                    key = min(key, key[6:] + key[3:6] + key[:3])
                    try:
                        n = face.normalAt(m)
                        n = (n.x, n.y, n.z)
                    except Exception:
                        n = None
                    if key in seen:
                        seen[key]["sides"].append((first + index, n))
                        continue
                    count = 2 if edge.geomType() == "LINE" else 97
                    points = edge.positions([i / (count - 1) for i in range(count)])
                    seen[key] = {"p": np.array([(p.x, p.y, p.z) for p in points]), "sides": [(first + index, n)]}
            for edge in seen.values():
                sides = edge.pop("sides")
                edge["kind"] = "sharp"
                if len(sides) == 2:
                    (fa, na), (fb, nb) = sides
                    if fa == fb:
                        edge["kind"] = "seam"
                    elif na is not None and nb is not None and abs(sum(x * y for x, y in zip(na, nb))) > 0.975:
                        edge["kind"] = "smooth"
                edge["faces"] = [f for f, _ in sides]
                edges.append(edge)
        out[part] = {"faces": faces, "edges": edges, "info": info}
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(cached(name, kernel), "wb") as f:
        pickle.dump(out, f)


def cached(name, kernel):
    """Where a scene's meshes are kept: a file named after the code that builds them."""
    import hashlib
    import inspect

    code = inspect.getsource(SCENES[name]) + inspect.getsource(dump_scene)
    return CACHE / ("%s-%s-%s.pickle" % (name, kernel, hashlib.sha1(code.encode()).hexdigest()[:10]))


def scene(name, kernel, rebuild=False):
    """A scene's meshes on a kernel ('kentron' or 'cadquery'), built in a worker if not cached."""
    path = cached(name, kernel)
    if rebuild or not path.exists():
        done = subprocess.run([sys.executable, __file__, "--dump", name, kernel], capture_output=True, text=True)
        if done.returncode != 0:
            raise RuntimeError("%s on %s:\n%s" % (name, kernel, done.stderr[-3000:]))
    with open(path, "rb") as f:
        return pickle.load(f)


# ---------------------------------------------------------------------------------------
# The renderer: an orthographic view of meshes, hidden surfaces removed in a depth buffer,
# whose regions and shading bands are then traced back into vector outlines.
# ---------------------------------------------------------------------------------------

SUPERSAMPLING = 4
LIGHT = (-0.42, 0.62, 0.66)                          # In the view's frame: right, up, to the viewer.
BANDS = [0.30 + 0.08 * i for i in range(9)]          # Each band reached lightens the face a step.


class View:
    """An orthographic view: `eye` points from the scene to the viewer, `span` model units
    fill the width, `centre` is in the middle."""

    def __init__(self, eye, centre, span, width, height, up=(0, 0, 1)):
        import numpy as np

        self.width, self.height = width, height
        self.scale = width / span
        self.centre = np.array(centre, dtype=float)
        e = np.array(eye, dtype=float)
        self.e = e / np.linalg.norm(e)
        r = np.cross(np.array(up, dtype=float), self.e)
        self.r = r / np.linalg.norm(r)
        self.u = np.cross(self.e, self.r)

    def project(self, points):
        """Model points -> (x, y, depth), in pixels of the view; depth grows to the viewer."""
        import numpy as np

        p = np.atleast_2d(np.array(points, dtype=float)) - self.centre
        return np.stack([self.width / 2 + p @ self.r * self.scale, self.height / 2 - p @ self.u * self.scale, p @ self.e * self.scale], axis=1)

    def to_frame(self, vectors):
        import numpy as np

        return np.stack([vectors @ self.r, vectors @ self.u, vectors @ self.e], axis=1)

    def zoomed(self, centre, span, width=None, height=None):
        """The same view closer in."""
        return View(self.e, centre, span, width or self.width, height or self.height, self.u)

    def window(self, other):
        """Where another view of the same direction lies in this one: x, y, width, height."""
        x, y, _ = self.project(other.centre)[0]
        w, h = other.width / other.scale * self.scale, other.height / other.scale * self.scale
        return x - w / 2, y - h / 2, w, h


def vertex_normals(vertices, triangles):
    import numpy as np

    normals = np.zeros_like(vertices)
    a, b, c = (vertices[triangles[:, i]] for i in range(3))
    weighted = np.cross(b - a, c - a)
    for i in range(3):
        np.add.at(normals, triangles[:, i], weighted)
    length = np.linalg.norm(normals, axis=1)
    return normals / np.where(length > 0, length, 1)[:, None]


def simplify(points, tolerance):
    """Douglas-Peucker; a closed ring (first point repeated last) stays closed."""
    import numpy as np

    n = len(points)
    if n < 3:
        return points
    keep = np.zeros(n, dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    if np.allclose(points[0], points[-1]):
        far = int(np.argmax(np.hypot(*(points - points[0]).T)))
        keep[far] = True
        stack = [(0, far), (far, n - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        chord = points[b] - points[a]
        length = math.hypot(*chord)
        rel = points[a + 1:b] - points[a]
        d = np.hypot(rel[:, 0], rel[:, 1]) if length < 1e-12 else np.abs(chord[0] * rel[:, 1] - chord[1] * rel[:, 0]) / length
        i = int(np.argmax(d))
        if d[i] > tolerance:
            keep[a + 1 + i] = True
            stack += [(a, a + 1 + i), (a + 1 + i, b)]
    return points[keep]


def number(x):
    text = "%.1f" % x
    return text[:-2] if text.endswith(".0") else text


def path_data(points, closed=False):
    return "M" + "L".join("%s %s" % (number(x), number(y)) for x, y in points) + ("Z" if closed else "")


def render(meshes, view, recolour=None, hide=(), hide_smooth=False, outline=True, reference=None, tolerance=2.0):
    """SVG elements (a string) drawing the meshes in the view. `recolour` gives materials to
    faces by their index (counted through the meshes), over the scene's own; the faces whose
    indices are in `hide` are left out. With `reference` meshes, what is seen farther than
    the tolerance (in pixels) from every face of theirs is drawn as WRONG."""
    import numpy as np
    from contourpy import FillType, contour_generator
    from scipy import ndimage

    ss = SUPERSAMPLING
    W, H = int(view.width * ss), int(view.height * ss)
    depth = np.full((H, W), -np.inf, dtype=np.float64)
    owner = np.full((H, W), -1, dtype=np.int32)
    normal = np.zeros((H, W, 3), dtype=np.float32)
    faces, edges, silhouettes = [], [], []
    for mesh in meshes:
        edges += [e for e in mesh["edges"] if not all(len(faces) + f in hide for f in e["faces"])]
        faces += mesh["faces"]
    for index, face in enumerate(faces):
        if len(face["t"]) == 0 or index in hide:
            continue
        p = view.project(face["v"]) * [ss, ss, 1]
        n = view.to_frame(vertex_normals(face["v"], face["t"]))
        for i0, i1, i2 in face["t"]:
            (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) = p[i0], p[i1], p[i2]
            left, right = max(int(min(x0, x1, x2)), 0), min(int(max(x0, x1, x2)) + 1, W)
            top, bottom = max(int(min(y0, y1, y2)), 0), min(int(max(y0, y1, y2)) + 1, H)
            den = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
            if left >= right or top >= bottom or abs(den) < 1e-12:
                continue
            X = np.arange(left, right)[None, :] + 0.5 - x2
            Y = np.arange(top, bottom)[:, None] + 0.5 - y2
            w0 = ((y1 - y2) * X + (x2 - x1) * Y) / den
            w1 = ((y2 - y0) * X + (x0 - x2) * Y) / den
            w2 = 1 - w0 - w1
            z = w0 * z0 + w1 * z1 + w2 * z2
            window = (slice(top, bottom), slice(left, right))
            hit = (w0 >= -1e-9) & (w1 >= -1e-9) & (w2 >= -1e-9) & (z > depth[window])
            if not hit.any():
                continue
            depth[window][hit] = z[hit]
            owner[window][hit] = index
            normal[window][hit] = (w0[..., None] * n[i0] + w1[..., None] * n[i1] + w2[..., None] * n[i2])[hit]
        if not face["plane"]:
            silhouettes += silhouette_chains(face, n[:, 2])

    inside = owner >= 0
    length = np.linalg.norm(normal, axis=2)
    unit = normal / np.where(length > 0, length, 1)[..., None]
    unit = np.where(unit[..., 2:3] < 0, -unit, unit)      # What is seen faces the viewer.
    light = np.array(LIGHT) / np.linalg.norm(LIGHT)
    shade = np.where(inside, 0.22 + 0.78 * np.clip(unit @ light, 0, 1), -1.0)

    def regions(field, level, cls, tolerance=0.35):
        polygons, offsets = contour_generator(z=field, fill_type=FillType.OuterOffset).filled(level, 1e9)
        data = []
        for polygon, offset in zip(polygons, offsets):
            for a, b in zip(offset[:-1], offset[1:]):
                ring = simplify((polygon[a:b] + 0.5) / ss, tolerance / ss * 2)
                if len(ring) > 3:
                    data.append(path_data(ring[:-1], closed=True))
        return '<path class="%s" fill-rule="evenodd" d="%s"/>\n' % (cls, "".join(data)) if data else ""

    out = []
    face_material = np.array([f["material"] for f in faces] + [BODY])
    for index, material in (recolour or {}).items():
        face_material[index] = material
    material = np.where(inside, face_material[owner], -1)
    if reference is not None:
        from scipy.spatial import cKDTree

        rows, columns = np.nonzero(inside)
        seen = (view.centre + np.outer(((columns + 0.5) / ss - view.width / 2) / view.scale, view.r)
                + np.outer((view.height / 2 - (rows + 0.5) / ss) / view.scale, view.u) + np.outer(depth[rows, columns] / view.scale, view.e))
        reach = tolerance / view.scale
        far = cKDTree(sampled(reference, 0.9 * reach, view)).query(seen, distance_upper_bound=2 * reach, workers=-1)[0] > reach
        material[rows[far], columns[far]] = WRONG
    # Slightly swollen, so that no background shows between a region and the next.
    out.append(regions(ndimage.gaussian_filter(inside.astype(np.float32), 0.8), 0.4, "body"))
    for level in BANDS:
        out.append(regions(shade, level, "lit"))
    for m in sorted(set(np.unique(material).tolist()) - {BODY, -1}):
        out.append(regions(ndimage.gaussian_filter((material == m).astype(np.float32), 0.6), 0.5, "tint " + MATERIAL_CLASSES[m]))

    # Lines: the faces' edges and the silhouettes of curved faces, where nothing hides them.
    behind = ndimage.minimum_filter(depth, size=2 * ss + 1)
    slack = 0.004 * view.width

    def strokes(points, cls):
        p = view.project(points)
        steps = np.hypot(*np.diff(p[:, :2], axis=0).T)
        along = np.concatenate([[0], np.cumsum(steps)])
        if along[-1] < 0.5:
            return ""
        t = np.linspace(0, along[-1], int(along[-1] / 0.5) + 2)
        x, y, z = (np.interp(t, along, p[:, i]) for i in range(3))
        column, row = (x * ss).astype(int), (y * ss).astype(int)
        shown = (column >= 0) & (column < W) & (row >= 0) & (row < H)
        shown[shown] = z[shown] >= behind[row[shown], column[shown]] - slack
        data, start = [], None
        for i in range(len(t) + 1):
            if i < len(t) and shown[i]:
                start = i if start is None else start
            elif start is not None:
                if i - start > 2:
                    data.append(path_data(simplify(np.stack([x[start:i], y[start:i]], axis=1), 0.07)))
                start = None
        return '<path class="%s" d="%s"/>\n' % (cls, "".join(data)) if data else ""

    lines = {"smooth": [], "sharp": [], "silhouette": []}
    for edge in edges:
        if edge["kind"] == "seam" or (edge["kind"] == "smooth" and hide_smooth):
            continue
        lines[edge["kind"]].append(strokes(edge["p"], edge["kind"]))
    for chain in silhouettes:
        lines["silhouette"].append(strokes(chain, "sharp"))
    out += lines["smooth"] + lines["silhouette"] + lines["sharp"]
    if outline:
        polygons = contour_generator(z=inside.astype(np.float32)).lines(0.5)
        data = "".join(path_data(simplify((p + 0.5) / ss, 0.12), closed=bool(np.allclose(p[0], p[-1]))) for p in polygons if len(p) > 8)
        out.append('<path class="outline" d="%s"/>\n' % data)
    return "".join(out)


def silhouette_chains(face, facing):
    """Where a curved face turns away from the viewer: the curves along which its normal
    (interpolated over the mesh) is square to the line of sight, as chains of points."""
    import numpy as np

    facing = np.where(facing == 0, 1e-12, facing)
    points, links = {}, {}
    for triangle in face["t"]:
        crossing = []
        for a, b in ((triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0])):
            if facing[a] * facing[b] < 0:
                key = (min(a, b), max(a, b))
                t = facing[a] / (facing[a] - facing[b])
                points[key] = face["v"][a] + t * (face["v"][b] - face["v"][a])
                crossing.append(key)
        if len(crossing) == 2:
            links.setdefault(crossing[0], []).append(crossing[1])
            links.setdefault(crossing[1], []).append(crossing[0])
    chains, used = [], set()
    for start in sorted(links, key=lambda k: len(links[k])):     # Ends of open chains first.
        if start in used:
            continue
        chain, at = [start], start
        used.add(start)
        while True:
            following = [k for k in links[at] if k not in used]
            if not following:
                break
            at = following[0]
            used.add(at)
            chain.append(at)
        if len(links[start]) == 2 and start in links[at] and len(chain) > 2:
            chain.append(start)
        if len(chain) > 1:
            chains.append(np.array([points[k] for k in chain]))
    return chains


def sampled(meshes, spacing, view=None):
    """Points on the meshes' faces, no farther apart than the spacing; with a view, of the
    triangles it shows only."""
    import numpy as np

    out = [np.zeros((1, 3))]
    triangles = np.concatenate([face["v"][face["t"]] for mesh in meshes for face in mesh["faces"]] + [np.zeros((0, 3, 3))])
    while len(triangles):
        if view is not None:
            p = view.project(triangles.reshape(-1, 3))[:, :2].reshape(-1, 3, 2)
            low, high = p.min(axis=1), p.max(axis=1)
            triangles = triangles[(high[:, 0] > -20) & (low[:, 0] < view.width + 20) & (high[:, 1] > -20) & (low[:, 1] < view.height + 20)]
        # Cut in two across the longest side until none is longer than the spacing: the
        # triangles' middles are then near enough to all of their points.
        sides = np.linalg.norm(np.roll(triangles, -1, axis=1) - triangles, axis=2)
        long = sides.max(axis=1) > spacing
        out.append(triangles[~long].mean(axis=1))
        first = sides[long].argmax(axis=1)
        a, b, c = (triangles[long][np.arange(len(first)), (first + i) % 3] for i in range(3))
        middle = (a + b) / 2
        triangles = np.concatenate([np.stack([a, middle, c], axis=1), np.stack([middle, b, c], axis=1)])
    return np.concatenate(out)


def traced(view, points, cls, closed=False):
    """A path through model points as the view shows them, hidden by nothing."""
    return '<path class="%s" d="%s"/>\n' % (cls, path_data(view.project(points)[:, :2], closed))


def marked(view, mesh, uses=lambda count: count != 2):
    """A mesh's edges bounding other than two faces, drawn through whatever hides them."""
    return "".join(traced(view, edge["p"], "mark") for edge in mesh["edges"] if uses(len(edge["faces"])))


def sliced(faces, axis, value):
    """Where the plane of one coordinate (`axis`) at a value cuts faces' meshes: segments,
    an array of (count, 2, 3)."""
    import numpy as np

    out = [np.zeros((0, 2, 3))]
    for face in faces:
        t = face["v"][face["t"]]
        d = t[:, :, axis] - value
        for i in range(3):
            a, b, c = d[:, i], d[:, (i + 1) % 3], d[:, (i + 2) % 3]
            alone = (a * b < 0) & (a * c < 0)
            pa, pb, pc = t[alone, i], t[alone, (i + 1) % 3], t[alone, (i + 2) % 3]
            tb, tc = (a[alone] / (a[alone] - b[alone]))[:, None], (a[alone] / (a[alone] - c[alone]))[:, None]
            out.append(np.stack([pa + tb * (pb - pa), pa + tc * (pc - pa)], axis=1))
    return np.concatenate(out)


def unmatched(mesh, reference, distance=1e-3, relative=1e-3):
    """The indices of a mesh's faces with no face of the same area and centre in a reference."""
    out = []
    for index, face in enumerate(mesh["faces"]):
        for other in reference["faces"]:
            if abs(face["area"] - other["area"]) <= relative * max(other["area"], 1e-9) and math.dist(face["centre"], other["centre"]) <= distance:
                break
        else:
            out.append(index)
    return out


# ---------------------------------------------------------------------------------------
# Figures: a card with panels, labels and callouts.
# ---------------------------------------------------------------------------------------

STYLE = """
.card{fill:#fcfcfb;stroke:#e3e2dd}
.panel{fill:#f3f2ee}
text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans",Helvetica,Arial,sans-serif;font-size:13px;fill:#0b0b0b}
.title{font-size:14px;font-weight:600}
.note{fill:#52514e}
.small{font-size:12px}
.body{fill:#8f8e86}
.lit{fill:#fff;fill-opacity:.13}
.tint{fill-opacity:.6}
.wrong{fill:#d03b3b}
.first{fill:#2a78d6}
.second{fill:#eb6834}
.third{fill:#1baf7a}
.fourth{fill:#eda100}
.sharp,.smooth,.outline,.leader,.rule,.mark,.frame,.glyph,.ink,.thin,.dash,.sketch,.ghost,.one,.two{fill:none;stroke-linecap:round;stroke-linejoin:round}
.sharp{stroke:#1d1d1b;stroke-width:1}
.smooth{stroke:#1d1d1b;stroke-width:.6;stroke-opacity:.38}
.outline{stroke:#1d1d1b;stroke-width:1.7}
.leader{stroke:#52514e;stroke-width:.9}
.rule{stroke:#e3e2dd}
.mark{stroke:#d03b3b;stroke-width:1.8}
.frame{stroke:#0b0b0b;stroke-width:1.2}
.ink{stroke:#1d1d1b;stroke-width:1.5}
.thin{stroke:#52514e;stroke-width:.9}
.dash{stroke:#52514e;stroke-width:.9;stroke-dasharray:4 3}
.sketch,.one{stroke:#2a78d6;stroke-width:1.6}
.two{stroke:#eb6834;stroke-width:1.6}
.wash{fill-opacity:.22}
.cut{fill:#8f8e86;fill-opacity:.4}
.node{fill:#fcfcfb;stroke:#52514e}
.group{fill:#e9e8e2;stroke:#52514e}
.off{stroke:#d03b3b;stroke-width:1.8}
.halo{paint-order:stroke;stroke:#f3f2ee;stroke-opacity:.85;stroke-width:3.5px;stroke-linejoin:round}
.ghost{stroke:#8f8e86;stroke-width:1.1;stroke-dasharray:5 4}
.dot{fill:#52514e}
.bad{fill:#d03b3b}
.good{fill:#0c8f0c}
.glyph{stroke:#fff;stroke-width:1.7}
@media (prefers-color-scheme:dark){
.card{fill:#1a1a19;stroke:#383835}
.panel{fill:#232321}
text{fill:#f0efec}
.note{fill:#c3c2b7}
.body{fill:#4b4a45}
.lit{fill-opacity:.07}
.tint{fill-opacity:.55}
.wrong{fill:#e04a4a}
.first{fill:#3987e5}
.second{fill:#d95926}
.third{fill:#199e70}
.fourth{fill:#c98500}
.sharp,.smooth,.outline,.ink{stroke:#f0efec}
.leader,.thin,.dash{stroke:#c3c2b7}
.rule{stroke:#383835}
.mark{stroke:#ff6b6b}
.frame{stroke:#f0efec}
.sketch,.one{stroke:#3987e5}
.two{stroke:#d95926}
.halo{stroke:#232321}
.cut{fill:#4b4a45;fill-opacity:.9}
.node{fill:#1a1a19;stroke:#c3c2b7}
.group{fill:#2e2e2b;stroke:#c3c2b7}
.off{stroke:#ff6b6b}
.ghost{stroke:#77766f}
.dot{fill:#c3c2b7}
.bad{fill:#e04a4a}
.good{fill:#22a722}
.glyph{stroke:#1a1a19}
}
"""


def escape(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Figure:
    def __init__(self, width, height, summary):
        self.width, self.height, self.summary = width, height, summary
        self.parts, self.clips = [], 0

    def add(self, svg):
        self.parts.append(svg)

    def text(self, x, y, text, cls="", anchor="start"):
        """Text at (x, y); lines after the first go below."""
        for i, line in enumerate(text.split("\n")):
            self.add('<text x="%s" y="%s"%s%s>%s</text>\n' % (number(x), number(y + 16 * i), ' class="%s"' % cls if cls else "", ' text-anchor="%s"' % anchor if anchor != "start" else "", escape(line)))

    def panel(self, x, y, width, height, drawing="", title=None, note=None, verdict=None):
        """A framed drawing, a title and a note above it, and a verdict below:
        ('bad' | 'good', text), shown with a cross or a tick."""
        self.clips += 1
        self.add('<clipPath id="c%d"><rect x="%s" y="%s" width="%s" height="%s" rx="6"/></clipPath>\n' % (self.clips, number(x), number(y), number(width), number(height)))
        self.add('<rect class="panel" x="%s" y="%s" width="%s" height="%s" rx="6"/>\n' % (number(x), number(y), number(width), number(height)))
        self.add('<g clip-path="url(#c%d)"><g transform="translate(%s %s)">\n%s</g></g>\n' % (self.clips, number(x), number(y), drawing))
        if title:
            self.text(x + 2, y - 27 if note else y - 10, title, "title")
        if note:
            self.text(x + 2, y - 9, note, "note small")
        if verdict:
            self.verdict(x + 2, y + height + 20, *verdict)

    def verdict(self, x, y, kind, text):
        glyph = "M-2.7 -2.7L2.7 2.7M2.7 -2.7L-2.7 2.7" if kind == "bad" else "M-3.2 .2L-1 2.5L3.2 -2.4"
        self.add('<g transform="translate(%s %s)"><circle class="%s" r="7"/><path class="glyph" d="%s"/></g>\n' % (number(x + 7), number(y - 4.5), kind, glyph))
        self.text(x + 20, y, text, "small")

    def callout(self, at, text, to, anchor="start", cls="note small halo"):
        """A label at `to` with a leader to the point `at`."""
        (ax, ay), (tx, ty) = at, to
        lines = text.split("\n")
        ex = tx - 5 if anchor == "start" else tx + 5 if anchor == "end" else tx
        ey = ty - 4 if anchor != "middle" else (ty - 15 if ay < ty else ty + 6 + 16 * (len(lines) - 1))
        self.add('<path class="leader" d="M%s %sL%s %s"/><circle class="dot" cx="%s" cy="%s" r="2.2"/>\n' % (number(ex), number(ey), number(ax), number(ay), number(ax), number(ay)))
        self.text(tx, ty, text, cls, anchor)

    def arrow(self, a, b, cls="thin"):
        """A line from a to b with a head at b."""
        (ax, ay), (bx, by) = a, b
        length = math.hypot(bx - ax, by - ay) or 1
        dx, dy = (bx - ax) / length, (by - ay) / length
        head = [(bx, by), (bx - 8 * dx + 3.2 * dy, by - 8 * dy - 3.2 * dx), (bx - 8 * dx - 3.2 * dy, by - 8 * dy + 3.2 * dx)]
        self.add('<path class="%s" d="%s"/><path class="dot" d="%s"/>\n' % (cls, path_data([a, (bx - 6 * dx, by - 6 * dy)]), path_data(head, closed=True)))

    def ring(self, at, radius=13):
        self.add('<circle class="mark" cx="%s" cy="%s" r="%s"/>\n' % (number(at[0]), number(at[1]), number(radius)))

    def scale(self, x, y, length, label):
        """A scale bar from (x, y) to the right, its length named above it."""
        self.add('<path class="thin" d="M%s %sv5h%sv-5"/>\n' % (number(x), number(y - 5), number(length)))
        self.text(x + length / 2, y - 9, label, "note small", "middle")

    def node(self, x, y, width, height, text, cls="node"):
        """A box of a diagram with its text in the middle."""
        self.add('<rect class="%s" x="%s" y="%s" width="%s" height="%s" rx="5"/>\n' % (cls, number(x), number(y), number(width), number(height)))
        self.text(x + width / 2, y + height / 2 + 4 - 8 * text.count("\n"), text, "small", "middle")

    def frame(self, x, y, width, height):
        self.add('<rect class="frame" x="%s" y="%s" width="%s" height="%s" rx="2"/>\n' % (number(x), number(y), number(width), number(height)))

    def svg(self):
        return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" role="img" aria-label="%s">\n<style>%s</style>\n'
                '<rect class="card" x=".5" y=".5" width="%d" height="%d" rx="10"/>\n%s</svg>\n'
                % (self.width, self.height, self.width, self.height, escape(self.summary).replace('"', "&quot;"), STYLE.replace("\n", ""), self.width - 1, self.height - 1, "".join(self.parts)))


class Chart:
    """A flat drawing for a panel: (x, y) to the panel's pixels, y upwards, x drawn
    `stretch` times larger than y."""

    def __init__(self, width, height, centre, scale, stretch=1.0):
        self.width, self.height, self.centre, self.scale, self.stretch = width, height, centre, scale, stretch

    def at(self, x, y, origin=(0, 0)):
        return (origin[0] + self.width / 2 + (x - self.centre[0]) * self.scale * self.stretch, origin[1] + self.height / 2 - (y - self.centre[1]) * self.scale)

    def path(self, points, cls, closed=False):
        """A path through points. Of an open one, only the stretches by the drawing are kept,
        thinned to what shows."""
        import numpy as np

        pixels = np.array([self.at(x, y) for x, y in points])
        if closed:
            return '<path class="%s" d="%s"/>\n' % (cls, path_data(pixels, True))
        near = np.concatenate([[False], (np.abs(pixels[:, 0] - self.width / 2) < self.width / 2 + 30) & (np.abs(pixels[:, 1] - self.height / 2) < self.height / 2 + 30), [False]])
        keep = near[1:-1] | near[:-2] | near[2:]
        starts = [i for i in range(len(keep)) if keep[i] and (i == 0 or not keep[i - 1])]
        ends = [i + 1 for i in range(len(keep)) if keep[i] and (i + 1 == len(keep) or not keep[i + 1])]
        data = "".join(path_data(simplify(pixels[a:b], 0.12)) for a, b in zip(starts, ends) if b - a > 1)
        return '<path class="%s" d="%s"/>\n' % (cls, data) if data else ""


def signed(text):
    """With true minus signs."""
    return text.replace("-", "\u2212")


def at(view, x, y, point):
    """A model point's place in a figure, for a view drawn in a panel at (x, y)."""
    px, py, _ = view.project(point)[0]
    return x + px, y + py


# The usual layout: what is built on the left, then kentron's result and the expected one.
WIDTH, MARGIN, GAP, TOP = 860, 18, 14, 58
SIDE, MAIN = 216, 290
COLUMNS = (MARGIN, MARGIN + SIDE + GAP, MARGIN + SIDE + GAP + MAIN + GAP)


def figure_saddle(rebuild):
    kentron, cadquery = (scene("saddle", k, rebuild) for k in ("kentron", "cadquery"))
    height = 300
    figure = Figure(WIDTH, TOP + height + 40, "Filleting every edge where three box arms meet: kentron's fillets run on past the mixed corners, CadQuery blends them with corner patches.")
    eye = (1, 0.62, 0.52)
    whole = View(eye, (4.2, 4.2, 4.6), 27, SIDE, height)
    close = whole.zoomed((2.6, 2.6, 3.0), 11.5, MAIN, height)
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([kentron["arms"]], whole), "Three arms", "every edge then filleted by 1")
    figure.frame(*(a + b for a, b in zip(whole.window(close), (COLUMNS[0], TOP, 0, 0))))
    patches = {i: FIRST for i, f in enumerate(cadquery["filleted"]["faces"]) if not f["plane"] and 3.3 < f["area"] < 4.8}
    for x, mesh, title, recolour, verdict in (
            (COLUMNS[1], kentron, "kentron", None, ("bad", "area 1003.0: the corners are mitred")),
            (COLUMNS[2], cadquery, "CadQuery", patches, ("good", "area 966.2: corner patches blend the fillets"))):
        info = mesh["filleted"]["info"][0]
        figure.panel(x, TOP, MAIN, height, render([mesh["filleted"]], close, recolour), title, "%d faces, volume %.3f" % (info["faces"], info["volume"]), verdict)
    for corner in ((3, -3, 3), (-3, 3, 3)):
        figure.ring(at(close, COLUMNS[1], TOP, corner), 26)
    return figure


def figure_sliver(rebuild):
    import numpy as np

    kentron, cadquery = (scene("sliver", k, rebuild)["shelled"] for k in ("kentron", "cadquery"))
    height = 300
    figure = Figure(WIDTH, TOP + height + 40, "A chamfered box shelled inwards: kentron's cavity has three sliver faces at each of its corners, CadQuery's has none.")
    above = View((0.45, 0.3, 1), (0, 0, -1.5), 31, SIDE, height)
    corner = View((-1, -0.75, -0.6), (-7, -7, -7), 0.8, MAIN, height)
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([kentron], above), "Chamfered box, shelled", "walls 3 thick, open at the top")
    figure.callout(at(above, COLUMNS[0], TOP, (-7, -7, -7)), "a corner of the cavity: seen from\ninside the wall in the next panels", (COLUMNS[0] + 10, TOP + height - 28))
    for x, mesh, title in ((COLUMNS[1], kentron, "kentron"), (COLUMNS[2], cadquery, "CadQuery")):
        faces, info = mesh["faces"], mesh["info"][0]
        walls = {i for i, f in enumerate(faces) if np.abs(f["v"]).max() > 7.06}
        slivers = {i: WRONG for i, f in enumerate(faces) if f["area"] < 0.01}
        verdict = ("bad", "%d sliver faces, three at each corner of the cavity" % len(slivers)) if slivers else ("good", "the cavity's faces alone meet at the corner")
        figure.panel(x, TOP, MAIN, height, render([mesh], corner, slivers, walls), title, "%d faces, volume %.3f" % (info["faces"], info["volume"]), verdict)
    return figure


def figure_upto(rebuild):
    kentron, cadquery = (scene("upto", k, rebuild) for k in ("kentron", "cadquery"))
    height = 300
    figure = Figure(WIDTH, TOP + height + 40, "A slot cut up to the spherical socket of a ball joint: kentron's cut runs on through the shaft and the far wall, where the slot's rays miss the socket; CadQuery's stops at the socket.")
    whole = View((0.6, -1, 0.75), (0, 0, 1.5), 41, SIDE, height)
    close = whole.zoomed((0, 0, 0.6), 21, MAIN, height)
    turn = [math.pi * i / 16 for i in range(17)]
    slot = [(6, math.cos(a), 6 + math.sin(a)) for a in turn] + [(6, -math.cos(a), -2 - math.sin(a)) for a in turn]
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([kentron["joint"]], whole) + traced(whole, slot, "sketch", closed=True), "Ball joint", "a slot is cut up to the socket")
    figure.callout(at(whole, COLUMNS[0], TOP, (6, 1, 5)), "the slot's sketch,\ncut inwards", (COLUMNS[0] + SIDE - 8, TOP + 34), "end")
    solid = kentron["joint"]["info"][0]["volume"]
    for x, mesh, reference, title, verdict in (
            (COLUMNS[1], kentron, [cadquery["cut"]], "kentron", ("bad", "removes %.2f: cut on past the socket")),
            (COLUMNS[2], cadquery, None, "CadQuery", ("good", "removes %.2f: cut only up to the socket"))):
        info = mesh["cut"]["info"][0]
        figure.panel(x, TOP, MAIN, height, render([mesh["cut"]], close, reference=reference), title, "%d faces, volume %.1f" % (info["faces"], info["volume"]), (verdict[0], verdict[1] % (solid - info["volume"])))
    return figure


def figure_rounded(rebuild):
    kentron, cadquery = (scene("rounded", k, rebuild) for k in ("kentron", "cadquery"))
    height = 250
    figure = Figure(WIDTH, TOP + height + 40, "Removing every fillet of a box rounded on all its edges: kentron leaves no faces, where the box is expected.")
    view = View((0.8, -1, 0.62), (0, 0, 0), 21, SIDE, height)
    larger = view.zoomed((0, 0, 0), 23, MAIN, height)
    box = kentron["box"]
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([kentron["rounded"]], view, hide_smooth=True), "Rounded box", "its 12 fillets and 8 corners are removed")
    silhouette = "".join(traced(larger, e["p"], "ghost") for e in box["edges"])
    figure.panel(COLUMNS[1], TOP, MAIN, height, silhouette, "kentron", "defeature::remove_faces", ("bad", "no faces are left"))
    figure.text(COLUMNS[1] + MAIN / 2, TOP + height - 16, "nothing: the box's place, dashed", "note small halo", "middle")
    figure.panel(COLUMNS[2], TOP, MAIN, height, render([box], larger), "Expected", "as SMLib, OCCT and sampling give", ("good", "the box: 6 faces, volume %.0f" % box["info"][0]["volume"]))
    return figure


def figure_twisted(rebuild):
    kentron = scene("twisted", "kentron", rebuild)
    height = 300
    figure = Figure(WIDTH, TOP + height + 56, "Two twisted prisms of a knurl pass through each other where their sections coincide: kentron's intersection of the pair is an open shell of four faces, not the small lens they share.")
    whole = View((0.866, -0.5, 0.42), (0, 0, 5), 62, SIDE, height)
    close = whole.zoomed((23.17, -13.38, 7.5), 3.6, MAIN, height)
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([kentron["sets"]], whole), "A knurl's prisms", "12 twisted by +40°, 12 by −40°")
    figure.frame(*(a + b for a, b in zip(whole.window(close), (COLUMNS[0], TOP, 0, 0))))
    figure.panel(COLUMNS[1], TOP, MAIN, height, render([kentron["a"], kentron["b"]], close), "A pair, where they cross", "the same section at z = 7.5")
    figure.text(COLUMNS[1] + 2, TOP + height + 20, "they share a small lens: volume 0.13 by sampling,\nwhere a prism's is 4.33", "note small")
    common, union = kentron["common"], kentron["union"]
    ghost = "".join(traced(close, e["p"], "ghost") for part in ("a", "b") for e in kentron[part]["edges"])
    figure.panel(COLUMNS[2], TOP, MAIN, height, ghost + render([common], close) + marked(close, common), "kentron: their intersection", signed("%d faces, volume %.2f" % (common["info"][0]["faces"], common["info"][0]["volume"])),
                 ("bad", "an open shell: the edges in red bound one face;\ntheir union: %d faces, not valid, volume %.2f" % (union["info"][0]["faces"], union["info"][0]["volume"])))
    return figure


def figure_lamp(rebuild):
    kentron, cadquery = (scene("lamp", k, rebuild) for k in ("kentron", "cadquery"))
    height = 300
    figure = Figure(WIDTH, TOP + height + 56, "A spiral arm swept through a hub cut by a second sweep along the same spiral: kentron's union of the two is not valid, CadQuery's is one sound solid.")
    whole = View((0.3, 0.8, 0.7), (0, 0, 0), 400, SIDE, height)
    close = whole.zoomed((36, 4, 68), 62, MAIN, height)
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([cadquery["arm"], cadquery["hub"]], whole), "Spiral arm and its hub", "two sweeps along one spiral")
    figure.frame(*(a + b for a, b in zip(whole.window(close), (COLUMNS[0], TOP, 0, 0))))
    for x, mesh, title, verdict in (
            (COLUMNS[1], kentron["union"], "kentron: their union", ("bad", "not valid: the edges in red, where the arm lies\nin the hub's channel, bound 1, 3, 4 or 6 faces")),
            (COLUMNS[2], cadquery["union"], "CadQuery: their union", ("good", "one valid solid"))):
        info = mesh["info"][0]
        figure.panel(x, TOP, MAIN, height, render([mesh], close) + (marked(close, mesh) if verdict[0] == "bad" else ""), title, "%d faces, volume %.0f" % (info["faces"], info["volume"]), verdict)
    return figure


def figure_mirrored(rebuild):
    import numpy as np

    kentron = scene("mirrored", "kentron", rebuild)
    union, plane = kentron["union"], 1.99
    height = 300
    figure = Figure(WIDTH, TOP + height + 56, "A free-form plate united with its mirror image: the plate's side face crosses the mirror plane at a small angle, and kentron's union keeps the parts of the two side faces that lie inside the other half.")
    whole = View((1, -0.7, 0.6), (plane, 1, 1), 5.2, SIDE, height)
    close = View((1, -0.62, 0.42), (plane, 0.12, 1.2), 1.05, MAIN, height)
    drawing = render([kentron["half"], kentron["image"]], whole)
    figure.panel(COLUMNS[0], TOP, SIDE, height, drawing, "A plate and its mirror image", "mirrored in the plane x = 1.99")
    figure.callout(at(whole, COLUMNS[0], TOP, (plane, 0, 1.2)), "where they meet:\nthe next panels", (COLUMNS[0] + 10, TOP + height - 30))
    # The fragments: small faces of the union lying by the mirror plane, inside the solid.
    fragments = [f for f in union["faces"] if 0.03 < f["area"] < 0.04 and np.abs(f["v"][:, 0] - plane).max() < 0.004]
    inside = {"faces": [dict(f, material=WRONG) for f in fragments], "edges": []}
    info = union["info"][0]
    figure.panel(COLUMNS[1], TOP, MAIN, height, render([union], close) + render([inside], close, outline=False), "kentron: their union", "%d faces, area %.3f" % (info["faces"], info["area"]),
                 ("bad", "%d fragments of the side faces are kept\ninside the solid: 0.14 of area too much" % len(fragments)))
    figure.callout(at(close, COLUMNS[1], TOP, (plane, 0.04, 1.3)), "a fragment, drawn through\nthe solid it lies inside", (COLUMNS[1] + 205, TOP + 254), "middle")

    # A section across the mirror plane, the distances from it drawn larger.
    level, reach, stretch = 1.25, 0.2, 30
    chart = Chart(MAIN, height, (0, reach / 2 - 0.022), (height - 80) / reach, stretch)
    # The plate's side face by the mirror plane: within 0.03 of it all over.
    side = sliced([f for f in kentron["half"]["faces"] if np.abs(f["v"][:, 0] - plane).max() < 0.03 and np.ptp(f["v"][:, 1]) > 0.3], 2, level).reshape(-1, 3)
    side = side[side[:, 1] < reach]
    side = side[np.argsort(side[:, 1])]
    curve = [(side[0, 0] - plane, 0.0)] + [(x - plane, y) for x, y, _ in side] + [(side[-1, 0] - plane, reach)]
    crossing = next(y for (a, _), (b, y) in zip(curve, curve[1:]) if a > 0 >= b)
    edge = 0.0042
    mirror = [(-x, y) for x, y in curve]
    kept = [p for p in curve if p[1] < crossing] + [(0.0, crossing)]
    drawing = (chart.path([(-edge, 0)] + curve + [(-edge, reach)], "wash first", True) + chart.path([(edge, 0)] + mirror + [(edge, reach)], "wash second", True)
               + chart.path([(0, -0.012), (0, reach)], "dash") + chart.path([(-edge, 0)] + curve, "one") + chart.path([(edge, 0)] + mirror, "two")
               + chart.path(kept, "mark") + chart.path([(-x, y) for x, y in kept], "mark"))
    x0, y0 = COLUMNS[2], TOP
    figure.panel(x0, y0, MAIN, height, drawing, "Why: a section across the plane", "distances from the plane drawn %d× larger" % stretch)
    figure.text(*chart.at(-edge * 0.56, reach * 0.8, (x0, y0)), "the plate", "small", "middle")
    figure.text(*chart.at(edge * 0.56, reach * 0.8, (x0, y0)), "its mirror image", "small", "middle")
    figure.text(*chart.at(0, -0.027, (x0, y0)), "mirror plane", "note small", "middle")
    figure.callout(chart.at(0, crossing, (x0, y0)), "the side\nfaces cross\nat %.2f rad" % (2 * math.atan2(curve[0][0], crossing)), (x0 + 236, y0 + 180), "middle")
    figure.text(x0 + 2, y0 + height + 20, "in red: the parts of the side faces that kentron\nkeeps, each inside the other half", "note small")
    return figure


def figure_gyroid(rebuild):
    kentron = scene("gyroid", "kentron", rebuild)
    height = 300
    corner = (10, 10, 10)
    figure = Figure(WIDTH, TOP + height + 56, "Four thin plates of a gyroid cell meet at a corner of the cell, their side strips a hundred-millionth apart: kentron's union of them is not valid, with an edge 3.5e-5 long used by seven faces at the corner.")
    whole = View((1, -0.7, 0.6), (8.5, 9.5, 9.6), 31, SIDE, height)
    close = whole.zoomed((9.75, 10, 10.1), 4.2, MAIN, height)
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([kentron["plates"]], whole), "Four plates of a gyroid cell", "each 0.1 thick")
    figure.frame(*(a + b for a, b in zip(whole.window(close), (COLUMNS[0], TOP, 0, 0))))
    figure.panel(COLUMNS[1], TOP, MAIN, height, render([kentron["plates"]], close), "Where they meet", "a corner of the cell, (10, 10, 10)")
    figure.callout(at(close, COLUMNS[1], TOP, (9.2, 10, 10.06)), "side strips of two plates,\n1e-8 apart", (COLUMNS[1] + 12, TOP + height - 30))
    info = kentron["union"]["info"][0]
    figure.panel(COLUMNS[2], TOP, MAIN, height, render([kentron["union"]], close), "kentron: their union", "%d faces" % info["faces"],
                 ("bad", "not valid: at the corner, an edge 3.5e-5 long\nis used by seven faces"))
    x, y = at(close, COLUMNS[2], TOP, corner)
    figure.ring((x, y), 9)
    # The edge seen end-on, as a sketch: a sound edge has two faces.
    cx, cy, r = COLUMNS[2] + MAIN - 62, TOP + height - 62, 46
    figure.add('<path class="leader" d="%s"/><circle class="panel" cx="%s" cy="%s" r="%d"/><circle class="thin" cx="%s" cy="%s" r="%d"/>\n' % (path_data([(x + 7, y + 6), (cx - 30, cy - 34)]), number(cx), number(cy), r, number(cx), number(cy), r))
    for angle in (12, 58, 97, 163, 201, 262, 318):
        a = math.radians(angle)
        figure.add('<path class="ink" d="%s"/>\n' % path_data([(cx, cy), (cx + (r - 9) * math.cos(a), cy - (r - 9) * math.sin(a))]))
    figure.add('<circle class="bad" cx="%s" cy="%s" r="3.6"/>\n' % (number(cx), number(cy)))
    figure.text(cx - r - 8, cy + 4, "the edge, end-on\n(a sketch): seven\nfaces, not two", "note small halo", "end")
    return figure


def figure_evolved(rebuild):
    kentron, occt = (scene("evolved", k, rebuild) for k in ("kentron", "cadquery"))
    height = 300
    figure = Figure(WIDTH, TOP + height + 56, "A profile leaning at 45 degrees to its spine, swept along it (OCCT's test AGT003): at each corner of the spine OCCT's solid has two faces filling in between the sweeps of the two sides, which kentron's lacks.")
    eye = (1, -0.45, 0.6)
    whole = View(eye, (10, 5, 4), 29, SIDE, height)
    close = whole.zoomed((0.5, -0.2, 6.6), 8.6, MAIN, height)
    spine, profile = kentron["spine"]["points"], kentron["profile"]["points"]
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([occt["result"]], whole) + traced(whole, spine, "sketch") + traced(whole, profile, "two"), "A profile swept along a spine", "OCCT's test evolved/voluved/AGT003")
    figure.frame(*(a + b for a, b in zip(whole.window(close), (COLUMNS[0], TOP, 0, 0))))
    figure.callout(at(whole, COLUMNS[0], TOP, (16.5, 0, 0)), "the spine", (COLUMNS[0] + SIDE - 10, TOP + height - 14), "end")
    figure.callout(at(whole, COLUMNS[0], TOP, (0, 1, 8.5)), "the profile, in a plane\nat 45° to the spine", (COLUMNS[0] + SIDE / 2 + 16, TOP + 26), "middle")
    fill = {i: FIRST for i, f in enumerate(occt["result"]["faces"]) if 3.0 < f["area"] < 3.6}
    for x, mesh, title, recolour, reference, verdict in (
            (COLUMNS[1], kentron, "kentron", None, [occt["result"]], ("bad", "area %.2f: no faces fill in at the corners;\nin red, where OCCT's join the walls")),
            (COLUMNS[2], occt, "OCCT", fill, None, ("good", "area %.2f: two faces at each corner fill in\nbetween what is left of the sides' sweeps"))):
        info = mesh["result"]["info"][0]
        figure.panel(x, TOP, MAIN, height, render([mesh["result"]], close, recolour, reference=reference), title, "%d faces, volume %.3f" % (info["faces"], info["volume"]), (verdict[0], verdict[1] % info["area"]))
    figure.callout(at(close, COLUMNS[2], TOP, (0.5, -0.25, 4.6)), "a plane of 3.125\nand a trimmed\nsurface of 3.509", (COLUMNS[2] + 66, TOP + 96), "middle")
    return figure


def figure_throat(rebuild):
    import numpy as np

    kentron, occt = (scene("throat", k, rebuild) for k in ("kentron", "cadquery"))
    height = 300
    figure = Figure(WIDTH, TOP + height + 56, "A chamfer of constant throat round a cylinder standing on a tilted box (OCCT's test bug28828_2): OCCT's chamfer is wider where the faces meet at a wider angle; kentron's has the distances of the edge's middle all round.")
    whole = View((1, -0.75, 0.5), (0, 0, 12), 92, SIDE, height)
    close = View((1, -0.2, 0.3), (0, -5, 19.5), 36, MAIN, height)
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([kentron["fused"]], whole), "A cylinder on a tilted box", "OCCT's test bug28828_2")
    figure.callout(at(whole, COLUMNS[0], TOP, (6, -13, 18.5)), "the edge chamfered,\nwith a throat of 2", (COLUMNS[0] + 10, TOP + height - 30))
    info = kentron["result"]["info"][0]
    plain = {i: BODY for i in range(info["faces"])}
    figure.panel(COLUMNS[1], TOP, MAIN, height, render([kentron["result"]], close, plain, reference=[occt["result"]]), "kentron", "%d faces, volume %.1f" % (info["faces"], info["volume"]),
                 ("bad", "volume %.1f for OCCT's %.1f; in red,\nwhere the chamfer is off OCCT's" % (info["volume"], occt["result"]["info"][0]["volume"])))

    # The sections in the plane x = 0, through the cylinder's axis: the edge's lowest and
    # highest points, side by side.
    tilt = math.tan(math.radians(20))
    half, scale = MAIN / 2, 23.0
    x0, y0 = COLUMNS[2], TOP
    drawing = '<path class="rule" d="M%s 10V%s"/>\n' % (number(half), number(height - 10))
    labels = []
    for column, wall, away, angle in ((0, -15.0, -1, 110), (1, 5.0, 1, 70)):
        ends = {}
        for name, mesh in (("kentron", kentron), ("OCCT", occt)):
            cut = sliced([f for f in mesh["result"]["faces"] if f["material"] == FIRST], 0, 0.0).reshape(-1, 3)[:, 1:]
            cut = cut[((cut[:, 0] - wall) * away > -1e-3) & (np.abs(cut[:, 0] - wall) < 6)]
            ends[name] = cut[np.abs(cut[:, 0] - wall).argmin()], cut[np.abs(cut[:, 0] - wall).argmax()]       # On the cylinder, on the box.
        on_cylinder, on_box = ends["OCCT"]
        corner = np.array([wall, on_box[1] + (wall - on_box[0]) * tilt])
        chart = Chart(half, height, (wall + away * 1.3, corner[1] + 1.0), scale)
        reach = 3.7
        far, top, inside = corner + away * reach * np.array([1, tilt]), corner + np.array([0, reach + 0.5]), wall - away * 1.5
        floor = corner[1] - reach * tilt - 1.5
        part = (chart.path([far, on_box, on_cylinder, top, (inside, top[1]), (inside, floor), (far[0], floor)], "cut", True) + chart.path([far, on_box], "ink") + chart.path([on_cylinder, top], "ink")
                + chart.path([on_box, corner, on_cylinder], "dash"))
        origin = (x0 + column * half, y0)
        outer = max(ends, key=lambda name: abs(ends[name][1][0] - wall))
        for name, cls in (("OCCT", "one"), ("kentron", "mark")):
            a, b = ends[name]
            along = (b - a) / np.linalg.norm(b - a)
            foot = a + along * ((corner - a) @ along)
            throat = float(np.linalg.norm(corner - foot))
            if name == outer:
                part += chart.path([corner, foot], "thin")
            part += chart.path([a, b], cls)
            text = "%s:\nthroat %s" % (name, ("%.2f" % throat).rstrip("0").rstrip("."))
            centre = origin[0] + half / 2 + away * 28
            labels.append((chart.at(*(a + (b - a) * (0.5 if name == outer else 0.8)), origin), text, (centre, y0 + 64) if name == outer else (centre, y0 + 262)))
        drawing += '<g transform="translate(%s 0)">%s</g>\n' % (number(column * half), part)
        labels.append((origin[0] + half / 2, y0 + 24, "faces at %d°" % angle))
    figure.panel(x0, y0, MAIN, height, drawing, "Why: sections across the edge", "at its lowest point and at its highest")
    for label in labels:
        if len(label) == 3 and isinstance(label[2], str):
            figure.text(label[0], label[1], label[2], "small", "middle")
        else:
            figure.callout(label[0], label[1], label[2], "middle")
    figure.text(x0 + 2, y0 + height + 20, "throat: the chamfer's distance from the edge;\nkentron's meets both faces 2.83 from it", "note small")
    return figure


def figure_revolved(rebuild):
    kentron, expected = (scene("revolved", k, rebuild) for k in ("kentron", "cadquery"))
    height = 270
    figure = Figure(WIDTH, TOP + height + 56, "A profile revolved by 270 degrees and hollowed by 1 inwards, open at one of its side faces (OCCT's test offset/faces_type_i/F3): kentron's offset fails, where a hollow solid with walls 1 thick is expected.")
    view = View((1, -0.62, 0.55), (0, 0, -1), 58, SIDE, height)
    larger = view.zoomed((0, 0, -1), 60, MAIN, height)
    source = kentron["source"]
    figure.panel(COLUMNS[0], TOP, SIDE, height, render([source], view), "A profile revolved by 270°", "OCCT's test offset/faces_type_i/F3")
    figure.callout(at(view, COLUMNS[0], TOP, (16, 0, -6)), "hollowed by 1 inwards,\nopen at this side", (COLUMNS[0] + SIDE - 10, TOP + height - 30), "end")
    ghost = "".join(traced(larger, e["p"], "ghost") for e in source["edges"] if e["kind"] == "sharp")
    figure.panel(COLUMNS[1], TOP, MAIN, height, ghost, "kentron", "offset::thick_solid", ("bad", "fails: “an offset leaving nothing”"))
    figure.text(COLUMNS[1] + MAIN / 2, TOP + height - 16, "nothing: the solid's place, dashed", "note small halo", "middle")
    info = expected["hollow"]["info"][0]
    figure.panel(COLUMNS[2], TOP, MAIN, height, render([expected["hollow"]], larger), "Expected", "the profile's exact offset, revolved and cut",
                 ("good", "walls 1 thick: %d faces, volume %.2f\n(OCCT's own offset builds no wall)" % (info["faces"], info["volume"])))
    return figure


def figure_wire(rebuild):
    import numpy as np
    import shapely
    from shapely.geometry import LineString, MultiLineString, Point
    from shapely.ops import nearest_points

    data = scene("wire", "cadquery", rebuild)
    wire, contour, boundary, distance = data["wire"], data["contour"], data["boundary"], data["distance"]
    tip = wire[np.argmax(wire[:, 0])]
    top, lower = 150, 220
    y1 = TOP + top + 62
    figure = Figure(WIDTH, y1 + lower + 56, "An open wire turning back at its tip far more tightly than the offset distance (OCCT's test wire offset_wire_058), offset by 2.5 per cent of its length: by the tip, kentron's contour cuts inside the true one by up to 0.68.")
    full = WIDTH - 2 * MARGIN
    chart = Chart(full, top, (-182, -6), full / 1290)
    figure.panel(MARGIN, TOP, full, top, chart.path(wire, "ink") + "".join(chart.path(piece, "one") for piece in contour), "An open wire and its offset",
                 "OCCT's test offset/wire_unclosed_outside_0_025/B6: the wire offset_wire_058, offset by 2.5%% of its length, %.1f" % distance)
    origin = (MARGIN, TOP)
    lowest = np.concatenate(contour)
    lowest = min((p for p in lowest if abs(p[0] + 80) < 2), key=lambda p: p[1])
    figure.callout(chart.at(*wire[np.abs(wire[:, 0] + 420).argmin()], origin), "the wire", (MARGIN + 200, TOP + 20), "end")
    figure.callout(chart.at(*lowest, origin), "its offset: the contour %.1f from it" % distance, (MARGIN + 500, TOP + top - 10))
    figure.callout(chart.at(*tip, origin), "the tip, where the wire turns back", (MARGIN + full - 110, TOP + 20), "end")
    figure.scale(MARGIN + 14, TOP + top - 10, 100 * chart.scale, "100")

    # The tip, and the contour by it.
    left, right = 300, full - 300 - GAP
    zoom = 10000
    close = Chart(left, lower, (tip[0] - 0.0105, tip[1] - 0.001), chart.scale * zoom)
    figure.panel(MARGIN, y1, left, lower, close.path(data["tip"], "ink"), "The tip", "drawn %s× larger" % format(zoom, ",").replace(",", " "))
    figure.scale(MARGIN + left - 14 - 0.01 * close.scale, y1 + lower - 12, 0.01 * close.scale, "0.01")
    figure.text(MARGIN + 2, y1 + lower + 20, "a corner and a turn 0.02 across,\nwhere the offset distance is %.1f" % distance, "note small")

    x2 = MARGIN + left + GAP
    centre = (400.6, -35.6)
    stray = Chart(right, lower, centre, right / 25)
    # The piece of kentron's contour farthest from the true one, and where.
    true = MultiLineString([LineString(loop) for loop in boundary])
    strays = [shapely.distance(true, shapely.points(piece)) for piece in contour]
    most = max(range(len(contour)), key=lambda i: strays[i].max())
    piece, off = contour[most], strays[most].max()
    worst = piece[strays[most].argmax()]
    nearest = np.array(nearest_points(true, Point(worst))[0].coords[0])
    others = max(s.max() for i, s in enumerate(strays) if i != most)
    drawing = "".join(stray.path(loop, "one") for loop in boundary) + "".join(stray.path(p, "mark") for p in contour) + stray.path([worst, nearest], "thin")
    figure.panel(x2, y1, right, lower, drawing, "kentron's contour by the tip", "drawn %.0f× larger" % (stray.scale / chart.scale),
                 ("bad", "length %.2f for %.3f: a straight piece cuts inside the round of the tip\nby up to %.2f; elsewhere the contour is within %.2f of the true one"
                  % (data["contour_length"], sum(np.linalg.norm(np.diff(loop, axis=0), axis=1).sum() for loop in boundary), off, others)))
    window = (*chart.at(centre[0] - 12.5, centre[1] + lower / stray.scale / 2, origin), 25 * chart.scale, lower / stray.scale * chart.scale)
    figure.frame(*window)
    figure.callout((window[0] + window[2] / 2, window[1] + window[3]), "by the tip: below", (MARGIN + full - 14, TOP + top - 10), "end")
    figure.scale(x2 + 14, y1 + lower - 12, stray.scale, "1")
    along = piece[-1] + 0.62 * (piece[np.linalg.norm(piece - piece[-1], axis=1).argmax()] - piece[-1])
    figure.callout(stray.at(*along, (x2, y1)), "kentron's contour:\na straight piece %.1f long" % np.linalg.norm(piece - piece[-1], axis=1).max(), (x2 + 104, y1 + 30), "middle")
    beyond = nearest_points(true, Point(centre[0] + 5.5, centre[1] - 1.5))[0].coords[0]
    figure.callout(stray.at(*beyond, (x2, y1)), "the true contour: the points\n%.1f from the wire" % distance, (x2 + right - 14, y1 + lower - 50), "end")
    mx, my = stray.at(*((worst + nearest) / 2), (x2, y1))
    figure.text(mx + 9, my + 12, "%.2f" % off, "note small halo")
    return figure


def figure_tolerance(rebuild):
    """A sketch, not a model: a fitted curve between its fit points."""
    import numpy as np

    width, height = WIDTH - 2 * MARGIN, 190
    figure = Figure(WIDTH, TOP + height + 40, "A sketch of a fitted section curve: it passes through its fit points on the surfaces' intersection, and between them strays from it by more than the tolerance recorded for the edge.")
    x0, y0 = MARGIN, TOP
    left, right, middle = 40, width - 40, 100
    strays, stray, tolerance = (1.0, 0.6, 0.8, 0.45), 5.75e-7, 1e-7
    unit = 62 / stray
    u = np.linspace(0, len(strays), 801)
    x = left + (right - left) * u / len(strays)
    y = middle - unit * stray * np.sin(math.pi * u) * np.interp(u, np.arange(len(strays)) + 0.5, strays)
    curve = np.stack([x, y], axis=1)
    normal = np.gradient(curve, axis=0)[:, ::-1] * (1, -1)
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    band = np.concatenate([curve + unit * tolerance * normal, (curve - unit * tolerance * normal)[::-1]])
    from scipy.spatial import cKDTree

    up = normal * -np.sign(normal[:, 1:])
    line = np.stack([x, np.full_like(x, middle)], axis=1)
    outside = cKDTree(curve).query(line)[0] > unit * tolerance
    drawing = '<path class="wash first" d="%s"/>\n' % path_data(band, True) + '<path class="ink" d="%s"/>\n' % path_data(line[[0, -1]])
    edges = np.flatnonzero(np.diff(np.concatenate([[False], outside, [False]]).astype(int)))
    for a, b in zip(edges[::2], edges[1::2]):
        drawing += '<path class="mark" d="%s"/>\n' % path_data(line[[a, b - 1]])
    drawing += '<path class="one" d="%s"/>\n' % path_data(simplify(curve, 0.1))
    peak = curve[np.argmin(y)]
    drawing += '<path class="thin" d="%s"/>\n' % path_data([(peak[0], middle), peak])
    for i in range(len(strays) + 1):
        drawing += '<circle class="dot" cx="%s" cy="%s" r="3.4"/>\n' % (number(left + (right - left) * i / len(strays)), number(middle))
    figure.panel(x0, y0, width, height, drawing, "A fitted section curve and the tolerance recorded for its edge", "a sketch with the teapot spout's numbers: the distances across drawn far larger than those along",
                 ("bad", "in red: where the surfaces are farther from the curve than its tolerance says; a later boolean may not find the edge on them"))
    span = (right - left) / len(strays)
    figure.text(x0 + peak[0] + 7, y0 + middle - 9, format(stray, ".2e").replace("e-0", "e-"), "note small")
    figure.callout((x0 + left + 0.8 * span, y0 + curve[160, 1]), "the fitted curve\nof the section edge", (x0 + left + 1.5 * span, y0 + 22), "middle")
    top = curve[560] + unit * tolerance * up[560]
    figure.callout((x0 + top[0], y0 + top[1]), "within the tolerance recorded:\n%s from the curve" % format(tolerance, ".0e").replace("e-0", "e-"), (x0 + left + 3.5 * span, y0 + 22), "middle")
    figure.callout((x0 + left + 0.38 * span, y0 + middle), "the surfaces' intersection", (x0 + left + 0.5 * span, y0 + middle + 56), "middle")
    figure.callout((x0 + left + 2 * span, y0 + middle + 2), "fit points,\non the intersection", (x0 + left + 2 * span + 30, y0 + middle + 56), "middle")
    return figure


def figure_copies(rebuild):
    """A diagram: the faces a solid with voids holds, where the voids are one box placed twice."""
    width, height = (WIDTH - 2 * MARGIN - GAP) / 2, 200
    figure = Figure(WIDTH, TOP + height + 40, "A solid with two voids made from one small box placed twice: OCCT's solid holds the placed faces themselves, kentron's holds copies of them, which the operation's history does not list.")
    for column, title, note, verdict in (
            (0, "OCCT: the solid holds the placed faces", "solid(faces, inner_faces, history=…)", ("good", "every face of the solid is among the history's images")),
            (1, "kentron: the solid holds copies", "make_solid copies a placed shell", ("bad", "the solid's 12 inner faces are not among the history's images"))):
        x, y = MARGIN + column * (width + GAP), TOP
        figure.panel(x, y, width, height, "", title, note, verdict)
        figure.node(x + 14, y + 74, 92, 52, "a small box:\n6 faces")
        rows = (y + 48, y + 118)
        if column == 0:
            figure.node(x + 176, y + 14, width - 190, height - 28, "", "group")
            figure.text(x + 176 + (width - 190) / 2, y + 34, "the solid", "small", "middle")
            for row, place in zip(rows, "AB"):
                figure.node(x + 190, row, width - 218, 50, "its faces, placed at %s:\ngiven, and in the history" % place)
                figure.arrow((x + 106, y + 100), (x + 190, row + 25))
        else:
            figure.node(x + 284, y + 14, width - 298, height - 28, "", "group")
            figure.text(x + 284 + (width - 298) / 2, y + 34, "the solid", "small", "middle")
            for row, place in zip(rows, "AB"):
                figure.node(x + 140, row, 112, 50, "placed at %s:\nin the history" % place)
                figure.arrow((x + 106, y + 100), (x + 140, row + 25))
                figure.node(x + 294, row, width - 318, 50, "a copy:\n6 new faces", "node off")
                figure.arrow((x + 252, row + 25), (x + 294, row + 25), "dash")
    return figure


FIGURES = {name[7:]: f for name, f in list(globals().items()) if name.startswith("figure_")}


def main(args):
    if args[:1] == ["--dump"]:
        return dump_scene(args[1], args[2])
    if "--list" in args:
        return print("\n".join(FIGURES))
    rebuild = "--rebuild" in args
    names = [a for a in args if not a.startswith("--")] or list(FIGURES)
    OUT.mkdir(parents=True, exist_ok=True)
    for name in names:
        path = OUT / (name.replace("_", "-") + ".svg")
        path.write_text(FIGURES[name](rebuild).svg())
        print("%s: %.1f kB" % (path.relative_to(REPO), path.stat().st_size / 1000))


if __name__ == "__main__":
    main(sys.argv[1:])

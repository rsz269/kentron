# Test cases in YAML

Tests that build shapes and check measurable facts about them are written as data, independent
of the modeling kernel that runs them: no API calls, and no names that depend on how a kernel
enumerates sub-shapes. `crates/occt-tests/cases/**/*.yaml` holds them; the `cases` test target
runs each against kentron.

A file holds one or more cases, as YAML documents separated by `---`. Cases ported from OCCT live
in `cases/occt`: `draw/` is generated from its DRAW scripts by `tools/draw2rs.py` (do not edit;
deviations go to `draw_deviations.json`, sub-shape positions to `draw_selectors.json`), and
`gtest/` was moved from the hand-translated GTests by `tools/gtest2yaml.py` (edit freely).

## A case

```yaml
case: blend_simple_f4                 # unique name (the test's name)
ported_from:                          # for tests ported from another project
  project: OCCT
  repo: https://github.com/Open-Cascade-SAS/OCCT
  git_hash: 3d097a0328e71b826377d4814ab05ec3c3d23871
  file: tests/blend/simple/F4       # for GTests also `test: Fixture.Name`
tags: [blend]                         # optional
params: { r: 10 }                     # optional named constants for expressions
skip: "reason"                        # optional: not run (missing data, known failure upstream)
expect: ok                            # ok (default) | not_done (a known gap; passing is reported)
timeout: 5                            # optional: its own limit on processor time, in seconds
build:                                # steps, in order
  - s: { box: { size: [100, 100, 10] } }
  - check: { valid: s }
  - result: { fillet: { shape: s, radius: r, edges: [e1] } }
checks:                               # all must hold after the build
  - area: result
    eq: 23914.2
    rel: 0.01
```

A build step is a one-key map: `name: operation` binds the operation's value to `name` (a name may
be bound again; later steps see the new value; an operation that is a bare name binds the same
value under another name), and `check: assertion` asserts there and then. `notes` holds free
text (deviations from the original test, upstream known failures).

## Values

Cases generated from OCCT's scripts inline the scripts' variables (such as the blend group's
`SCALE1`): arithmetic on numbers is written as its value, and only expressions with `pi` or
functions stay expressions.

- Numbers: YAML numbers, or expression strings over `params`: `"1*SCALE1"`, `"pi/2"`,
  `"2*r + 1"`, with `+ - * /`, `^` or `**`, parentheses, `pi`, `sqrt sin cos tan asin acos atan
  abs exp log`, and `atan2 pow min max` of two arguments.
- Points and vectors: `[x, y, z]`; 2D: `[x, y]`. Directions: vectors (normalized) or `x`, `-y`,
  ... Angles are in degrees; curve parameters (edge ranges) are plain numbers.
- A frame: `{origin: P, z: D, x: D}`, `x` optional (chosen as OCCT's `gp_Ax2` does), or the name of
  a bound plane.
- A name refers to a bound value; where a list is expected, a single item may stand alone.

## Operations

Primitives (`frame` optional, default the world frame):

| operation | arguments |
|---|---|
| `box` | `size: [dx, dy, dz]`, `corner: P` (default origin) |
| `wedge` | `size: [dx, dy, dz]`, and `top_dx` or `top: [xmin, zmin, xmax, zmax]` |
| `cylinder` | `radius`, `height`, `angle` |
| `cone` | `r1`, `r2`, `height`, `angle` |
| `sphere` | `radius`, `latitudes: [a1, a2]`, `angle` |
| `torus` | `major`, `minor`, `tube: [a1, a2]`, `angle` |
| `plane` | `frame` (a plane, also usable as a frame) |
| `half_space` | `face`, `toward: P` |

Curves (3D, or 2D given 2D points):

| operation | arguments |
|---|---|
| `circle` | `frame` (or `center`, `normal`), `radius` |
| `ellipse` | `frame` (or `center`, `normal`), `major`, `minor` |
| `parabola` | `frame` (or `center`, `normal`), `focal` |
| `hyperbola` | `frame` (or `center`, `normal`), `major`, `minor` |
| `line` | `origin`, `direction` |
| `bezier` | `poles`, `weights` |
| `bspline` | `degree`, `knots: [[k, multiplicity], ...]`, `poles`, `weights`, `periodic` |

Topology:

| operation | arguments |
|---|---|
| `vertex` | `at: P` |
| `edge` | `curve` (with `range: [t1, t2]`), or `from`, `to` (points or vertices) |
| `polyline` | `points` (closed if the last is the first) |
| `wire` | `edges` (edges or wires, in any order) |
| `face` | `wire` (planar) |
| `compound` | `of: [...]` |
| `empty` | `solid`, `shell`, `compound`, ... |
| `add` | `to`, `shape` (adds to a container) |
| `load` | `brep: file`: a shape in OCCT's BRep format, or `step: file`: one in a STEP file, from OCCT's test data or a BRL-CAD source checkout (external downloads: see `crates/occt-tests/README.md`; a case whose file is missing is skipped) |
| `select` | a sub-shape by position (below) |

Transforms (each makes a new value): `translate` (`shape`, `by`), `rotate` (`shape`, `axis:
{origin, direction}`, `angle`), `scale` (`shape`, `center`, `factor`), `reverse` (`shape`),
`orient` (`shape`, `as: forward | reversed | internal | external`), `copy` (`shape`).

Sweeps: `extrude` (`profile`, `by`), `revolve` (`profile`, `axis`, `angle`), `pipe` (`profile`,
`spine`). Conversions: `to_nurbs` (`shape`).

Sketches: `profile` (`frame`, `steps`, `outcome`) draws with a pen in the frame's XY plane. The
pen starts at the origin heading along X; each step is a one-key map:

| step | does |
|---|---|
| `start: [x, y]` | puts the pen there (before anything is drawn); the profile closes back to it |
| `origin: P` | moves the plane's origin, keeping its axes |
| `plane: {normal: V, x: V}` | turns the plane, keeping its origin |
| `on: f` | draws on the plane of a planar face (or a bound plane) |
| `along_x: l`, `along_y: l` | a line along X (Y) by `l`; the pen then heads that way |
| `to_x: x`, `to_y: y` | a line along X (Y) to a coordinate |
| `forward: l` | a line ahead (negative: back, then heading back) |
| `forward_to_x: x`, `forward_to_y: y` | a line ahead until the pen reaches a coordinate |
| `by: [dx, dy]`, `to: [x, y]` | a line by a vector, or to a point; the pen heads along it |
| `turn: a`, `heading: a`, `direction: [dx, dy]` | turns the pen by `a` degrees (left positive), to `a` degrees from X, or to a direction |
| `arc: {radius: r, angle: a}` | an arc tangent to the heading, its centre on the left (right for a negative radius), through `a` degrees (the other way for a negative angle) |

`outcome` is `face` (default: closed back to the start, a face on the plane), `closed_wire` or
`open_wire` (the edges as drawn).

Booleans: `union`, `difference` (`object` minus `tool`), `intersection` (`object`, `tool`, each
a shape or a list of shapes), `section` (the same arguments: the edges where the operands meet,
and the points where they touch), `partition` (`of: [...]`: every argument split by the others),
`split` (`objects`, `tools`).

Lofts: `loft` (`sections: [...]` wires or vertices, `solid`, `ruled`, `matching`: how the
sections' edges are matched up, `by_rays` as OCCT's thrusections by default, `as_given` for its
`-N`, or `by_length`).

Wire offsets: `offset_wire` (`shape`, and `distance` or `length_fraction`, the distance as a fraction of
the wire's length).

Shape offsets: `offset_shape` (`shape`, `distance`, `join`: `arc` or `intersection`) and `thick_solid`
(`shape`, `faces` removed, `thickness`, `join`).

Sewing: `sew` (`shapes`, `tolerance` (default 1e-6), `min_tolerance`, `max_tolerance`, and the
switches of OCCT's `BRepBuilderAPI_Sewing`, each true or false: `sewing`, `analysis`, `cutting`,
`non_manifold`, `same_parameter`, `floating_edges`, `faces`). kentron's `sew::sew` takes `tolerance`,
`cutting` and `non_manifold` (not implemented); the other arguments are accepted and ignored.

Fillets:

| operation | arguments |
|---|---|
| `fillet` | `shape`, `edges: [{edge, radius}]` or `edges` with one `radius`, or `edges: [{edge, law: [[u, r], ...]}]` for a radius varying along the edge |
| `filleted_union`, `filleted_difference` | `object`, `tool`, `radius`: the union or difference, the edges where the operands met filleted |

## Selecting sub-shapes

```yaml
e1: { select: { of: s, edge: [160, 0, 50] } }               # the edge whose centroid is there
e2: { select: { of: s, edge: [0, 0, 0], length: 125.66 } }  # several share it: by length
f1: { select: { of: s, face: [5, 5, 10], area: 100 } }      # faces: by area
v1: { select: { of: s, vertex: [0, 0, 10] } }
s1: { select: { of: a, solid: [0, 0, 5], volume: 1000 } }   # also wire (length), shell (area)
```

The centroid is that of the length of an edge or wire, the area of a face or shell, the volume
of a solid; it must be within 1e-5 (1 + |P|) of the point, and exactly one sub-shape must fit
(the measure tells apart those sharing a centroid), else the case fails.

## Assertions

| assertion | meaning |
|---|---|
| `valid: s` | the shape passes validity checking |
| `faults: s` | the number of faulty sub-shapes (the shape included), with a comparison, or with `source: t`: no more than `t` has |
| `exists: name` | the name is bound |
| `volume: s`, `area: s`, `length: s` | with `eq: v` and `rel: r` (default 1e-6) or `abs: a` (both: their sum); `eq: 0` is exact |
| `tolerance: s` | the largest tolerance of the shape's vertices, edges and faces, with a comparison. Tolerances are a kernel's own: ported checks that a tolerance equals the original kernel's value become upper bounds |
| `kind: s` | with `is: compound`, `solid`, `shell`, `face`, ... |
| `closed: s` | each shell of the shape is closed |
| `free_edges: s` | with `bounds: closed` (default) or `open` and a comparison: the number of edges of exactly one face (degenerate ones aside) lying on closed (or open) chains of such edges |
| `free_vertices: s` | the number of vertices met an odd number of times along the edges, with a comparison |
| `count: s` | with `eq: {vertices, edges, wires, faces, shells, solids, compsolids, compounds, shapes}` (or `lt`, `le`, `gt`, `ge` with such a map): distinct sub-shapes, the shape included (`shapes` sums them), one placed twice counting once; with `placements: true` distinct placed sub-shapes; with `occurrences: true` each use (an edge of two faces twice) |
| `centroid: s` | `eq: P`, within `abs` (default 1e-7) per coordinate: the centre of the solids' volume, else of the faces' area, else of the edges' length |

Assertions may also use `lt`, `le`, `gt`, `ge` in place of `eq`.

## Expected outcomes

`expect: not_done` marks a case known to fail for want of a feature; the runner reports it as
expected to fail, and reports it if it passes. `skip: reason` keeps a case from running.
A case fails once it has used `KENTRON_CASE_TIMEOUT` seconds of processor time (40 by default),
or `timeout` seconds where it sets that (from a `{timeout, reason}` entry in
`draw_deviations.json`: a case known to run long and fail); an explicit `KENTRON_CASE_TIMEOUT`
overrides both.

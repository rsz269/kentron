# Known gaps

Major gaps and errors found during development that were judged too complex or too
time-consuming to fix when they were found, and remain to be fixed. Each entry says what
fails, how to reproduce it, the evidence, what was tried, and why it was left. Remove an
entry when it is fixed (the commit fixing it says so).

Smaller differences on the public corpus are tracked by the corpus runs themselves; this
log is for the root causes behind them that need real work.

## Fillets

### Saddle vertices: three concave and three convex edges at one corner

*Found 2026-11-08.*

Where six edges meet at a vertex, three concave and three convex (the inside corner where
three box arms meet in 3D), filleting every edge gives a valid solid with the wrong shape at
every mixed corner of that solid.

- **Repro:** three 6×6 arms along +x, +y and +z from the origin
  (`box(15, 6, 6).translate((7.5, 0, 0))` and the same along y and z), unioned, then
  `edges().fillet(1)`. Its saddle vertex is (3, 3, 3).
- **Evidence:** the volume is 1440.990 against CadQuery's 1441.091; the area is 1003.0
  against 966.2, which CadQuery gets with tori and 8 B-spline corner patches that kentron
  lacks. Against a brute-force rolling ball (the morphological closing then opening of the
  solid by radius 1, on a 0.02 voxel grid round the mixed corner at (3, −3, 3)), 144 of
  60 000 points are wrong on kentron, about 0.15 of volume at that corner, and none on
  CadQuery.
- **Cause:** `fillet_solid_edges_inner` rounds the concave edges first, so that the
  convex ones go round their ends as tori, as OCCT does. After that first stage, the three
  concave fills meet at the saddle vertex only along lines: six short non-manifold edges,
  each bounding four faces. The second stage rejects them ("a fillet edge must bound two
  faces, it bounds 4"), and the whole solid is filleted at once instead, every mixed corner
  mitred.
- **Needed:** a corner patch for saddle vertices (OCCT fills them with B-spline patches),
  or a first stage that leaves the saddle manifold. Several days' work.
- **Not affected:** planar L-shapes and other convex/concave corners without a saddle
  vertex match CadQuery exactly.

## Offsets and shells

### Sliver faces in shells found by distance

*Found 2026-11-08.*

Shells and offsets that fall back on the points-within-the-distance construction
(`offset/within.rs`, used where faces vanish) have the right volume and are valid, but
their faces can come out split into thin slivers near the corners of the cavity.

- **Repro:** `box(20, 20, 20).edges().chamfer(2).faces(">Z").shell(-3)` gives 123 faces
  against 43 (volume 4482.848, exact), and a 60×40×25 block with vertical edges filleted 5
  and top edges filleted 2, `faces("<Z").shell(-3)`, gives 43 faces against 27.
- **Cause:** the cell containing the cavity leaks into the corner regions in thin slivers
  (about 0.015 thick) where the moved faces' sheets nearly meet. The slivers' faces lie on
  the right planes, folded back on themselves, and do not merge.
- **Tried:** deciding tiny cells by the neighbours they share faces with, and growing every
  planar sheet past its face. Neither changed anything.
- **Matters because:** selectors such as `faces(">Z")` on the result can pick up slivers.

## Features

### Extruding up to a curved face that not every ray meets

*Found 2026-11-08 (`anniedoris_synthetic_cadquery/docs_examples/example24.py`).*

`cutBlind(face)` / `extrude(face)` up to a curved face (pykentron's `Solid.dprism`) splits
the profile's prism, extruded far, by the face's whole surface, and keeps the piece touching
the profile. Where some of the profile's rays miss the face, that piece runs on along them
for the whole length, cutting whatever it passes. OCCT cuts only along the rays that meet
the face (the face itself, as trimmed), up to it, and leaves the others alone.

- **Repro:** example24: a slot 2 wide from z −3 to 7 at x = 6, cut up to the spherical
  cavity's face (r 5, trimmed at the box's top, z = 3). kentron removes 41.35 (through the
  shaft above the box), CadQuery 15.43.
- **Evidence:** a Monte Carlo count of the solid's points along the slot's rays that meet
  the face, between x = 6 and it, gives 15.54 ± 0.09.
- **Needed:** the face's shadow on the profile's plane (its edges' projections and its
  silhouette), to keep only the rays meeting it; a face that rays cross twice (as here, a
  sphere) rules out sweeping the face back instead.

### Removing all the fillets of a rounded box

*Found 2026-11-07 (usd-brep's `TestSmRemoveFilletNurbsPlanes`, the `expect: not_done` case
`usd_brep::cpp_fillets::remove_fillets_of_a_nurbs_box`).*

`defeature::remove_faces` of every fillet of a box rounded on all its edges leaves no faces,
where SMLib, OCCT and sampling all give the box back.

- **Cause:** all the fillets are one feature, connected at the corners, and `defeature.rs`
  carries the faces round a feature over it as cells: with the feature everywhere, there is
  nothing round it to carry.
- **Needed:** removing a feature that surrounds the faces it should extend.

## Booleans

### Gyroid lattice cells meeting at a corner

*Found 2026-11-08.*

LatticeQuery's gyroid unit cell is a union of eight thickened plates, some a hundred-
millionth apart. The plates' pairwise unions are valid now (7d15022). The whole cell is
not: where four plates meet at a corner of the cell, their side strips meet 1e-8 apart, and
the union leaves an edge 3.5e-5 long used by seven faces.

- **Repro:** the plates `g000`, `g110`, `g001` and `g011` of `lq/topologies/gyroid.py`
  (see `test_gyroid_plates_a_hundred_millionth_apart` in
  `crates/pykentron/tests/test_against_cadquery.py` for building them), unioned in turn.
  The edge is at (10, 10, 10).
- **Evidence:** it is the same before 7d15022, which only exposed it (the self-crossing
  faces were reported first).
- **Also:** CadQuery's own unions of these plates are wrong: a brute-force count puts the
  first two plates' common part at 0.0024 ± 0.0004, not 0.18.
- **Needed:** coincidence handling for several faces within tolerance of one another at a
  point.

### Sweeps along separately fitted copies of one curve

*Found 2026-10-24; still failing 2026-11-08 (`jhnphm_spiral_lamp/lamp.py`).*

Two sweeps along B-splines fitted separately to the same parametric spiral, over
different ranges, have faces that should coincide but differ by 1e-9 to 1e-4 (the frames'
error times the 30 mm profile). Their union is wrong (61 faces and 2 solids against 37 and
1), and the spiral lamp's later fillet fails on the invalid solid ("a fillet edge must bound
two faces, it bounds 4").

- **Tried:** a geometric test of partial coincidence (grid points with their feet inside
  the other face all within 1e-6). It catches only some pairs, and accepting 1e-4 would be
  unsafe. Reverted.

### Twisted prisms passing through each other with the same section

*Found 2026-11-08 (`mryndzionek_cadquery_plaground/screw_cap_container.py`).*

A knurl tool is the union of two sets of 12 triangular prisms extruded with a twist of +40°
and −40° (`placeSketch(...).twistExtrude`). A prism of one set and one of the other, whose
twists bring them to the same angle at some height, have identical sections there and pass
through each other in a small lens around it. kentron's booleans of such a pair are wrong:
their intersection comes out −6.78 (0.0022 from the same BREPs in Rust), and the union of
the two sets is invalid with 224.4 of volume, against 103.9 for the two sets apart. So the
knurled cap is cut to 252 of volume, invalid, against CadQuery's 4804.

- **Repro:** `knurl()` in the script's `util.py` with height 10, radius 26.75, cut angle
  120, depth 0.5, twist 40 and 12 prisms. The smallest case is two prisms: the one starting
  at −60° and the one at 0°, which coincide at z = 7.5.
- **Evidence:** their section is a closed network of 9 edges about 1 long round the
  coinciding section. CadQuery calls the pair touching (common volume 0, distance 1e-15),
  but its own union of the two sets is empty, so it is no reference either.
- **Needed:** handling of faces that coincide only along a section (a curve of tangency
  sweeping through), on B-spline ruled surfaces.

### Mirrored free-form faces meeting the mirror plane almost square

*Found 2026-11-08.*

A union of a solid with its mirror image, where a free-form face meets the mirror plane at
nearly a right angle, keeps fragments of faces inside the other solid. Their intersection
curve lies in the mirror plane, and there the face and its mirror image cross at only 1e-2 to
1e-3 rad.

- **Repro:** LatticeQuery's `lattice_scripts/hetero-shwarz.py` (one Schwarz-P cell, size
  3.98, thickness 0.9). Its first step is enough: `s000 = schwartz_p_000(0.9, 1.99)` (an
  `interpPlate` face thickened by ±0.45 and unioned) placed at (0.995, 0.995, 0.995), then
  `s000.union(s000.mirror("ZY", (1.99, 0, 0)))`.
- **Evidence:** the union's area is 31.732 on kentron (41 faces), against 31.589 from a
  brute force on kentron's own two halves (a fine mesh of one half, each triangle near the
  plane classified by OCCT against the other half). Classifying every face of the result
  against both halves finds four fragments, 0.036 and 0.034 of area each side, lying up to
  3.6e-3 inside the other half; that accounts for the whole 0.14 surplus. CadQuery's result
  is wrong there too (31.583, with 0.14 of face inside and some missing). The whole cell,
  after four more such unions, comes to 119.288 on kentron against CadQuery's 116.379. The
  volumes are within 0.1%.
- **Cause, as far as found:** the halves' offset faces (A0 and its mirror image B0) cross
  in the plane x = 1.99, along y from 0 to 0.094 near (1.99, y, 1.444). That stretch is
  nearly tangent; it ends at a point where the faces are tangent and a second, slanting
  branch of A0 ∩ B0 crosses it. `march::trace` does trace that branch in the plane: a
  318-point chain at x = 1.99, kept as a section, its tangent stretches at its ends.
  But `boolean::section` of the two faces yields only a 1.85e-5 piece of it at the
  corner. So the curve is lost after marching, where the traced sections are cut to the
  faces. Face A0 is then not split, and `clean_components` joins A0's regions inside
  and outside B0 into one group of 21. Its majority vote (Outside) then overrides the four
  regions classified Inside on their own.
- **Next:** follow chain 318 from `intersect.rs` (`Section::Sampled` from `marched()`)
  through `graph.intersect_faces` to where it is cut down to the faces' common part, with
  the two faces exported alone (`section_probe` on them reproduces it).

### Section edges under-report their tolerance

*Found 2026-10-23; not re-checked since (bug32945 still passes as it stands).*

Fitted B-spline section curves stray from their surfaces between fit points by more than
the tolerance recorded for them (the teapot spout: 5.75e-7 against 1e-7), so a later
boolean may not find them on the surface.

- **Tried:** recording the true deviation by dense sampling. It is correct, but it fails
  OCCT's tolerance checks in bug32945 and buc60555_1, whose real gaps are larger.
- **Needed:** tighter section fitting, so that the true deviation is small.

## pykentron

### Solids from placed shells keep copies of their faces

*Found 2026-11-08 (CadQuery's `tests/test_free_functions.py::test_solid`, the one test of
its suite in `crates/pykentron/tests/cadquery_baseline.txt`).*

`solid(faces, inner_faces, history=...)` with inner faces of placed shells (boxes `moved` to
two places): the history's images do not include the solid's faces from the inner shells,
because those faces are copies, not the inputs.

- **Cause:** kentron's solids hold their shells unplaced (`make_solid` copies a placed shell),
  where OCCT's hold placed ones, so its solid's faces are the placed inputs themselves.
- **Needed:** placed shells in solids (a core invariant of the model), or a history whose
  images are not the inputs' own faces. Not worth it for this alone.

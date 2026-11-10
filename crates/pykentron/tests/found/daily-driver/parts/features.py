# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Reusable mechanical primitives — established convention, authored ONCE.

These are solved-by-practice mechanical features: screw bosses, posts, fillets,
threads, snaps. They are NOT taste and NOT regenerated per part by AI or
re-derived by hand each time. Author each one here from established practice,
then REUSE it across parts. Taste (proportions, vent form, where things go)
stays in params.py + the parametric geometry in parts/*.py. See
docs/design-pipeline.md, "Taste vs. convention".

Convention captured so far:
- boss()       — cylindrical boss for a brass heat-set insert
- screw_post() — post for a socket-head fastener (clearance/pilot bore)

Both stand on a host's horizontal plane, merge into the host, and get a fillet
at the base junction so they're structurally tied and printable (no floating
features). Parts call these instead of inlining boss/post geometry.
"""

import cadquery as cq


def _is_horizontal_circle(edge: "cq.Edge", z: float, radius: float,
                          tol: float = 0.25) -> bool:
    """True for a circular edge of ~`radius` lying on the plane z=`z`."""
    try:
        if edge.geomType() != "CIRCLE":
            return False
        if abs(edge.Center().z - z) > tol:
            return False
        return abs(edge.radius() - radius) < tol
    except Exception:  # noqa: BLE001 — non-circular / degenerate edges
        return False


def _studs(host: cq.Workplane, points, *, floor_z: float, outer_diameter: float,
           bore_diameter: float, height: float, bore_depth: float,
           base_fillet: float) -> cq.Workplane:
    """Shared body for boss()/screw_post().

    Builds cylinders (outer_diameter x height) standing on the plane z=floor_z at
    each (x, y) in `points`, bores each from the top (bore_diameter x bore_depth),
    unions them into `host`, and fillets the base junction circles into the host
    wall (base_fillet). Returns the merged host.
    """
    studs = (
        cq.Workplane("XY")
        .workplane(offset=floor_z)
        .pushPoints(list(points))
        .circle(outer_diameter / 2)
        .extrude(height)
    )
    if bore_depth > 0 and bore_diameter > 0:
        bores = (
            cq.Workplane("XY")
            .workplane(offset=floor_z + height)
            .pushPoints(list(points))
            .circle(bore_diameter / 2)
            .extrude(-bore_depth)
        )
        studs = studs.cut(bores)

    result = host.union(studs)

    if base_fillet and base_fillet > 0:
        base_edges = result.edges().filter(
            lambda e: _is_horizontal_circle(e, floor_z, outer_diameter / 2)
        )
        if base_edges.vals():
            try:
                result = base_edges.fillet(base_fillet)
            except Exception as e:  # noqa: BLE001 — report, don't mask
                print(
                    f"  [warn] features: base fillet skipped ({e}). "
                    "Boss is still merged into the host (tie is just not rounded)."
                )
    return result


def boss(host: cq.Workplane, points, *, floor_z: float, outer_diameter: float,
         bore_diameter: float, height: float, bore_depth: float,
         base_fillet: float) -> cq.Workplane:
    """Heat-set-insert boss(es).

    A cylinder (outer_diameter x height) with a central blind bore
    (bore_diameter x bore_depth) sized for a BRASS HEAT-SET INSERT, standing on
    the host's horizontal plane at z=floor_z at each (x, y) in `points`, merged
    into the host with a base_fillet at the junction. Returns the merged host.

    Convention, reused — not taste. The insert is installed in the open (top) end.
    """
    return _studs(
        host, points, floor_z=floor_z, outer_diameter=outer_diameter,
        bore_diameter=bore_diameter, height=height, bore_depth=bore_depth,
        base_fillet=base_fillet,
    )


def screw_post(host: cq.Workplane, points, *, floor_z: float,
               outer_diameter: float, bore_diameter: float, height: float,
               bore_depth: float, base_fillet: float) -> cq.Workplane:
    """Socket-head fastener post(s).

    Like boss(), but `bore_diameter` is a screw CLEARANCE or PILOT bore (the
    screw passes through or thread-forms into it) rather than an insert bore.
    Pass bore_depth >= height for a through clearance hole. Merged into the host
    with a base_fillet at the junction. Returns the merged host.

    Convention, reused — not taste. Pairs with a clearance hole on the mating part.
    """
    return _studs(
        host, points, floor_z=floor_z, outer_diameter=outer_diameter,
        bore_diameter=bore_diameter, height=height, bore_depth=bore_depth,
        base_fillet=base_fillet,
    )

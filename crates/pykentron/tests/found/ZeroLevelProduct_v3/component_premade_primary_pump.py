"""
Same pump assembly, but with the left/right nozzle symmetry fixed.

Previous bug: the left elbow was built via `mirror('YZ')` then rotated +90°
about Z. The mirror flipped the elbow's handedness (a right-hand bend
became a left-hand bend in the local frame), and the subsequent rotation
then mapped that flipped frame to world coords in a slightly-off way.

Fix: don't mirror at all. Build the elbow once (right-hand bend, inlet
along +Y, curves toward +X). Place the right nozzle with one rotation,
and place the left nozzle with a DIFFERENT rotation that takes the same
unmirrored elbow to the +X-curving-toward-+Y... wait, that's wrong too.

The clean approach: build the elbow once, place the right copy by rotating
the local +Y axis to point along world +X. For the left copy, we want a
mirror image — but we achieve it by mirroring AFTER placement, not before.
So place the right elbow first, then mirror the world-frame result across
the YZ plane to produce the left copy.

Because mirroring is the LAST operation, the final positions of the two
nozzles are exact mirror images of each other — no compounded transform
errors.
"""

from __future__ import annotations
import math
import cadquery as cq

from profile_built_in_2D_sketch import build_2D_sketch
from utils import sweep_profile


# ──────────────────────────────────────────────────────────────────────
# Elbow centerline + sweep (unchanged)
# ──────────────────────────────────────────────────────────────────────
def _elbow_path_wire(
    R_bend: float, arc_deg: float, L_inlet: float, L_leg: float,
    overshoot: float = 0.0,
) -> cq.Wire:
    arc_rad = math.radians(arc_deg)

    P_arc_start = (0.0, L_inlet, 0.0)
    P_arc_end   = (R_bend * (1.0 - math.cos(arc_rad)),
                   L_inlet + R_bend * math.sin(arc_rad),
                   0.0)
    tan_end     = (math.sin(arc_rad), math.cos(arc_rad), 0.0)

    leg_total = L_leg + overshoot
    P_end     = (P_arc_end[0] + leg_total * tan_end[0],
                 P_arc_end[1] + leg_total * tan_end[1],
                 0.0)

    half      = arc_rad / 2.0
    P_arc_mid = (R_bend * (1.0 - math.cos(half)),
                 L_inlet + R_bend * math.sin(half),
                 0.0)

    edges = [
        cq.Edge.makeLine(cq.Vector(0.0, 0.0, 0.0),
                         cq.Vector(*P_arc_start)),
        cq.Edge.makeThreePointArc(cq.Vector(*P_arc_start),
                                  cq.Vector(*P_arc_mid),
                                  cq.Vector(*P_arc_end)),
    ]
    if leg_total > 0:
        edges.append(cq.Edge.makeLine(cq.Vector(*P_arc_end),
                                      cq.Vector(*P_end)))

    return cq.Wire.assembleEdges(edges)


def _build_elbow_outer_inner(
    r_pipe: float, wall_t: float,
    R_bend: float, arc_deg: float,
    L_inlet: float, L_leg: float,
    inner_overshoot: float,
):
    outer_path = _elbow_path_wire(R_bend, arc_deg, L_inlet, L_leg, overshoot=0.0)
    inner_path = _elbow_path_wire(R_bend, arc_deg, L_inlet, L_leg, overshoot=inner_overshoot)

    profile_outer = build_2D_sketch({"obj_type": "circle", "radius": r_pipe},          sketch_plane="XY")
    profile_inner = build_2D_sketch({"obj_type": "circle", "radius": r_pipe - wall_t}, sketch_plane="XY")

    outer = sweep_profile(profile_outer, outer_path, isFrenet=True)
    inner = sweep_profile(profile_inner, inner_path, isFrenet=True)
    return outer, inner


def _place_right_nozzle(elbow: cq.Workplane,
                        barrel_radius: float, overshoot: float,
                        nozzle_z: float) -> cq.Workplane:
    """
    Place the right-hand nozzle: rotate so local +Y → world +X
    (inlet points radially outward to the right), then translate.
    """
    return (elbow
            .rotate((0, 0, 0), (0, 0, 1), -90)
            .translate((barrel_radius - overshoot, 0, nozzle_z)))


# ──────────────────────────────────────────────────────────────────────
# Pump assembly
# ──────────────────────────────────────────────────────────────────────
def create_primary_pump(
    barrel_radius:  float,
    barrel_wall_t:  float,
    barrel_height:  float,
    nozzle_r_pipe:  float,
    nozzle_wall_t:  float,
    nozzle_L_leg:   float,
    nozzle_R_bend:  float,
    nozzle_arc_deg: float,
    nozzle_L_inlet: float,
    nozzle_z:       float,
    flange_width:   float,
    flange_height:  float,
    flange_depth:   float,
    z_bottom:       float = 0.0,
    flange_z_top:   float | None = None,
) -> cq.Workplane:

    if flange_z_top is None:
        flange_z_top = barrel_height - 0.5

    overshoot       = barrel_wall_t  # * 1.5
    inner_overshoot = nozzle_wall_t * 2

    # Barrel
    barrel_outer = cq.Workplane("XY").circle(barrel_radius).extrude(barrel_height)

    # Build the elbow once
    elbow_out, elbow_in = _build_elbow_outer_inner(
        r_pipe=nozzle_r_pipe, wall_t=nozzle_wall_t,
        R_bend=nozzle_R_bend, arc_deg=nozzle_arc_deg,
        L_inlet=nozzle_L_inlet, L_leg=nozzle_L_leg,
        inner_overshoot=inner_overshoot,
    )

    # ─ Right-hand nozzle: place via single rotation + translation ─
    j_right_out = _place_right_nozzle(elbow_out, barrel_radius, overshoot, nozzle_z)
    j_right_in  = _place_right_nozzle(elbow_in,  barrel_radius, overshoot, nozzle_z)

    # ─ Left-hand nozzle: mirror the already-placed right nozzle across YZ ─
    # This is the ONLY mirror, applied as the last step in world coordinates.
    # Result is a guaranteed exact mirror image of the right nozzle.
    j_left_out = j_right_out.mirror("YZ")
    j_left_in  = j_right_in.mirror("YZ")

    envelope = barrel_outer.union(j_right_out).union(j_left_out)

    # Flange on +Y side, near the top
    flange = (cq.Workplane("XY")
              .workplane(offset=flange_z_top - flange_height)
              .moveTo(0, barrel_radius + flange_width / 2 - overshoot)
              .rect(flange_depth, flange_width)
              .extrude(flange_height))
    envelope = envelope.union(flange)

    # Cut bores
    barrel_bore = (cq.Workplane("XY")
                   .workplane(offset=barrel_wall_t)
                   .circle(barrel_radius - barrel_wall_t)
                   .extrude(barrel_height - 2 * barrel_wall_t))

    solid = envelope.cut(barrel_bore).cut(j_right_in).cut(j_left_in)

    #if z_bottom != 0.0:
    #    solid = solid.translate((0, 0, z_bottom))

    return solid.clean()


# ──────────────────────────────────────────────────────────────────────
# Demo
# ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    from ocp_vscode import show

    pump = create_primary_pump(
        barrel_radius  = 1.350 / 2,
        barrel_wall_t  = 0.040,
        barrel_height  = 12.000,
        nozzle_r_pipe  = 0.460 / 2,
        nozzle_wall_t  = 0.025,



        # nozzle_L_leg   = 0.645,    # was 0.600 — +4.5cm to reach the tips
        # nozzle_R_bend  = 0.478,    # was 0.460 — +1.8cm to push outer edge out
        # nozzle_arc_deg = 104.0,    # was 105.0 — 1° less, exits slightly more outward
        # nozzle_L_inlet = 0.050,    # unchanged

        nozzle_L_leg   = 0.600, #0.650,
        nozzle_R_bend  = 0.460,    # was 0.430 then 0.550 — go slightly smaller, not bigger
        nozzle_arc_deg = 105.0,
        nozzle_L_inlet = 0.050, 

        # nozzle_L_leg   = 0.650,    # was 0.650 — taller exits to reach the black outline
        # nozzle_R_bend  = 0.430,
        # nozzle_arc_deg = 105.0,  # was 107
        # nozzle_L_inlet = 0.050,
        nozzle_z       = 0.450,
        flange_width   = 0.548,
        flange_height  = 0.900,
        flange_depth   = 0.500,
    )
    show(pump)























# """
# Parametric primary pump — ZLP V3 (with corrected J-bend elbow geometry).

# Two hollow elbow nozzles attach to the sides of a vertical hollow barrel.
# Each elbow consists of:
#   - straight inlet stub
#   - circular arc bend (sweep angle = arc_deg)
#   - straight exit leg

# The elbow's local frame has the inlet starting at (0,0,0) along +Y,
# curving toward +X for the right-hand version. The whole elbow is then
# rotated −90° about Z so the inlet points along +X (radially outward from
# the barrel center). The left-hand version mirrors across YZ and rotates
# +90° about Z so its inlet points along −X.

# Operation order:
#   1. Build filled barrel + filled elbow outer shells, unioned → envelope.
#   2. Cut the barrel bore and both elbow bores from the envelope.
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_built_in_2D_sketch import build_2D_sketch
# from utils import sweep_profile


# # ──────────────────────────────────────────────────────────────────────
# # Elbow centerline + sweep (matches the standalone elbow script).
# # ──────────────────────────────────────────────────────────────────────
# def _elbow_path_wire(
#     R_bend: float, arc_deg: float, L_inlet: float, L_leg: float,
#     overshoot: float = 0.0,
# ) -> cq.Wire:
#     """
#     Build the elbow centerline as a cq.Wire (line + arc + line).
#     Starts at (0,0,0), initial tangent +Y, curves clockwise toward +X.
#     Optional `overshoot` extends the final straight leg by that amount
#     (used for the bore wire so the cut breaks cleanly through the end).
#     """
#     arc_rad = math.radians(arc_deg)

#     P_arc_start = (0.0, L_inlet, 0.0)
#     P_arc_end   = (R_bend * (1.0 - math.cos(arc_rad)),
#                    L_inlet + R_bend * math.sin(arc_rad),
#                    0.0)
#     tan_end     = (math.sin(arc_rad), math.cos(arc_rad), 0.0)

#     leg_total = L_leg + overshoot
#     P_end     = (P_arc_end[0] + leg_total * tan_end[0],
#                  P_arc_end[1] + leg_total * tan_end[1],
#                  0.0)

#     half      = arc_rad / 2.0
#     P_arc_mid = (R_bend * (1.0 - math.cos(half)),
#                  L_inlet + R_bend * math.sin(half),
#                  0.0)

#     edges = [
#         cq.Edge.makeLine(cq.Vector(0.0, 0.0, 0.0),
#                          cq.Vector(*P_arc_start)),
#         cq.Edge.makeThreePointArc(cq.Vector(*P_arc_start),
#                                   cq.Vector(*P_arc_mid),
#                                   cq.Vector(*P_arc_end)),
#     ]
#     if leg_total > 0:
#         edges.append(cq.Edge.makeLine(cq.Vector(*P_arc_end),
#                                       cq.Vector(*P_end)))

#     return cq.Wire.assembleEdges(edges)


# def _build_elbow_outer_inner(
#     r_pipe: float, wall_t: float,
#     R_bend: float, arc_deg: float,
#     L_inlet: float, L_leg: float,
#     inner_overshoot: float,
# ):
#     """
#     Sweep two solids: outer (r_pipe) and inner (r_pipe - wall_t).
#     Inner path overshoots the outer endpoint by `inner_overshoot` so the
#     bore cut breaks cleanly through the exit face.
#     """
#     outer_path = _elbow_path_wire(R_bend, arc_deg, L_inlet, L_leg, overshoot=0.0)
#     inner_path = _elbow_path_wire(R_bend, arc_deg, L_inlet, L_leg, overshoot=inner_overshoot)

#     profile_outer = build_2D_sketch({"obj_type": "circle", "radius": r_pipe},          sketch_plane="XY")
#     profile_inner = build_2D_sketch({"obj_type": "circle", "radius": r_pipe - wall_t}, sketch_plane="XY")

#     outer = sweep_profile(profile_outer, outer_path, isFrenet=True)
#     inner = sweep_profile(profile_inner, inner_path, isFrenet=True)
#     return outer, inner


# # ──────────────────────────────────────────────────────────────────────
# # Pump assembly
# # ──────────────────────────────────────────────────────────────────────
# def create_primary_pump(
#     # ── Barrel ────────────────────────────────────────────────────────
#     barrel_radius:  float,
#     barrel_wall_t:  float,
#     barrel_height:  float,

#     # ── Nozzle elbows ─────────────────────────────────────────────────
#     nozzle_r_pipe:  float,
#     nozzle_wall_t:  float,
#     nozzle_L_leg:   float,
#     nozzle_R_bend:  float,
#     nozzle_arc_deg: float,    # arc sweep (e.g. 110°)
#     nozzle_L_inlet: float,    # short straight stub between barrel wall and bend
#     nozzle_z:       float,    # axial position of the elbow centerline on the barrel

#     # ── Flange ────────────────────────────────────────────────────────
#     flange_width:   float,
#     flange_height:  float,
#     flange_depth:   float,

#     # ── Positioning ───────────────────────────────────────────────────
#     z_bottom:       float = 0.0,
#     flange_z_top:   float | None = None,
# ) -> cq.Workplane:

#     if flange_z_top is None:
#         flange_z_top = barrel_height - 0.5

#     # how far the elbow inlet pokes into the barrel wall for a clean union
#     overshoot       = barrel_wall_t * 1.5
#     # how far the elbow bore pokes past the outer end face for a clean cut
#     inner_overshoot = nozzle_wall_t * 2

#     # ── 1. FILLED ENVELOPE ──────────────────────────────────────────────
#     barrel_outer = cq.Workplane("XY").circle(barrel_radius).extrude(barrel_height)

#     # Build the elbow once in its local frame (inlet along +Y, curves to +X).
#     elbow_out, elbow_in = _build_elbow_outer_inner(
#         r_pipe=nozzle_r_pipe, wall_t=nozzle_wall_t,
#         R_bend=nozzle_R_bend, arc_deg=nozzle_arc_deg,
#         L_inlet=nozzle_L_inlet, L_leg=nozzle_L_leg,
#         inner_overshoot=inner_overshoot,
#     )

#     # ─ Right-hand nozzle ─
#     # Rotate the elbow −90° about Z so its inlet (originally +Y) points along +X,
#     # then translate so the inlet start sits at the barrel wall on the +X side.
#     j_right_out = (elbow_out
#                    .rotate((0, 0, 0), (0, 0, 1), -90)
#                    .translate((barrel_radius - overshoot, 0, nozzle_z)))
#     j_right_in  = (elbow_in
#                    .rotate((0, 0, 0), (0, 0, 1), -90)
#                    .translate((barrel_radius - overshoot, 0, nozzle_z)))

#     # ─ Left-hand nozzle ─
#     # Mirror across YZ first (curves toward −X), then rotate +90° about Z so
#     # its inlet points along −X, then translate to the −X side of the barrel.
#     j_left_out = (elbow_out.mirror("YZ")
#                   .rotate((0, 0, 0), (0, 0, 1), 90)
#                   .translate((-(barrel_radius - overshoot), 0, nozzle_z)))
#     j_left_in  = (elbow_in.mirror("YZ")
#                   .rotate((0, 0, 0), (0, 0, 1), 90)
#                   .translate((-(barrel_radius - overshoot), 0, nozzle_z)))

#     envelope = barrel_outer.union(j_right_out).union(j_left_out)

#     # ─ Flange on the +Y side, near the top ─
#     flange = (cq.Workplane("XY")
#               .workplane(offset=flange_z_top - flange_height)
#               .moveTo(0, barrel_radius + flange_width / 2 - overshoot)
#               .rect(flange_depth, flange_width)
#               .extrude(flange_height))
#     envelope = envelope.union(flange)

#     # ── 2. CUT BORES ────────────────────────────────────────────────────
#     barrel_bore = (cq.Workplane("XY")
#                    .workplane(offset=barrel_wall_t)
#                    .circle(barrel_radius - barrel_wall_t)
#                    .extrude(barrel_height - 2 * barrel_wall_t))

#     solid = envelope.cut(barrel_bore).cut(j_right_in).cut(j_left_in)

#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid.clean()


# # ──────────────────────────────────────────────────────────────────────
# # Demo
# # ──────────────────────────────────────────────────────────────────────
# if __name__ == "__main__":
#     from ocp_vscode import show

#     pump = create_primary_pump(
#         # barrel
#         barrel_radius  = 1.350 / 2,
#         barrel_wall_t  = 0.040,
#         barrel_height  = 12.000,
#         # nozzles — matches the working standalone elbow
#         nozzle_r_pipe  = 0.460 / 2,
#         nozzle_wall_t  = 0.025,
#         nozzle_L_leg   = 0.650,    # was 0.500 — longer exits
#         nozzle_R_bend  = 0.430,    # original
#         nozzle_arc_deg = 107.0,    # was 110 — +3° to bring exits closer to vertical
#         nozzle_L_inlet = 0.050,    # original
#         # nozzle_L_leg   = 0.500,
#         # nozzle_R_bend  = 0.430,
#         # nozzle_arc_deg = 110.0,
#         # nozzle_L_inlet = 0.050,
#         nozzle_z       = 0.450,
#         # flange
#         flange_width   = 0.548,
#         flange_height  = 0.900,
#         flange_depth   = 0.500,
#     )
#     show(pump)



























# # GOOD BUT WITH WRONG BEND GEOMETRIES

# """
# Parametric primary pump (simplified) — ZLP V3.

# Two hollow J-bend nozzles inserted into a vertical hollow barrel.
# The nozzle bores connect to the barrel interior so fluid can flow.

# All geometry dimensions are required — no defaults. This follows the same
# convention as create_reactor_vessel and create_strongback, keeping
# reactor-specific values in the assembly file rather than the component builder.

# Operation order:
#   1. Build solid (filled) barrel + filled J-bend outer shells, unioned.
#      This gives the full external solid envelope of the pump.
#   2. Cut the barrel inner bore + the J-bend inner bores from this envelope.
#      What remains is the wall everywhere, with all fluid paths hollow.
# """

# from __future__ import annotations
# import cadquery as cq
# from profile_built_in_2D_sketch import build_2D_sketch
# from utils import sweep_profile


# def _jbend_outer_and_inner(
#     r_pipe: float, wall_t: float, L_leg: float, R_bend: float, clockwise: bool,
#     arc_deg: float = 90.0,
#     inner_stub: float = 0.0,
# ):
#     """
#     Returns (outer_solid, inner_solid) for a J-bend along an XY-plane path.
#       arc_deg: total sweep of the arc in degrees (>90 makes the leg angle inward)
#       outer_solid / inner_solid: filled cylinders along the J path (and inner stub)
#     """
#     import math
#     if clockwise:
#         sign = +1
#     else:
#         sign = -1
#     cx = sign * R_bend
#     cy = L_leg

#     a = math.radians(arc_deg)
#     a_eff = -a if clockwise else a

#     rx0 = -sign * R_bend
#     ry0 = 0
#     rx1 = rx0 * math.cos(a_eff) - ry0 * math.sin(a_eff)
#     ry1 = rx0 * math.sin(a_eff) + ry0 * math.cos(a_eff)
#     end = (cx + rx1, cy + ry1)

#     if clockwise:
#         tip_dir = (ry1, -rx1)
#     else:
#         tip_dir = (-ry1, rx1)
#     tlen = math.hypot(*tip_dir)
#     tip_dir = (tip_dir[0] / tlen, tip_dir[1] / tlen)

#     rad = sign * R_bend

#     outer_path = (
#         cq.Workplane("XY")
#         .moveTo(0, 0)
#         .lineTo(0, L_leg)
#         .radiusArc(end, rad)
#         .wire()
#         .val()
#     )

#     if inner_stub > 0:
#         ext_end = (end[0] + tip_dir[0] * inner_stub,
#                    end[1] + tip_dir[1] * inner_stub)
#         inner_path = (
#             cq.Workplane("XY")
#             .moveTo(0, 0)
#             .lineTo(0, L_leg)
#             .radiusArc(end, rad)
#             .lineTo(*ext_end)
#             .wire()
#             .val()
#         )
#     else:
#         inner_path = outer_path

#     profile_outer = build_2D_sketch({"obj_type": "circle", "radius": r_pipe},          sketch_plane="XY")
#     profile_inner = build_2D_sketch({"obj_type": "circle", "radius": r_pipe - wall_t}, sketch_plane="XY")

#     outer = sweep_profile(profile_outer, outer_path, isFrenet=True)
#     inner = sweep_profile(profile_inner, inner_path, isFrenet=True)
#     return outer, inner


# def create_primary_pump(
#     # ── Barrel (required) ─────────────────────────────────────────────────
#     barrel_radius:  float,
#     barrel_wall_t:  float,
#     barrel_height:  float,

#     # ── Nozzle J-bends (required) ─────────────────────────────────────────
#     nozzle_r_pipe:  float,
#     nozzle_wall_t:  float,
#     nozzle_L_leg:   float,
#     nozzle_R_bend:  float,
#     nozzle_arc_deg: float,   # >90° makes legs angle inward
#     nozzle_z:       float,   # axial position of nozzle centre on barrel

#     # ── Flange (required) ─────────────────────────────────────────────────
#     flange_width:   float,   # radial protrusion
#     flange_height:  float,   # axial height
#     flange_depth:   float,   # tangential depth

#     # ── Positioning ───────────────────────────────────────────────────────
#     z_bottom:       float = 0.0,
#     flange_z_top:   float | None = None,   # defaults to barrel_height - 0.5
# ) -> cq.Workplane:
#     """
#     Build a primary pump: hollow barrel with two J-bend nozzles and a flange.

#     Args:
#         barrel_radius:  outer radius of the barrel [m]
#         barrel_wall_t:  barrel wall thickness [m]
#         barrel_height:  total axial height of the barrel [m]
#         nozzle_r_pipe:  outer radius of each nozzle pipe [m]
#         nozzle_wall_t:  nozzle wall thickness [m]
#         nozzle_L_leg:   straight leg length of each J-bend [m]
#         nozzle_R_bend:  bend radius of each J-bend [m]
#         nozzle_arc_deg: arc sweep in degrees (90 = right angle, >90 = legs angle inward)
#         nozzle_z:       z of nozzle connection on barrel (above barrel base) [m]
#         flange_width:   radial protrusion of rectangular flange [m]
#         flange_height:  axial height of rectangular flange [m]
#         flange_depth:   tangential depth of rectangular flange [m]
#         z_bottom:       translate so barrel base sits at this z [m]
#         flange_z_top:   z of flange top face; defaults to barrel_height - 0.5 [m]

#     Returns:
#         cq.Workplane — single solid
#     """
#     if flange_z_top is None:
#         flange_z_top = barrel_height - 0.5

#     overshoot = barrel_wall_t * 1.5

#     # ── 1. Build the FILLED OUTER ENVELOPE (no bores yet) ────────────────
#     barrel_outer = cq.Workplane("XY").circle(barrel_radius).extrude(barrel_height)

#     inner_stub = barrel_wall_t * 4
#     j1_out, j1_in = _jbend_outer_and_inner(nozzle_r_pipe, nozzle_wall_t,
#                                            nozzle_L_leg, nozzle_R_bend,
#                                            clockwise=False,
#                                            arc_deg=nozzle_arc_deg,
#                                            inner_stub=inner_stub)
#     shift_r = (barrel_radius + nozzle_R_bend - overshoot,
#                -(nozzle_L_leg + nozzle_R_bend), nozzle_z)
#     j1_out = j1_out.translate(shift_r)
#     j1_in  = j1_in.translate(shift_r)

#     j2_out, j2_in = _jbend_outer_and_inner(nozzle_r_pipe, nozzle_wall_t,
#                                            nozzle_L_leg, nozzle_R_bend,
#                                            clockwise=True,
#                                            arc_deg=nozzle_arc_deg,
#                                            inner_stub=inner_stub)
#     shift_l = (-(barrel_radius + nozzle_R_bend - overshoot),
#                -(nozzle_L_leg + nozzle_R_bend), nozzle_z)
#     j2_out = j2_out.translate(shift_l)
#     j2_in  = j2_in.translate(shift_l)

#     envelope = barrel_outer.union(j1_out).union(j2_out)

#     flange = (cq.Workplane("XY")
#               .workplane(offset=flange_z_top - flange_height)
#               .moveTo(0, barrel_radius + flange_width / 2 - overshoot)
#               .rect(flange_depth, flange_width)
#               .extrude(flange_height))
#     envelope = envelope.union(flange)

#     # ── 2. Cut all bores from the envelope ───────────────────────────────
#     # barrel_bore = cq.Workplane("XY").circle(barrel_radius - barrel_wall_t).extrude(barrel_height)
#     barrel_bore = (
#     cq.Workplane("XY")
#     .workplane(offset=barrel_wall_t)
#     .circle(barrel_radius - barrel_wall_t)
#     .extrude(barrel_height - 2 * barrel_wall_t)
#     )
#     solid = envelope.cut(barrel_bore).cut(j1_in).cut(j2_in)

#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid.clean()


# if __name__ == "__main__":
#     from ocp_vscode import show
#     pump = create_primary_pump(
#         barrel_radius  = 1.350 / 2,
#         barrel_wall_t  = 0.040,
#         barrel_height  = 12.000,
#         nozzle_r_pipe  = 0.460 / 2,
#         nozzle_wall_t  = 0.025,
#         nozzle_L_leg   = 0.800,
#         nozzle_R_bend  = 0.800,
#         nozzle_arc_deg = 112.5,
#         nozzle_z       = 0.450,
#         flange_width   = 0.548,
#         flange_height  = 0.900,
#         flange_depth   = 0.500,
#     )
#     show(pump)

















































"""
Parametric primary pump (simplified) — ZLP V3.

Two hollow J-bend nozzles inserted into a vertical hollow barrel.
The nozzle bores connect to the barrel interior so fluid can flow.

Operation order:
  1. Build solid (filled) barrel + filled J-bend outer shells, unioned.
     This gives the full external solid envelope of the pump.
  2. Cut the barrel inner bore + the J-bend inner bores from this envelope.
     What remains is the wall everywhere, with all fluid paths hollow.
"""
"""
from __future__ import annotations
import cadquery as cq
from profile_built_in_2D_sketch import build_2D_sketch
from utils import sweep_profile


def _jbend_outer_and_inner(
    r_pipe: float, wall_t: float, L_leg: float, R_bend: float, clockwise: bool,
    arc_deg: float = 90.0,
    inner_stub: float = 0.0,
):
    """
    # Returns (outer_solid, inner_solid) for a J-bend along an XY-plane path.
    #  arc_deg: total sweep of the arc in degrees (>90 makes the leg angle inward)
    #  outer_solid / inner_solid: filled cylinders along the J path (and inner stub)
"""
    import math
    # Arc starts at (0, L_leg) tangent to +Y. Sweeps `arc_deg` clockwise (or CCW).
    # Centre of arc at (sign*R_bend, L_leg).
    # End point of arc: rotate (0, L_leg) by arc_deg around centre.
    if clockwise:
        sign = +1
    else:
        sign = -1
    cx = sign * R_bend
    cy = L_leg

    # Start point relative to centre: (-sign*R_bend, 0)
    # Rotate by arc_deg in the appropriate direction
    a = math.radians(arc_deg)
    if clockwise:
        # CW rotation: (x,y) -> (x cos - y sin, x sin + y cos)? Use -a for CW
        a_eff = -a
    else:
        a_eff = a

    rx0 = -sign * R_bend
    ry0 = 0
    rx1 = rx0 * math.cos(a_eff) - ry0 * math.sin(a_eff)
    ry1 = rx0 * math.sin(a_eff) + ry0 * math.cos(a_eff)
    end = (cx + rx1, cy + ry1)

    # Tangent at arc end: perpendicular to radius (rx1, ry1), rotated 90° in sweep direction
    if clockwise:
        tip_dir = (ry1, -rx1)
    else:
        tip_dir = (-ry1, rx1)
    # Normalise
    tlen = math.hypot(*tip_dir)
    tip_dir = (tip_dir[0] / tlen, tip_dir[1] / tlen)

    # radiusArc accepts the endpoint and signed radius (sign chooses arc side)
    rad = sign * R_bend

    outer_path = (
        cq.Workplane("XY")
        .moveTo(0, 0)
        .lineTo(0, L_leg)
        .radiusArc(end, rad)
        .wire()
        .val()
    )

    if inner_stub > 0:
        ext_end = (end[0] + tip_dir[0] * inner_stub,
                   end[1] + tip_dir[1] * inner_stub)
        inner_path = (
            cq.Workplane("XY")
            .moveTo(0, 0)
            .lineTo(0, L_leg)
            .radiusArc(end, rad)
            .lineTo(*ext_end)
            .wire()
            .val()
        )
    else:
        inner_path = outer_path

    profile_outer = build_2D_sketch({"obj_type": "circle", "radius": r_pipe},          sketch_plane="XY")
    profile_inner = build_2D_sketch({"obj_type": "circle", "radius": r_pipe - wall_t}, sketch_plane="XY")

    outer = sweep_profile(profile_outer, outer_path, isFrenet=True)
    inner = sweep_profile(profile_inner, inner_path, isFrenet=True)
    return outer, inner


def create_primary_pump(
    barrel_radius:  float = 1.350 / 2,
    barrel_wall_t:  float = 0.040,
    barrel_height:  float = 12.000,

    nozzle_r_pipe:  float = 0.460 / 2,
    nozzle_wall_t:  float = 0.025,
    nozzle_L_leg:   float = 0.800,
    nozzle_R_bend:  float = 0.800,
    nozzle_arc_deg: float = 112.5,    # >90° makes legs angle inward (135° between legs)
    nozzle_z:       float = 0.450,

    # Top-right rectangular flange (from drawing)
    flange_width:   float = 0.548,    # 548 mm — radial protrusion
    flange_height:  float = 0.900,    # 900 mm — axial height
    flange_depth:   float = 0.500,    # tangential depth (estimated; not given on drawing)
    flange_z_top:   float = 11.500,   # top of flange aligns with top of barrel area

    z_bottom:       float = 0.0,
) -> cq.Workplane:
    overshoot = barrel_wall_t * 1.5

    # ── 1. Build the FILLED OUTER ENVELOPE (no bores yet) ────────────────
    barrel_outer = cq.Workplane("XY").circle(barrel_radius).extrude(barrel_height)

    # Right J outer & inner sweeps. Inner sweep is extended past the arc tip
    # by `inner_stub` so the cutter cleanly punches through the barrel wall.
    inner_stub = barrel_wall_t * 4
    j1_out, j1_in = _jbend_outer_and_inner(nozzle_r_pipe, nozzle_wall_t,
                                           nozzle_L_leg, nozzle_R_bend,
                                           clockwise=False,
                                           arc_deg=nozzle_arc_deg,
                                           inner_stub=inner_stub)
    shift_r = (barrel_radius + nozzle_R_bend - overshoot,
               -(nozzle_L_leg + nozzle_R_bend), nozzle_z)
    j1_out = j1_out.translate(shift_r)
    j1_in  = j1_in.translate(shift_r)

    # Left J outer & inner sweeps
    j2_out, j2_in = _jbend_outer_and_inner(nozzle_r_pipe, nozzle_wall_t,
                                           nozzle_L_leg, nozzle_R_bend,
                                           clockwise=True,
                                           arc_deg=nozzle_arc_deg,
                                           inner_stub=inner_stub)
    shift_l = (-(barrel_radius + nozzle_R_bend - overshoot),
               -(nozzle_L_leg + nozzle_R_bend), nozzle_z)
    j2_out = j2_out.translate(shift_l)
    j2_in  = j2_in.translate(shift_l)

    # Union all outer solids → the complete external envelope
    envelope = barrel_outer.union(j1_out).union(j2_out)

    # ── Top-right rectangular flange (from drawing, near top of barrel) ─
    # Top view shows it on the opposite side of the barrel from the two U-bends.
    # Push slightly inward (overshoot) so it overlaps the barrel wall instead of
    # just touching at a single point.
    flange = (cq.Workplane("XY")
              .workplane(offset=flange_z_top - flange_height)
              .moveTo(0, barrel_radius + flange_width / 2 - overshoot)
              .rect(flange_depth, flange_width)
              .extrude(flange_height))
    envelope = envelope.union(flange)

    # ── 2. Cut all bores from the envelope ───────────────────────────────
    # Barrel inner bore: full-height vertical cylinder
    barrel_bore = cq.Workplane("XY").circle(barrel_radius - barrel_wall_t).extrude(barrel_height)

    solid = envelope.cut(barrel_bore).cut(j1_in).cut(j2_in)

    if z_bottom != 0.0:
        solid = solid.translate((0, 0, z_bottom))

    return solid.clean()


if __name__ == "__main__":
    from ocp_vscode import show
    pump = create_primary_pump()
    show(pump)
"""





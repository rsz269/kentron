# SPDX-FileCopyrightText: 2026 Jamey Warren
# SPDX-License-Identifier: MIT

"""
Assembly — one side of the Daily Driver in its correct relationship (v0.3).

Chain: cup → baffle (front-mounted) → fork-yoke (pivoted to the cup) → slider
(rides the fork's adjustment post) → bow (REFERENCE body). The cup/baffle/yoke/slider
poses
are the real kinematic chain; the bow is posed representatively so one end sits in
the slider channel and the arc sweeps over toward the other (mirror) side — the
full head-size kinematics are TBD from the measured bow (ESTIMATES, see params).

View it:  python show.py        (OCP CAD Viewer)
Export:   build.py writes output/assembly.step
"""

import math

import cadquery as cq

from params import P
from parts.cup import make_cup
from parts.baffle import make_baffle
from parts.yoke import make_yoke
from parts.slider import make_slider
from parts.bow import make_bow
from parts.headband_pad import make_headband_pad
from parts.driver import make_driver
from parts.driver_clamp import make_driver_clamp
from parts.earpad import make_earpad
from parts.headband_clamp import make_headband_clamp


# Sub-assembly groups for the manual's interactive parts viewer. The node NAMES
# below are a PUBLIC CONTRACT: the website viewer toggles/isolates parts by them
# (they match the asm.add(name=...) calls in make_assembly). Renaming a part there
# means updating this too (and the manual). build.py emits this to
# docs/models/daily-driver.groups.json next to the GLB; the viewer fetches it.
SUBASSEMBLIES = {
    "groups": [
        {"id": "earcup", "label": "Earcup",
         "nodes": ["cup_R", "cup_L", "baffle_R", "baffle_L",
                   "driver_R", "driver_L", "driver_clamp_R", "driver_clamp_L"]},
        {"id": "earpad", "label": "Earpads",
         "nodes": ["earpad_R", "earpad_L"]},
        {"id": "acoustic", "label": "Felt + gasket + paper",
         "nodes": ["damping_R", "damping_L", "gasket_R", "gasket_L",
                   "paper_R", "paper_L"]},
        {"id": "gimbal", "label": "Gimbal",
         "nodes": ["yoke_R", "yoke_L", "yoke_rod_R", "yoke_rod_L",
                   "insert_p_R", "insert_p_L", "insert_m_R",
                   "insert_m_L", "screw_p_R", "screw_p_L", "screw_m_R", "screw_m_L"]},
        {"id": "headband", "label": "Headband",
         "nodes": ["bow_ref", "slider_R", "slider_L", "thumbscrew_R", "thumbscrew_L",
                   "slider_shoe_R", "slider_shoe_L", "headband_clamp_R", "headband_clamp_L"]},
        {"id": "headband_pad", "label": "Headband pad",
         "nodes": ["headband_pad"]},
        {"id": "head", "label": "Reference head", "nodes": ["head_ref"]},
    ],
    "bought": ["bow_ref", "earpad_R", "earpad_L", "paper_R", "paper_L"],
    # The reference head is translucent worn-fit CONTEXT: the viewer shows it OFF by default and holds
    # it OUT of the explode motion (context, not a part). The whole POSE is fitted to it per size — the
    # viewer swaps GLBs (daily-driver-{s,m,l}.glb) to re-fit, rather than nudging one group. Public contract.
    "reference_context": ["head_ref"],
    "head_sizes": ["s", "m", "l"],   # the viewer offers these; each loads daily-driver-<size>.glb
}


def make_assembly(worn_head: str = "m") -> cq.Assembly:
    """Both ears + the shared headband, posed like a WORN headphone.

    Head frame (global): X = inter-ear (right ear at +Xe, pad facing −X inward),
    Z = up (bow over the crown), Y = front-back (the tilt-pivot axis). The cup and
    yoke are co-designed with pad ∥ arch (both +Z), which is why earlier the cup sat
    "face-up". Here they're mounted at the correct 90° relative clocking:
      • cup-group  T_cup : pad (cup +Z) → −X (inward), pivot (cup ±X) → ±Y,
                           up (cup +Y) → +Z;  pivot centre → (Xe, 0, 0).
      • yoke-group T_yoke: eyes (yoke ±X) → ±Y, arch (yoke +Z) stays +Z (up);
                           pivot centre → (Xe, 0, 0)  (90° clocked about Y vs the cup).
    The two pivots coincide at (±Xe, 0, 0); the LEFT ear is the mirror across x=0.
    The bow + crown pad are shared parts arcing between the two sliders, posed
    FLEXED to the worn radius (bow_worn_radius); ear spacing = where the flexed
    band's ends land (~156 mm). The bow's relaxed dims are MEASURED off the Beyer
    part; the worn radius (head fit) is still ESTIMATE (see params/bow).

    HEADBAND JUNCTION (v0.7, offset-outer): the metal band rides INSIDE (head-side) and
    the slider's post-bore tube rides OUTSIDE it, so the cups step outboard by
    (barrel R + clamp-plate depth) to keep the post under the band's clamp. Band end →
    recess in the slider's clamp plate → cover (inner) + 2 screws; rib through the bow
    channel; thumbscrew locks the post in the tube.
    """
    CHARCOAL = cq.Color(0.30, 0.32, 0.35)
    ORANGE = cq.Color(0.92, 0.45, 0.10)
    YOKE_C = cq.Color(0.36, 0.38, 0.42)
    SLIDER_C = cq.Color(0.42, 0.44, 0.48)
    STEEL = cq.Color(0.75, 0.76, 0.78)
    BRASS = cq.Color(0.80, 0.68, 0.30)
    SCREW_C = cq.Color(0.55, 0.57, 0.60)
    PAD_C = cq.Color(0.13, 0.13, 0.15)   # near-black foam/velour
    DRIVER_C = cq.Color(0.10, 0.10, 0.12)  # driver mockup (black)

    pbz = P.pivot_boss_z
    # HEAD-DRIVEN worn pose (per `worn_head` ∈ {s,m,l}): the headphone is FITTED to the chosen head.
    # The FULL earpad sits on the head, COMPRESSING into it by (earpad_depth − earpad_worn_depth) ≈ 10 mm
    # = the worn CLAMP contact (maker: show the full pad with real clamping force). earpad_worn_depth is the
    # seated cup-front→head gap, so it sets the cup spacing.
    # WIDTH is then taken up by the bow SPRING (its flex), HEIGHT by the slider POST (travel) — the
    # maker's mechanism. The bow is a PURE ARC (no bends). For its end to lie FLAT on the slider clamp
    # face (maker: "the ends of the band are parallel to the slider face, no bends other than the arc"),
    # the whole yoke+post+slider CHAIN TILTS by the band's end-tangent angle `psi` about the cup pivot —
    # the cup swivels to stay flat on the ear. Tilting also pulls the band ends INBOARD, so the band
    # stays TIGHT (near its 5 in rest) instead of flexing flat/extended (the maker's other note).
    eh_map = {"s": P.head_s_ear_half, "m": P.head_ref_ear_half, "l": P.head_l_ear_half}
    eh = eh_map.get(worn_head, P.head_ref_ear_half)
    crown = P.head_ref_z + P.head_ref_height_half * (eh / P.head_ref_ear_half)
    # The metal band rides a PAD-THICKNESS above the crown (the headband pad fills the gap, resting on the
    # head) — so the band APEX targets crown + pad. This both seats the pad correctly (it was being buried)
    # and lifts the slider → the yoke posts show real EXTENSION (cups hang on visible rod), per the maker.
    apex_target = crown + P.headband_pad_thickness

    # Cup pivot x = head ear surface (eh) + seated pad contact (earpad_worn_depth) + the cup's FRONT→PIVOT
    # offset (= cup_total_height − pbz; the pivot now sits forward of mid by pivot_boss_forward, so this is
    # < pbz and pulls the whole junction inboard). Cup stays ON the ear; the full pad compresses in.
    Xe_cup = eh + P.earpad_worn_depth + (P.cup_total_height - pbz)

    # ---- Solve the TILT (psi) self-consistently --------------------------------------------------
    # A pure arc of the measured developed length L has end-tangent angle == its end PARAM angle, so
    # tying the post tilt to the band tangent gives  arc_worn = 180 − 2·psi  and  R = L/(2·rad(90−psi)).
    # The chain tilt rotates the slider clamp hole by psi about the cup pivot; `psi` is the value where
    # the bow's prong-tip hole (inset from the tip) lands on that rotated clamp hole with the band APEX
    # on the crown. seat_y = the recess-floor depth (clamp-hole x-offset from the post centre).
    L = P.bow_developed_length
    seat_y = (-P.slider_collar_diameter / 2 - P.slider_clamp_standoff) + (P.bow_thickness + 0.4)
    post_top = P.yoke_fork_height + 4 + P.yoke_post_length
    post_base = P.yoke_fork_height + 4                      # hub top = the barrel's bottom stop
    sz_hi = post_top - P.slider_collar_height / 2           # barrel at the post top (extended)
    sz_lo = post_base + P.slider_collar_height / 2          # barrel at the hub (retracted)

    def _R_of_psi(psi_deg):
        return L * 90.0 / (math.pi * (90.0 - psi_deg))

    def _slider_z_ideal(psi_deg, R, zh_):
        p_ = math.radians(psi_deg)
        return (zh_ + apex_target - R - seat_y * math.sin(p_)) / math.cos(p_)

    def _resid(psi_deg):                                    # band prong-hole x  −  rotated clamp-hole x
        p_ = math.radians(psi_deg)
        R = _R_of_psi(psi_deg)
        ah_ = math.radians(psi_deg + math.degrees(P.bow_endtab_hole_inset / R))
        zh_ = R * math.sin(ah_)
        sl_ = _slider_z_ideal(psi_deg, R, zh_)
        return R * math.cos(ah_) - (Xe_cup + seat_y * math.cos(p_) - sl_ * math.sin(p_))

    lo, hi = 1.0, 80.0
    flo = _resid(lo)
    for _ in range(120):                                    # bisection on the sign change
        mid = (lo + hi) / 2
        fm = _resid(mid)
        if (flo < 0) == (fm < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    psi = (lo + hi) / 2
    R_worn = _R_of_psi(psi)
    arc_worn = 180.0 - 2.0 * psi
    p = math.radians(psi)
    a_hole = math.radians(psi + math.degrees(P.bow_endtab_hole_inset / R_worn))
    zh = R_worn * math.sin(a_hole)                          # prong-hole height in the bow frame
    slider_z = max(sz_lo, min(sz_hi, _slider_z_ideal(psi, R_worn, zh)))   # post extension (height)

    def T_cup(w):    # pad → −X, pivot → ±Y, up → +Z; pivot centre → (Xe_cup,0,0)
        return (w.rotate((0, 0, 0), (0, 1, 0), -90)
                 .rotate((0, 0, 0), (1, 0, 0), 90)
                 .translate((Xe_cup + pbz, 0, 0)))

    def T_yoke(w):   # eyes → ±Y, arch stays +Z (up); pivot centre → (Xe_cup,0,0)
        return w.rotate((0, 0, 0), (0, 0, 1), -90).translate((Xe_cup, 0, 0))

    def tilt(w):     # tilt the post chain inboard by psi about the cup pivot, so the slider clamp
        return w.rotate((Xe_cup, 0, 0), (Xe_cup, 1, 0), -psi)   # face is parallel to the band end

    def mirror_L(w):  # right ear → left ear (true mirror across the head centre)
        return w.mirror("YZ", (0, 0, 0))

    # ---- RIGHT ear ----
    cup = T_cup(make_cup())
    baffle = T_cup(make_baffle().translate((0, 0, P.baffle_seat_z)))
    # Driver (mockup) seated in the baffle's back recess, firing forward; the clamp
    # ring behind it retains the rear rim. ledge = where the flange seats; rear_rim =
    # the driver back where the clamp shoulder bears.
    ledge_z = P.baffle_seat_z + P.driver_recess_depth
    rear_rim_z = ledge_z - P.driver_body_depth
    driver = T_cup(make_driver().translate((0, 0, ledge_z)))
    driver_clamp = T_cup(make_driver_clamp().translate((0, 0, rear_rim_z)))
    # Earpad (mockup) on the cup front rim, ear opening facing the head (cup +Z → −X). Shown at FULL
    # depth; the cup is placed so it compresses ~10 mm into the head = the worn clamp contact.
    earpad = T_cup(make_earpad(P.earpad_depth).translate((0, 0, P.cup_total_height)))
    # Yoke + post + slider all TILT (post chain) so the slider clamp face is parallel to the band.
    yoke = tilt(T_yoke(make_yoke()))
    # Bought Ø6 adjustment ROD — epoxy-bonded into the fork socket, rising as the post. Built from
    # its socket floor (z=0); shift up so it seats in the socket, then ride (tilted) with the yoke.
    # The post is a bought ISO 7379 shoulder screw: M5 thread into the fork boss, Ø6 shoulder is the
    # post, head is the top stop. Built from z=0 = the shoulder seat (boss top); ride with the yoke.
    from parts.yoke_rod import make_yoke_rod
    rod = tilt(T_yoke(make_yoke_rod().translate((0, 0, P.yoke_fork_height + 4))))
    # Slider rides the (tilted) post at slider_z (the height, solved above); the post slides + swivels
    # the full barrel height. slider_z = the clamp/barrel centre along the tilted post.
    slider = tilt(T_yoke(make_slider().translate((0, 0, slider_z))))
    # Pressure SHOE — rides in the slider's +Y pocket, saddle cradling the post (the thumbscrew
    # presses it, not the post). Built at the origin, shifted +Y so its saddle is post-coaxial.
    from parts.slider_shoe import make_slider_shoe, shoe_offset_y
    shoe = tilt(T_yoke(make_slider_shoe().translate((0, shoe_offset_y(), slider_z))))

    # ---- Shared headband: bow + crown pad, a PURE ARC arcing between the two sliders ----
    # The bow is symmetric about x=0 (apex on the crown). Its end tangent already matches the post
    # tilt, so a pure Z-shift lands the prong-tip HOLE on the (tilted) slider clamp hole: the hole
    # rides at global z = Hz = seat_y·sin(psi) + slider_z·cos(psi); the bow's hole sits at zh.
    _bow0 = make_bow(radius=R_worn, arc_degrees=arc_worn)
    Hz = seat_y * math.sin(p) + slider_z * math.cos(p)
    bow_xf = (0, 0, Hz - zh)
    bow = _bow0.translate(bow_xf)                            # the flexed bow, posed
    pad = make_headband_pad(radius=R_worn, arc_degrees=arc_worn).translate(bow_xf)
    # Headband CLAMP cover (INNER head-side piece) — built in the slider frame, posed (tilted) with
    # the slider so its bolts/slot align with the slider's clamp.
    cover = tilt(T_yoke(make_headband_clamp().translate((0, 0, slider_z))))

    asm = cq.Assembly(name="daily_driver")
    for nm, solid, col in (("cup", cup, CHARCOAL), ("baffle", baffle, ORANGE),
                           ("driver", driver, DRIVER_C), ("driver_clamp", driver_clamp, STEEL),
                           ("yoke", yoke, YOKE_C), ("slider", slider, SLIDER_C)):
        asm.add(solid, name=f"{nm}_R", color=col)
        asm.add(mirror_L(solid), name=f"{nm}_L", color=col)

    asm.add(rod, name="yoke_rod_R", color=STEEL)           # bought shoulder screw (post + head top-stop)
    asm.add(mirror_L(rod), name="yoke_rod_L", color=STEEL)
    # Translucent worn-fit REFERENCE head — the head this pose is FITTED to (worn_head). The headphone
    # is already posed to land on it (band flexed, cups spread, earpads flush); the viewer swaps the
    # whole pose when a different size is picked. Shown OFF by default + held out of the explode.
    from parts.head_reference import make_head_reference
    head = make_head_reference(eh).translate((0, 0, P.head_ref_z))
    asm.add(head, name="head_ref", color=cq.Color(0.55, 0.70, 0.90, 0.28))
    asm.add(bow, name="bow_ref", color=STEEL)              # shared headband (REF)
    asm.add(pad, name="headband_pad", color=PAD_C)         # shared crown cushion
    asm.add(earpad, name="earpad_R", color=PAD_C)          # round pad mockup (bring your own)
    asm.add(mirror_L(earpad), name="earpad_L", color=PAD_C)
    asm.add(cover, name="headband_clamp_R", color=SLIDER_C)   # outer clamp plate
    asm.add(mirror_L(cover), name="headband_clamp_L", color=SLIDER_C)
    asm.add(shoe, name="slider_shoe_R", color=ORANGE)         # pressure pad (screw → shoe → post)
    asm.add(mirror_L(shoe), name="slider_shoe_L", color=ORANGE)

    # ACOUSTIC soft goods (VIZ) — shown so the damping + seal interfaces READ (internal; isolate or
    # explode to see them). Built in the cup/baffle local frame, posed with the cup (T_cup):
    #  • DAMPING felt disc — sits in the cup's damping ring, over the grille (⌀damping_felt × thickness).
    #  • Front-seal GASKET — a foam ring on the driver frame rim against the baffle seat (compressed).
    #  • Front acoustic PAPER/MESH — N "hot-dog" ARC STRIPS glued into the baffle's front depressions over
    #    the vent holes; the paper (not the hole size) sets the back→front resistance (bought soft-good,
    #    GRADE measurement-gated). Built in the baffle local frame (+baffle_seat_z), posed with the cup.
    FOAM = cq.Color(0.38, 0.52, 0.50)   # muted teal-grey: reads as acoustic foam/felt, distinct from pads
    PAPER = cq.Color(0.82, 0.76, 0.62)  # warm paper/mesh tan, distinct from the foam goods
    damping = T_cup(cq.Workplane("XY").workplane(offset=P.cup_interior_floor_z)
                    .circle(P.damping_felt_diameter / 2).extrude(P.damping_felt_thickness))
    gasket = T_cup(cq.Workplane("XY").workplane(offset=ledge_z - P.front_gasket_compressed)
                   .circle(P.driver_recess_diameter / 2)
                   .circle(P.driver_recess_diameter / 2 - P.front_gasket_width)
                   .extrude(P.front_gasket_compressed))

    def _paper_strip(zc):                # one acoustic-paper hot-dog (arc sector) in a baffle depression
        vin, vout = P.baffle_vent_inner_r - 0.5, P.baffle_vent_outer_r + 0.5
        ah = P.baffle_vent_strip_half + 2.0
        a0, a1, am = zc - ah, zc + ah, zc
        z0 = P.baffle_seat_z + P.baffle_ring_thickness - P.baffle_paper_recess_depth

        def pt(rad, deg):
            a = math.radians(deg)
            return (rad * math.cos(a), rad * math.sin(a))
        return (cq.Workplane("XY").workplane(offset=z0)
                .moveTo(*pt(vin, a0)).lineTo(*pt(vout, a0))
                .threePointArc(pt(vout, am), pt(vout, a1))
                .lineTo(*pt(vin, a1))
                .threePointArc(pt(vin, am), pt(vin, a0))
                .close().extrude(P.baffle_paper_thickness))
    paper = None
    for s in range(P.baffle_vent_strip_count):
        strip = _paper_strip(s * 360.0 / P.baffle_vent_strip_count)
        paper = strip if paper is None else paper.union(strip)
    paper = T_cup(paper)
    asm.add(damping, name="damping_R", color=FOAM)
    asm.add(mirror_L(damping), name="damping_L", color=FOAM)
    asm.add(gasket, name="gasket_R", color=FOAM)
    asm.add(mirror_L(gasket), name="gasket_L", color=FOAM)
    asm.add(paper, name="paper_R", color=PAPER)
    asm.add(mirror_L(paper), name="paper_L", color=PAPER)

    # Pivot hardware on both ears (viz), riding with the cup group. Guarded.
    try:
        from parts.hardware import make_shoulder_screw, make_heatset_insert
        for sign in (+1, -1):
            tag = "p" if sign > 0 else "m"
            insert = T_cup(make_heatset_insert()
                           .rotate((0, 0, 0), (0, 1, 0), -90 * sign)
                           .translate((sign * P.pivot_boss_outer_radius, 0, pbz)))
            screw = T_cup(make_shoulder_screw()
                          .rotate((0, 0, 0), (0, 1, 0), 90 * sign)
                          .translate((sign * (P.pivot_boss_outer_radius - P.yoke_arm_thickness),
                                      0, pbz)))
            asm.add(insert, name=f"insert_{tag}_R", color=BRASS)
            asm.add(mirror_L(insert), name=f"insert_{tag}_L", color=BRASS)
            asm.add(screw, name=f"screw_{tag}_R", color=SCREW_C)
            asm.add(mirror_L(screw), name=f"screw_{tag}_L", color=SCREW_C)
    except Exception as e:  # noqa: BLE001 — viz only; never block the build
        print(f"  [warn] assembly: pivot hardware skipped ({e}).")

    # Thumbscrew (8-32 knurled knob) — the height lock, shown so the post+SHOE+thumbscrew mechanism
    # reads. Rides with the slider on the +Y OUTBOARD boss; its tip presses the SHOE (not the
    # post). T_yoke maps local +Y → global +X, so in the worn pose the knurled head faces
    # straight out the side of the head (the natural two-finger reach with the phones ON).
    try:
        from parts.hardware import make_thumbscrew
        from parts.slider_shoe import shoe_offset_y
        shoe_face_y = shoe_offset_y() + P.slider_shoe_thickness / 2     # +Y face of the shoe
        ts = (make_thumbscrew()
              .rotate((0, 0, 0), (1, 0, 0), -90)                       # shaft → +Y (outboard), tip at origin
              .translate((0, shoe_face_y, 0))                          # tip on the shoe's +Y face
              .translate((0, 0, slider_z + P.slider_thumbscrew_boss_z)))  # ride with the boss
        ts_R = tilt(T_yoke(ts))
        asm.add(ts_R, name="thumbscrew_R", color=SCREW_C)
        asm.add(mirror_L(ts_R), name="thumbscrew_L", color=SCREW_C)
    except Exception as e:  # noqa: BLE001 — viz only
        print(f"  [warn] assembly: thumbscrew skipped ({e}).")
    return asm


if __name__ == "__main__":
    make_assembly().export("output/assembly.step")
    print("wrote output/assembly.step")

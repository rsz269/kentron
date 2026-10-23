# Generated with YoCAD — https://yotown.com · 2026-07-19
"""Compact Mars ice-drilling scout — a grounded chassis on four rolling wheels,
four deployable stabilizer legs that fold down to plant, a central vertical mast
that slides up out of its collar, a drill string that rotates while feeding down
into the sample zone and retracts, a conceptual sample container that transfers
along a deck rail into a rotating six-slot storage carousel, then the whole
system stows. Parts are built in their own local frames and placed by the
deterministic connector assembler (yc_assemble); every joint frame is derived
with yc_pt / yc_dir so slider edits keep the machine assembled. Conceptual motion
only — no ice-extraction performance is implied."""
import math
import cadquery as cq

# ── YoCAD assembler (injected — do not edit) ──
def _yc_norm(v):
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n < 1e-12:
        return [0.0, 0.0, 1.0]
    return [v[0] / n, v[1] / n, v[2] / n]


def _yc_cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _yc_frame(c):
    """Connector dict -> [R(3x3 rows, columns = x,y,z basis), t]."""
    z = _yc_norm(c["z"])
    x = _yc_norm(c["x"])
    y = _yc_norm(_yc_cross(z, x))
    x = _yc_cross(y, z)
    R = [[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]]
    o = c["origin"]
    return [R, [float(o[0]), float(o[1]), float(o[2])]]


def _yc_mul(A, B):
    Ra = A[0]
    Rb = B[0]
    R = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    i = 0
    while i < 3:
        j = 0
        while j < 3:
            R[i][j] = Ra[i][0] * Rb[0][j] + Ra[i][1] * Rb[1][j] + Ra[i][2] * Rb[2][j]
            j = j + 1
        i = i + 1
    tb = B[1]
    ta = A[1]
    t = [Ra[0][0] * tb[0] + Ra[0][1] * tb[1] + Ra[0][2] * tb[2] + ta[0],
         Ra[1][0] * tb[0] + Ra[1][1] * tb[1] + Ra[1][2] * tb[2] + ta[1],
         Ra[2][0] * tb[0] + Ra[2][1] * tb[1] + Ra[2][2] * tb[2] + ta[2]]
    return [R, t]


def _yc_inv(A):
    R = A[0]
    t = A[1]
    Rt = [[R[0][0], R[1][0], R[2][0]], [R[0][1], R[1][1], R[2][1]], [R[0][2], R[1][2], R[2][2]]]
    ti = [-(Rt[0][0] * t[0] + Rt[0][1] * t[1] + Rt[0][2] * t[2]),
          -(Rt[1][0] * t[0] + Rt[1][1] * t[1] + Rt[1][2] * t[2]),
          -(Rt[2][0] * t[0] + Rt[2][1] * t[1] + Rt[2][2] * t[2])]
    return [Rt, ti]


def _yc_axis_angle(R):
    tr = R[0][0] + R[1][1] + R[2][2]
    c = (tr - 1.0) / 2.0
    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0
    ang = math.acos(c)
    if ang < 1e-9:
        return [[0.0, 0.0, 1.0], 0.0]
    if ang > math.pi - 1e-6:
        xx = math.sqrt(max((R[0][0] + 1.0) / 2.0, 0.0))
        yy = math.sqrt(max((R[1][1] + 1.0) / 2.0, 0.0))
        zz = math.sqrt(max((R[2][2] + 1.0) / 2.0, 0.0))
        if xx >= yy and xx >= zz and xx > 1e-9:
            if R[0][1] < 0.0:
                yy = -yy
            if R[0][2] < 0.0:
                zz = -zz
            return [_yc_norm([xx, yy, zz]), 180.0]
        if yy >= zz and yy > 1e-9:
            if R[0][1] < 0.0:
                xx = -xx
            if R[1][2] < 0.0:
                zz = -zz
            return [_yc_norm([xx, yy, zz]), 180.0]
        if R[0][2] < 0.0:
            xx = -xx
        if R[1][2] < 0.0:
            yy = -yy
        return [_yc_norm([xx, yy, zz]), 180.0]
    s = 2.0 * math.sin(ang)
    ax = [(R[2][1] - R[1][2]) / s, (R[0][2] - R[2][0]) / s, (R[1][0] - R[0][1]) / s]
    return [_yc_norm(ax), ang * 180.0 / math.pi]


def yc_assemble(parts, connectors, mates, root):
    """Place parts by aligning mated connector frames. Returns [placed, xf]:
    placed = dict part -> transformed solid; xf = dict part -> world [R, t]."""
    flip = [[[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], [0.0, 0.0, 0.0]]
    ident = [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0]]
    xf = {}
    xf[root] = ident
    queue = [root]
    while len(queue) > 0:
        cur = queue.pop(0)
        k = 0
        while k < len(mates):
            a_id = mates[k][0]
            b_id = mates[k][1]
            ca = connectors[a_id]
            cb = connectors[b_id]
            pa = ca["part"]
            pb = cb["part"]
            if pa == cur and (pb in xf) == False:
                w = _yc_mul(_yc_mul(_yc_mul(xf[pa], _yc_frame(ca)), flip), _yc_inv(_yc_frame(cb)))
                xf[pb] = w
                queue.append(pb)
            if pb == cur and (pa in xf) == False:
                w = _yc_mul(_yc_mul(_yc_mul(xf[pb], _yc_frame(cb)), flip), _yc_inv(_yc_frame(ca)))
                xf[pa] = w
                queue.append(pa)
            k = k + 1
    placed = {}
    for name in parts:
        solid = parts[name]
        if name in xf:
            aa = _yc_axis_angle(xf[name][0])
            if aa[1] != 0.0:
                solid = solid.rotate((0, 0, 0), (aa[0][0], aa[0][1], aa[0][2]), aa[1])
            t = xf[name][1]
            solid = solid.translate((t[0], t[1], t[2]))
        placed[name] = solid
    return [placed, xf]


def yc_pt(xf, part, p):
    """Part-local point -> world (for JOINTS origins)."""
    A = xf[part]
    R = A[0]
    t = A[1]
    return [R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2] + t[0],
            R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2] + t[1],
            R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2] + t[2]]


def yc_dir(xf, part, v):
    """Part-local direction -> world (for JOINTS axes)."""
    R = xf[part][0]
    return _yc_norm([R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
                     R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
                     R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]])
# ── end YoCAD assembler ──


# ── Parameters (live sliders) ────────────────────────────────────────────────
BODY_L = 400.0        # chassis length [300:520]
BODY_W = 300.0        # chassis width [220:400]
BODY_H = 120.0        # chassis height [90:180]
WHEEL_R = 60.0        # wheel radius [40:90]
WHEEL_W = 70.0        # wheel width [40:100]
STAB_LEN = 150.0      # stabilizer leg length [120:200]
MAST_H = 300.0        # drill mast height [220:380]
MAST_RAISE = 180.0    # mast raise travel [80:220]
DRILL_LEN = 260.0     # drill string length [180:320]
DRILL_FEED = 220.0    # drill feed-down travel [120:260]
CAROUSEL_R = 80.0     # sample carousel radius [55:120]

# ── Fixed proportions / derived placement (tied to the params above) ─────────
WHEEL_INSET = 40.0                       # wheel hub inset from the body corner
STAB_INSET = 20.0                        # stabilizer pivot inset from the body edge
STAB_Z = BODY_H + 30.0                   # stabilizer pivot height, just above the deck
STAB_OUT = 65.0                          # outrigger reach beyond the wheels
MAST_X = 80.0                            # mast axis, forward of centre
MAST_BASE_Z = BODY_H + 40.0              # mast base sits just above the deck
CAROUSEL_X = -120.0                      # carousel axis, aft of centre
CAROUSEL_Z = 170.0                       # carousel hub height
TRANSFER_Z = 180.0                       # transfer-rail / container axis height
CONT_R = 20.0                            # sample container radius
CONT_LEN = 60.0                          # sample container body length
TRANSFER_TRAVEL = MAST_X - CAROUSEL_X    # drill-zone -> carousel = 200 mm
SHAFT_R = 15.0                           # drill string radius
SLOT_R = CAROUSEL_R * 0.62               # carousel slot pitch radius

# ── Colors (from the plan) ───────────────────────────────────────────────────
COL_CHASSIS = cq.Color(0.76, 0.71, 0.60)
COL_WHEEL = cq.Color(0.15, 0.15, 0.17)
COL_STAB = cq.Color(0.58, 0.60, 0.66)
COL_MAST = cq.Color(0.40, 0.46, 0.52)
COL_DRILL = cq.Color(0.88, 0.82, 0.20)
COL_CAROUSEL = cq.Color(0.20, 0.52, 0.82)
COL_CONTAINER = cq.Color(0.90, 0.32, 0.20)


# ── Local-frame part builders ────────────────────────────────────────────────
def make_chassis():
    body = cq.Workplane("XY").box(BODY_L, BODY_W, BODY_H).translate((0.0, 0.0, BODY_H / 2.0))
    try:
        body = body.edges("|Z").fillet(10.0)
    except Exception:
        pass
    # short mast collar the mast slides out of (top 2 mm below the mast base)
    collar_h = MAST_BASE_Z - 2.0 - BODY_H
    body = body.union(
        cq.Workplane("XY").box(100.0, 100.0, collar_h)
        .translate((MAST_X, 0.0, BODY_H + collar_h / 2.0)))
    # carousel support post rising to the hub
    body = body.union(
        cq.Workplane("XY").circle(34.0).extrude(CAROUSEL_Z - BODY_H)
        .translate((CAROUSEL_X, 0.0, BODY_H)))
    # transfer rail spanning drill zone -> carousel; its top carries the container
    rail_len = TRANSFER_TRAVEL + 40.0
    rail_top = TRANSFER_Z - CONT_R
    rail_h = rail_top - BODY_H
    body = body.union(
        cq.Workplane("XY").box(rail_len, 80.0, rail_h)
        .translate(((MAST_X + CAROUSEL_X) / 2.0, 0.0, BODY_H + rail_h / 2.0)))
    # stabilizer outriggers: each corner gets a vertical post from the deck, a
    # horizontal arm reaching out past the wheel, and a knuckle boss at the pivot
    # so every leg seats on solid chassis structure (no floating gap)
    stab_px = BODY_L / 2.0 - STAB_INSET
    stab_py = BODY_W / 2.0
    piv_y = BODY_W / 2.0 + STAB_OUT
    post_h = STAB_Z - BODY_H + 30.0
    arm_len = piv_y - stab_py + 24.0
    for _sx in (stab_px, -stab_px):
        for _sgn in (1.0, -1.0):
            _sy = _sgn * stab_py
            _py = _sgn * piv_y
            body = body.union(
                cq.Workplane("XY").box(56.0, 44.0, post_h)
                .translate((_sx, _sy, BODY_H + post_h / 2.0)))
            body = body.union(
                cq.Workplane("XY").box(48.0, arm_len, 40.0)
                .translate((_sx, (_sy + _py) / 2.0, STAB_Z)))
            body = body.union(
                cq.Workplane("XY").box(56.0, 52.0, 50.0)
                .translate((_sx, _py, STAB_Z)))
    return body


def make_wheel():
    # spin axis = local Z (the hub connector); disc lies in local XY
    wheel = (cq.Workplane("XY").circle(WHEEL_R).extrude(WHEEL_W)
             .translate((0.0, 0.0, -WHEEL_W / 2.0)))
    # tread groove for a readable rolling surface
    try:
        groove = (cq.Workplane("XY").circle(WHEEL_R + 1.0).circle(WHEEL_R - 6.0)
                  .extrude(WHEEL_W * 0.5).translate((0.0, 0.0, -WHEEL_W * 0.25)))
        wheel = wheel.cut(groove)
    except Exception:
        pass
    # rim hub for depth
    wheel = wheel.union(cq.Workplane("XY").circle(WHEEL_R * 0.35).extrude(WHEEL_W + 6.0)
                        .translate((0.0, 0.0, -(WHEEL_W + 6.0) / 2.0)))
    return wheel


def make_stabilizer(inward_sign):
    # leg hangs straight DOWN from the pivot at the local origin, so joint 0 is a
    # planted outrigger sitting on the ground; the deploy joint only lifts/stows it.
    # inward_sign kept for call symmetry (the leg is left/right symmetric).
    leg = cq.Workplane("XY").box(30.0, 30.0, STAB_LEN).translate((0.0, 0.0, -STAB_LEN / 2.0))
    foot = cq.Workplane("XY").box(72.0, 72.0, 22.0).translate((0.0, 0.0, -STAB_LEN))
    knuckle = cq.Workplane("XY").box(50.0, 60.0, 46.0)  # overlaps the chassis pivot boss
    leg = leg.union(foot).union(knuckle)
    try:
        leg = leg.edges("|X").fillet(4.0)
    except Exception:
        pass
    return leg


def make_mast():
    mast = cq.Workplane("XY").box(80.0, 80.0, MAST_H).translate((0.0, 0.0, MAST_H / 2.0))
    # spindle cap at the top where the drill string mounts
    mast = mast.union(cq.Workplane("XY").box(96.0, 96.0, 18.0)
                      .translate((0.0, 0.0, MAST_H - 9.0)))
    try:
        mast = mast.edges("|Z").fillet(6.0)
    except Exception:
        pass
    # drill guide collar at the mast foot: the drill string runs through this bore,
    # anchoring the rod to the mast instead of dangling free
    try:
        collar = (cq.Workplane("XY").circle(SHAFT_R + 13.0).circle(SHAFT_R + 3.0)
                  .extrude(26.0).translate((0.0, 0.0, -13.0)))
        mast = mast.union(collar)
    except Exception:
        pass
    return mast


def make_drill():
    # shaft along local Z, 0..DRILL_LEN; the top connector is at z = DRILL_LEN
    drill = cq.Workplane("XY").circle(SHAFT_R).extrude(DRILL_LEN)
    # conical bit tapering to a point at the bottom (z = 0)
    try:
        bit = cq.Workplane("XY").circle(SHAFT_R).extrude(-26.0, taper=58.0)
        drill = drill.union(bit)
    except Exception:
        drill = drill.union(cq.Workplane("XY").circle(SHAFT_R).extrude(-14.0))
    # drive collar near the top
    drill = drill.union(cq.Workplane("XY").circle(SHAFT_R + 6.0).extrude(14.0)
                        .translate((0.0, 0.0, DRILL_LEN - 14.0)))
    return drill


def make_carousel():
    disc = cq.Workplane("XY").circle(CAROUSEL_R).extrude(40.0).translate((0.0, 0.0, -20.0))
    # six conceptual sample slots around the rim
    slot_pts = [(SLOT_R * math.cos(2.0 * math.pi * i / 6.0),
                 SLOT_R * math.sin(2.0 * math.pi * i / 6.0)) for i in range(6)]
    try:
        disc = disc.faces(">Z").workplane().pushPoints(slot_pts).hole(2.0 * (CONT_R + 2.0), 30.0)
    except Exception:
        pass
    # central hub
    disc = disc.union(cq.Workplane("XY").circle(SLOT_R * 0.35).extrude(48.0)
                      .translate((0.0, 0.0, -24.0)))
    return disc


def make_container():
    # capsule along local Z (the dock connector axis); rounded ends
    body = (cq.Workplane("XY").circle(CONT_R).extrude(CONT_LEN)
            .translate((0.0, 0.0, -CONT_LEN / 2.0)))
    body = body.union(cq.Workplane("XY").sphere(CONT_R).translate((0.0, 0.0, CONT_LEN / 2.0)))
    body = body.union(cq.Workplane("XY").sphere(CONT_R).translate((0.0, 0.0, -CONT_LEN / 2.0)))
    return body


PARTS = {
    "chassis": make_chassis(),
    "wheel_fl": make_wheel(),
    "wheel_fr": make_wheel(),
    "wheel_rl": make_wheel(),
    "wheel_rr": make_wheel(),
    "stabilizer_fl": make_stabilizer(-1.0),
    "stabilizer_fr": make_stabilizer(1.0),
    "stabilizer_rl": make_stabilizer(-1.0),
    "stabilizer_rr": make_stabilizer(1.0),
    "drill_mast": make_mast(),
    "drill_string": make_drill(),
    "sample_carousel": make_carousel(),
    "sample_container": make_container(),
}

# ── Connector frames (origins tied to the named parameters) ──────────────────
WX = BODY_L / 2.0 - WHEEL_INSET   # 160 by default
WY = BODY_W / 2.0                 # 150 by default
SX = BODY_L / 2.0 - STAB_INSET    # 180 by default
PIV_Y = BODY_W / 2.0 + STAB_OUT   # stabilizer pivot, outboard of the wheels

CONNECTORS = {
    "c_chassis_wfl": {"part": "chassis", "origin": [WX, WY, WHEEL_R], "z": [0, 1, 0], "x": [1, 0, 0]},
    "c_wfl_hub":     {"part": "wheel_fl", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_chassis_wfr": {"part": "chassis", "origin": [WX, -WY, WHEEL_R], "z": [0, -1, 0], "x": [1, 0, 0]},
    "c_wfr_hub":     {"part": "wheel_fr", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_chassis_wrl": {"part": "chassis", "origin": [-WX, WY, WHEEL_R], "z": [0, 1, 0], "x": [1, 0, 0]},
    "c_wrl_hub":     {"part": "wheel_rl", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_chassis_wrr": {"part": "chassis", "origin": [-WX, -WY, WHEEL_R], "z": [0, -1, 0], "x": [1, 0, 0]},
    "c_wrr_hub":     {"part": "wheel_rr", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},

    "c_chassis_sfl": {"part": "chassis", "origin": [SX, PIV_Y, STAB_Z], "z": [1, 0, 0], "x": [0, 0, 1]},
    "c_sfl_pivot":   {"part": "stabilizer_fl", "origin": [0, 0, 0], "z": [1, 0, 0], "x": [0, 0, 1]},
    "c_chassis_sfr": {"part": "chassis", "origin": [SX, -PIV_Y, STAB_Z], "z": [1, 0, 0], "x": [0, 0, 1]},
    "c_sfr_pivot":   {"part": "stabilizer_fr", "origin": [0, 0, 0], "z": [1, 0, 0], "x": [0, 0, 1]},
    "c_chassis_srl": {"part": "chassis", "origin": [-SX, PIV_Y, STAB_Z], "z": [1, 0, 0], "x": [0, 0, 1]},
    "c_srl_pivot":   {"part": "stabilizer_rl", "origin": [0, 0, 0], "z": [1, 0, 0], "x": [0, 0, 1]},
    "c_chassis_srr": {"part": "chassis", "origin": [-SX, -PIV_Y, STAB_Z], "z": [1, 0, 0], "x": [0, 0, 1]},
    "c_srr_pivot":   {"part": "stabilizer_rr", "origin": [0, 0, 0], "z": [1, 0, 0], "x": [0, 0, 1]},

    "c_chassis_mast": {"part": "chassis", "origin": [MAST_X, 0, MAST_BASE_Z], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_mast_base":    {"part": "drill_mast", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_mast_spindle": {"part": "drill_mast", "origin": [0, 0, MAST_H], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_drill_top":    {"part": "drill_string", "origin": [0, 0, DRILL_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},

    "c_chassis_carousel": {"part": "chassis", "origin": [CAROUSEL_X, 0, CAROUSEL_Z], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_carousel_hub":     {"part": "sample_carousel", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},

    "c_chassis_transfer": {"part": "chassis", "origin": [MAST_X, 0, TRANSFER_Z], "z": [1, 0, 0], "x": [0, 0, 1]},
    "c_container_dock":   {"part": "sample_container", "origin": [0, 0, 0], "z": [1, 0, 0], "x": [0, 0, 1]},
}

MATES = [
    ["c_chassis_wfl", "c_wfl_hub"],
    ["c_chassis_wfr", "c_wfr_hub"],
    ["c_chassis_wrl", "c_wrl_hub"],
    ["c_chassis_wrr", "c_wrr_hub"],
    ["c_chassis_sfl", "c_sfl_pivot"],
    ["c_chassis_sfr", "c_sfr_pivot"],
    ["c_chassis_srl", "c_srl_pivot"],
    ["c_chassis_srr", "c_srr_pivot"],
    ["c_chassis_mast", "c_mast_base"],
    ["c_mast_spindle", "c_drill_top"],
    ["c_chassis_carousel", "c_carousel_hub"],
    ["c_chassis_transfer", "c_container_dock"],
]

# ── Deterministic placement ──────────────────────────────────────────────────
_asm = yc_assemble(PARTS, CONNECTORS, MATES, "chassis")
PLACED = _asm[0]
XF = _asm[1]

_COLORS = {
    "chassis": COL_CHASSIS,
    "wheel_fl": COL_WHEEL, "wheel_fr": COL_WHEEL, "wheel_rl": COL_WHEEL, "wheel_rr": COL_WHEEL,
    "stabilizer_fl": COL_STAB, "stabilizer_fr": COL_STAB,
    "stabilizer_rl": COL_STAB, "stabilizer_rr": COL_STAB,
    "drill_mast": COL_MAST, "drill_string": COL_DRILL,
    "sample_carousel": COL_CAROUSEL, "sample_container": COL_CONTAINER,
}

result = cq.Assembly()
_order = ["chassis", "wheel_fl", "wheel_fr", "wheel_rl", "wheel_rr",
          "stabilizer_fl", "stabilizer_fr", "stabilizer_rl", "stabilizer_rr",
          "drill_mast", "drill_string", "sample_carousel", "sample_container"]
for _n in _order:
    result = result.add(PLACED[_n], name=_n, color=_COLORS[_n])

# ── Joints (every origin/axis derived from the placed connector frames) ──────
def joint_origin(cid):
    c = CONNECTORS[cid]
    o = c["origin"]
    return yc_pt(XF, c["part"], (o[0], o[1], o[2]))


def joint_axis(part, axis):
    return yc_dir(XF, part, (axis[0], axis[1], axis[2]))


_SPECS = [
    {"name": "wheel_fl_spin", "type": "revolute", "child": "wheel_fl", "parent": "chassis",
     "at": "c_chassis_wfl", "axis": [0, 1, 0], "limits": [-3600.0, 3600.0]},
    {"name": "wheel_fr_spin", "type": "revolute", "child": "wheel_fr", "parent": "chassis",
     "at": "c_chassis_wfr", "axis": [0, 1, 0], "limits": [-3600.0, 3600.0]},
    {"name": "wheel_rl_spin", "type": "revolute", "child": "wheel_rl", "parent": "chassis",
     "at": "c_chassis_wrl", "axis": [0, 1, 0], "limits": [-3600.0, 3600.0]},
    {"name": "wheel_rr_spin", "type": "revolute", "child": "wheel_rr", "parent": "chassis",
     "at": "c_chassis_wrr", "axis": [0, 1, 0], "limits": [-3600.0, 3600.0]},

    {"name": "stabilizer_fl_deploy", "type": "revolute", "child": "stabilizer_fl", "parent": "chassis",
     "at": "c_chassis_sfl", "axis": [1, 0, 0], "limits": [0.0, 90.0]},
    {"name": "stabilizer_fr_deploy", "type": "revolute", "child": "stabilizer_fr", "parent": "chassis",
     "at": "c_chassis_sfr", "axis": [1, 0, 0], "limits": [-90.0, 0.0]},
    {"name": "stabilizer_rl_deploy", "type": "revolute", "child": "stabilizer_rl", "parent": "chassis",
     "at": "c_chassis_srl", "axis": [1, 0, 0], "limits": [0.0, 90.0]},
    {"name": "stabilizer_rr_deploy", "type": "revolute", "child": "stabilizer_rr", "parent": "chassis",
     "at": "c_chassis_srr", "axis": [1, 0, 0], "limits": [-90.0, 0.0]},

    {"name": "mast_raise", "type": "slider", "child": "drill_mast", "parent": "chassis",
     "at": "c_chassis_mast", "axis": [0, 0, 1], "limits": [0.0, MAST_RAISE]},

    {"name": "drill_feed_rotate", "type": "cylindrical", "child": "drill_string", "parent": "drill_mast",
     "at": "c_mast_spindle", "axis": [0, 0, 1], "limits": [-DRILL_FEED, 0.0]},

    {"name": "carousel_index", "type": "revolute", "child": "sample_carousel", "parent": "chassis",
     "at": "c_chassis_carousel", "axis": [0, 0, 1], "limits": [0.0, 360.0]},

    {"name": "sample_transfer", "type": "slider", "child": "sample_container", "parent": "chassis",
     "at": "c_chassis_transfer", "axis": [1, 0, 0], "limits": [-TRANSFER_TRAVEL, 0.0]},
]

JOINTS = []
for s in _SPECS:
    JOINTS.append({
        "name": s["name"], "type": s["type"], "parts": [s["child"], s["parent"]],
        "origin": joint_origin(s["at"]), "axis": joint_axis(s["parent"], s["axis"]),
        "limits": s["limits"],
    })

LINKS = []

# ── Animation ────────────────────────────────────────────────────────────────
# Sequence: roll in & stop -> deploy 4 stabilizers -> raise mast -> feed the
# rotating drill down and retract it -> transfer a sample container into the
# carousel and index it -> stow everything -> drive off. Every excursion is
# anchored back to rest so nothing drifts across the gaps.
_wheel_keys = [[0.0, 0.0], [2.0, -720.0], [15.5, -720.0], [18.0, -1080.0]]
ANIMATION = {"fps": 30, "duration": 18.0, "tracks": [
    {"joint": "wheel_fl_spin", "unit": "deg", "keys": _wheel_keys},
    {"joint": "wheel_fr_spin", "unit": "deg", "keys": _wheel_keys},
    {"joint": "wheel_rl_spin", "unit": "deg", "keys": _wheel_keys},
    {"joint": "wheel_rr_spin", "unit": "deg", "keys": _wheel_keys},

    {"joint": "stabilizer_fl_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [2.5, 0.0], [3.3, 70.0], [4.2, 0.0], [15.5, 0.0], [16.3, 70.0], [17.1, 0.0], [18.0, 0.0]]},
    {"joint": "stabilizer_rl_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [2.5, 0.0], [3.3, 70.0], [4.2, 0.0], [15.5, 0.0], [16.3, 70.0], [17.1, 0.0], [18.0, 0.0]]},
    {"joint": "stabilizer_fr_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [2.5, 0.0], [3.3, -70.0], [4.2, 0.0], [15.5, 0.0], [16.3, -70.0], [17.1, 0.0], [18.0, 0.0]]},
    {"joint": "stabilizer_rr_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [2.5, 0.0], [3.3, -70.0], [4.2, 0.0], [15.5, 0.0], [16.3, -70.0], [17.1, 0.0], [18.0, 0.0]]},

    {"joint": "mast_raise", "unit": "mm",
     "keys": [[0.0, 0.0], [4.0, 0.0], [6.0, MAST_RAISE], [14.5, MAST_RAISE],
              [16.5, 0.0], [18.0, 0.0]]},

    # cylindrical drill: slide feeds down then retracts; rotate spins throughout
    {"joint": "drill_feed_rotate", "unit": "mm", "dof": "slide",
     "keys": [[0.0, 0.0], [6.5, 0.0], [9.0, -DRILL_FEED], [10.5, -DRILL_FEED],
              [12.5, 0.0], [18.0, 0.0]]},
    {"joint": "drill_feed_rotate", "unit": "deg", "dof": "rotate",
     "keys": [[0.0, 0.0], [6.5, 0.0], [12.5, -2880.0], [18.0, -2880.0]]},

    {"joint": "sample_transfer", "unit": "mm",
     "keys": [[0.0, 0.0], [13.0, 0.0], [14.5, -TRANSFER_TRAVEL], [15.5, -TRANSFER_TRAVEL],
              [16.5, 0.0], [18.0, 0.0]]},
    {"joint": "carousel_index", "unit": "deg",
     "keys": [[0.0, 0.0], [14.5, 0.0], [15.5, 120.0], [17.0, 120.0], [18.0, 0.0]]},
]}
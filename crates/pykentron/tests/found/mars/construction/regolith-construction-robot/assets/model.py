# Generated with YoCAD — https://yotown.com · 2026-07-19 · run 019f7b04-d214-7c7d-aa2d-b51adea9a122
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



# ── Parameters (live sliders) ─────────────────────────────────────────────────
BASE_H = 420.0          # chassis height [300:520]
UPPER_H = 450.0         # turret height [300:600]
BOOM_LEN = 900.0        # boom pivot-to-pivot length [600:1100]
STICK_LEN = 650.0       # stick pivot-to-pivot length [450:850]
BUCKET_W = 500.0        # bucket width [350:650]
HOPPER_H = 500.0        # hopper height [350:650]
STAB_LEN = 380.0        # stabilizer leg length [250:500]
STAB_DEPLOY = 80.0      # stabilizer deploy angle, deg [40:90]
TURRET_SWING_DEG = 150.0  # turret slew toward hopper, deg [60:180]

# ── Fixed proportions derived from the plan ───────────────────────────────────
BASE_L = 1300.0
BASE_W = 900.0
TRACK_H = 200.0
TRACK_W = 220.0
TRACK_Y = BASE_W / 2.0 - TRACK_W / 2.0
DECK_L = 1180.0

UPPER_L = 800.0
UPPER_W = 700.0
# Boom pivot on the FRONT-TOP shoulder of the turret house, well above the deck
# so the arm's arc clears the tracks and the rear hopper at both extremes.
BOOM_MOUNT_X = UPPER_L * 0.42      # ~336 → front shoulder of the turret
BOOM_MOUNT_Z = UPPER_H * 0.72      # ~324 → high shoulder pivot

BOOM_RISE = 150.0                  # boom tip rises 150 mm over its length
BOOM_W = 170.0
BOOM_TH = 230.0

STICK_DROP = -60.0                 # stick tip drops 60 mm
STICK_W = 140.0
STICK_TH = 190.0

BUCKET_L = 350.0
BUCKET_D = 320.0

HOPPER_L = 600.0
HOPPER_W = 700.0
HOPPER_MOUNT_X = -640.0            # pushed back to clear the turret slew circle
HOPPER_MOUNT_Z = BASE_H - 100.0    # rear-top hinge, follows chassis height

# Stabilizer hinge geometry — ONE set of numbers shared by the chassis bracket,
# the leg boss, and the joint origin so they physically coincide.
STAB_X = 520.0
STAB_Y = 420.0
STAB_Z = 90.0
STAB_BOSS_R = 74.0                 # hinge-pin boss radius (leg + bracket share it)
BRK_W = 210.0                      # chassis hinge bracket size
BRK_D = 180.0
BRK_H = 170.0

# Shared arm pivot-boss radii so mated bosses always interpenetrate (≤2 mm clr)
STICK_END_R = STICK_TH * 0.58      # stick tip boss
BUCKET_HUB_R = STICK_TH * 0.58     # bucket hub — same radius, concentric on axis

# ── Colors (from the plan) ────────────────────────────────────────────────────
COL_BASE = cq.Color(0.28, 0.29, 0.32)
COL_UPPER = cq.Color(0.82, 0.45, 0.20)
COL_ARM = cq.Color(0.86, 0.60, 0.28)
COL_BUCKET = cq.Color(0.55, 0.56, 0.60)
COL_HOPPER = cq.Color(0.75, 0.38, 0.16)
COL_STAB = cq.Color(0.35, 0.36, 0.40)


# ── Helpers ───────────────────────────────────────────────────────────────────
def cyl_y(r, length, center):
    """Cylinder whose axis runs along Y, centered at `center`."""
    c = (cq.Workplane("XY").circle(r).extrude(length)
         .translate((0.0, 0.0, -length / 2.0))
         .rotate((0, 0, 0), (1, 0, 0), 90))
    return c.translate((center[0], center[1], center[2]))


def bar_xz(p0, p1, width, height):
    """Rectangular bar (width along Y, height along Z) joining two points that
    lie in the XZ plane — used for the boom and stick segments."""
    dx = p1[0] - p0[0]
    dz = p1[2] - p0[2]
    length = math.hypot(dx, dz)
    pitch = math.degrees(math.atan2(dz, dx))
    b = (cq.Workplane("XY").box(length, width, height)
         .translate((length / 2.0, 0.0, 0.0))
         .rotate((0, 0, 0), (0, 1, 0), -pitch))
    return b.translate((p0[0], p0[1], p0[2]))


def try_fillet(wp, selector, radius):
    try:
        return wp.edges(selector).fillet(radius)
    except Exception:
        return wp


# ── track_base: tracked chassis (grounded, permanent static root) ─────────────
def make_track_base():
    # two tracked bogies straddling the deck, resting on the ground (z = 0)
    part = None
    for sy in (1.0, -1.0):
        trk = (cq.Workplane("XY").box(BASE_L, TRACK_W, TRACK_H)
               .translate((0.0, sy * TRACK_Y, TRACK_H / 2.0)))
        trk = try_fillet(trk, "|Y", 60.0)
        # a few road-wheel bosses for read
        for wx in (-420.0, -140.0, 140.0, 420.0):
            trk = trk.union(cyl_y(TRACK_H * 0.42, TRACK_W + 8.0,
                                  (wx, sy * TRACK_Y, TRACK_H * 0.42)))
        part = trk if part is None else part.union(trk)

    # main deck spanning the bogies, up to the chassis top
    deck = (cq.Workplane("XY").box(DECK_L, BASE_W - 30.0, BASE_H - TRACK_H)
            .translate((0.0, 0.0, TRACK_H + (BASE_H - TRACK_H) / 2.0)))
    deck = try_fillet(deck, "|Z", 40.0)
    part = part.union(deck)

    # slew ring the turret sits on (centered on the vertical swing axis)
    part = part.union(cq.Workplane("XY").circle(255.0).extrude(28.0)
                      .translate((0.0, 0.0, BASE_H - 28.0)))

    # rear mount pad for the hopper hinge
    part = part.union(cq.Workplane("XY").box(120.0, HOPPER_W + 40.0, 120.0)
                      .translate((HOPPER_MOUNT_X, 0.0, HOPPER_MOUNT_Z)))

    # four stabilizer hinge brackets on the chassis corners — each straddles the
    # pivot point (STAB_X, STAB_Y, STAB_Z) and carries a boss on the hinge axis
    # so the leg boss below plunges into it with only a hair of clearance.
    for sx, sy in ((STAB_X, STAB_Y), (STAB_X, -STAB_Y),
                   (-STAB_X, STAB_Y), (-STAB_X, -STAB_Y)):
        part = part.union(cq.Workplane("XY").box(BRK_W, BRK_D, BRK_H)
                          .translate((sx, sy, STAB_Z)))
        part = part.union(cyl_y(STAB_BOSS_R + 6.0, BRK_D + 30.0, (sx, sy, STAB_Z)))
    return part


# ── upper_body: slewing turret house — flush counterweight + cab, no floaters ─
def make_upper_body():
    body = (cq.Workplane("XY").box(UPPER_L, UPPER_W, UPPER_H)
            .translate((0.0, 0.0, UPPER_H / 2.0)))
    body = try_fillet(body, "|Z", 50.0)

    # rear counterweight: a distinct block that PROTRUDES off the back of the
    # house (overlapping it by ~50 mm for a solid bond) so it reads as a
    # counterweight rather than an interior box.
    cw = (cq.Workplane("XY").box(180.0, UPPER_W + 30.0, UPPER_H * 0.78)
          .translate((-UPPER_L / 2.0 - 40.0, 0.0, UPPER_H * 0.39)))
    cw = try_fillet(cw, "|Z", 24.0)
    body = body.union(cw)

    # operator cab: a smaller box seated flush on the front-top of the house,
    # embedded 30 mm into the roof (bottom below UPPER_H) so it is unambiguously
    # attached and never floats or interpenetrates the counterweight.
    cab = (cq.Workplane("XY").box(300.0, 250.0, 300.0)
           .translate((150.0, UPPER_W / 2.0 - 150.0, UPPER_H + 120.0)))
    cab = try_fillet(cab, "|Z", 28.0)
    body = body.union(cab)

    # boom pivot boss at the plan's front-shoulder mount point, fed by a stout
    # column rising from the deck so the pivot is fully supported.
    body = body.union(cyl_y(95.0, 260.0, (BOOM_MOUNT_X, 0.0, BOOM_MOUNT_Z)))
    body = body.union(cq.Workplane("XY").box(210.0, 300.0, BOOM_MOUNT_Z)
                      .translate((BOOM_MOUNT_X, 0.0, BOOM_MOUNT_Z / 2.0)))
    return body


# ── boom: first arm segment ──────────────────────────────────────────────────
def make_boom():
    b = bar_xz((0.0, 0.0, 0.0), (BOOM_LEN, 0.0, BOOM_RISE), BOOM_W, BOOM_TH)
    # concept "elbow" — a raised mid-hump so the profile isn't a plain bar
    b = b.union(bar_xz((BOOM_LEN * 0.5, 0.0, BOOM_RISE * 0.5 + 60.0),
                       (BOOM_LEN, 0.0, BOOM_RISE), BOOM_W, BOOM_TH * 0.7))
    b = b.union(cyl_y(BOOM_TH * 0.55, BOOM_W + 40.0, (0.0, 0.0, 0.0)))
    b = b.union(cyl_y(BOOM_TH * 0.45, BOOM_W + 40.0, (BOOM_LEN, 0.0, BOOM_RISE)))
    return b


# ── stick: second arm segment ────────────────────────────────────────────────
def make_stick():
    s = bar_xz((0.0, 0.0, 0.0), (STICK_LEN, 0.0, STICK_DROP), STICK_W, STICK_TH)
    s = s.union(cyl_y(STICK_TH * 0.55, STICK_W + 36.0, (0.0, 0.0, 0.0)))
    # tip boss on the bucket hinge axis — large enough that the bucket hub
    # engulfs it through the whole curl (no gap can open).
    s = s.union(cyl_y(STICK_END_R, STICK_W + 56.0, (STICK_LEN, 0.0, STICK_DROP)))
    return s


# ── bucket: scoop hinged at its top-back edge ────────────────────────────────
def make_bucket():
    cx = BUCKET_L / 2.0 - 40.0
    outer = (cq.Workplane("XY").box(BUCKET_L, BUCKET_W, BUCKET_D)
             .translate((cx, 0.0, -BUCKET_D / 2.0)))
    cavity = (cq.Workplane("XY").box(BUCKET_L - 80.0, BUCKET_W - 80.0, BUCKET_D)
              .translate((cx + 20.0, 0.0, -BUCKET_D / 2.0 + 60.0)))
    bkt = outer.cut(cavity)
    bkt = try_fillet(bkt, "|Y", 40.0)
    # hub on the pivot axis, concentric with (and overlapping) the stick tip
    # boss so the two solids stay physically joined through the full curl.
    bkt = bkt.union(cyl_y(BUCKET_HUB_R, BUCKET_W + 30.0, (0.0, 0.0, 0.0)))
    return bkt


# ── hopper: rear processing bin that tilts to dump ───────────────────────────
def make_hopper():
    outer = (cq.Workplane("XY").box(HOPPER_L, HOPPER_W, HOPPER_H)
             .translate((-HOPPER_L / 2.0, 0.0, HOPPER_H / 2.0)))
    outer = try_fillet(outer, "|Z", 40.0)
    cavity = (cq.Workplane("XY").box(HOPPER_L - 80.0, HOPPER_W - 80.0, HOPPER_H)
              .translate((-HOPPER_L / 2.0, 0.0, HOPPER_H / 2.0 + 40.0)))
    hop = outer.cut(cavity)
    # hinge boss at the front-bottom edge (on the dump axis)
    hop = hop.union(cyl_y(60.0, HOPPER_W + 20.0, (0.0, 0.0, 0.0)))
    return hop


# ── stabilizer: outrigger leg + foot pad (one geometry, two instances) ───────
def make_stabilizer():
    leg = (cq.Workplane("XY").box(STAB_LEN, 120.0, 100.0)
           .translate((STAB_LEN / 2.0, 0.0, 0.0)))
    leg = try_fillet(leg, "|Y", 24.0)
    # hinge boss sits on the pivot origin and plunges into the chassis bracket
    # boss (same shared axis point) so the leg is bolted to the chassis at rest
    # and through the whole deploy swing.
    leg = leg.union(cyl_y(STAB_BOSS_R, 150.0, (0.0, 0.0, 0.0)))
    foot = (cq.Workplane("XY").box(150.0, 200.0, 120.0)
            .translate((STAB_LEN, 0.0, -40.0)))                    # planting pad
    foot = try_fillet(foot, "|Z", 20.0)
    return leg.union(foot)


PARTS = {
    "track_base": make_track_base(),
    "upper_body": make_upper_body(),
    "boom": make_boom(),
    "stick": make_stick(),
    "bucket": make_bucket(),
    "hopper": make_hopper(),
    "stabilizer_front": make_stabilizer(),
    "stabilizer_rear": make_stabilizer(),
}

# ── Connector frames (origins tied to parameters). The base↔upper mate is a
#    pure Z-stack: base connector Z=+up, upper connector Z=-up, so the turret
#    seats level on the turntable and the swing is a flat yaw about vertical. ─
CONNECTORS = {
    "c_base_swing":  {"part": "track_base", "origin": [0.0, 0.0, BASE_H],
                      "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_upper_swing": {"part": "upper_body", "origin": [0.0, 0.0, 0.0],
                      "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_upper_boom":  {"part": "upper_body", "origin": [BOOM_MOUNT_X, 0.0, BOOM_MOUNT_Z],
                      "z": [0, 1, 0], "x": [1, 0, 0]},
    "c_boom_root":   {"part": "boom", "origin": [0.0, 0.0, 0.0],
                      "z": [0, -1, 0], "x": [1, 0, 0]},
    "c_boom_stick":  {"part": "boom", "origin": [BOOM_LEN, 0.0, BOOM_RISE],
                      "z": [0, 1, 0], "x": [1, 0, 0]},
    "c_stick_root":  {"part": "stick", "origin": [0.0, 0.0, 0.0],
                      "z": [0, -1, 0], "x": [1, 0, 0]},
    "c_stick_bucket": {"part": "stick", "origin": [STICK_LEN, 0.0, STICK_DROP],
                       "z": [0, 1, 0], "x": [1, 0, 0]},
    "c_bucket_root": {"part": "bucket", "origin": [0.0, 0.0, 0.0],
                      "z": [0, -1, 0], "x": [1, 0, 0]},
    "c_base_hopper": {"part": "track_base", "origin": [HOPPER_MOUNT_X, 0.0, HOPPER_MOUNT_Z],
                      "z": [0, 1, 0], "x": [1, 0, 0]},
    "c_hopper_root": {"part": "hopper", "origin": [0.0, 0.0, 0.0],
                      "z": [0, -1, 0], "x": [1, 0, 0]},
    "c_base_stab_front": {"part": "track_base", "origin": [STAB_X, STAB_Y, STAB_Z],
                          "z": [0, 1, 0], "x": [1, 0, 0]},
    "c_stab_front_root": {"part": "stabilizer_front", "origin": [0.0, 0.0, 0.0],
                          "z": [0, -1, 0], "x": [1, 0, 0]},
    "c_base_stab_rear": {"part": "track_base", "origin": [-STAB_X, -STAB_Y, STAB_Z],
                         "z": [0, 1, 0], "x": [1, 0, 0]},
    "c_stab_rear_root": {"part": "stabilizer_rear", "origin": [0.0, 0.0, 0.0],
                         "z": [0, -1, 0], "x": [1, 0, 0]},
}

MATES = [
    ["c_base_swing", "c_upper_swing"],
    ["c_upper_boom", "c_boom_root"],
    ["c_boom_stick", "c_stick_root"],
    ["c_stick_bucket", "c_bucket_root"],
    ["c_base_hopper", "c_hopper_root"],
    ["c_base_stab_front", "c_stab_front_root"],
    ["c_base_stab_rear", "c_stab_rear_root"],
]

_asm = yc_assemble(PARTS, CONNECTORS, MATES, "track_base")
PLACED = _asm[0]
XF = _asm[1]

# ── Assembly ──────────────────────────────────────────────────────────────────
_COLORS = {
    "track_base": COL_BASE, "upper_body": COL_UPPER, "boom": COL_ARM,
    "stick": COL_ARM, "bucket": COL_BUCKET, "hopper": COL_HOPPER,
    "stabilizer_front": COL_STAB, "stabilizer_rear": COL_STAB,
}
_ORDER = ["track_base", "upper_body", "boom", "stick", "bucket", "hopper",
          "stabilizer_front", "stabilizer_rear"]

result = cq.Assembly()
for _n in _ORDER:
    result = result.add(PLACED[_n], name=_n, color=_COLORS[_n])

# ── Joints (origins/axes derived from the placed frames). track_base is the
#    PARENT of every joint, so it is fixed to ground and can never tilt; only
#    the cab, arm, hopper and legs move relative to it. ────────────────────────
JOINTS = [
    # Vertical yaw: axis = base +Z (up), origin centered on the turntable.
    {"name": "turret_swing", "type": "revolute",
     "parts": ["upper_body", "track_base"],
     "origin": yc_pt(XF, "track_base", [0.0, 0.0, BASE_H]),
     "axis": yc_dir(XF, "track_base", [0, 0, 1]), "limits": [-180.0, 180.0]},
    {"name": "boom_lift", "type": "revolute",
     "parts": ["boom", "upper_body"],
     "origin": yc_pt(XF, "upper_body", [BOOM_MOUNT_X, 0.0, BOOM_MOUNT_Z]),
     "axis": yc_dir(XF, "upper_body", [0, -1, 0]), "limits": [-20.0, 60.0]},
    {"name": "stick_curl", "type": "revolute",
     "parts": ["stick", "boom"],
     "origin": yc_pt(XF, "boom", [BOOM_LEN, 0.0, BOOM_RISE]),
     "axis": yc_dir(XF, "boom", [0, 1, 0]), "limits": [-100.0, 20.0]},
    {"name": "bucket_curl", "type": "revolute",
     "parts": ["bucket", "stick"],
     "origin": yc_pt(XF, "stick", [STICK_LEN, 0.0, STICK_DROP]),
     "axis": yc_dir(XF, "stick", [0, 1, 0]), "limits": [-140.0, 30.0]},
    {"name": "hopper_dump", "type": "revolute",
     "parts": ["hopper", "track_base"],
     "origin": yc_pt(XF, "track_base", [HOPPER_MOUNT_X, 0.0, HOPPER_MOUNT_Z]),
     "axis": yc_dir(XF, "track_base", [0, 1, 0]), "limits": [0.0, 55.0]},
    {"name": "stabilizer_front_deploy", "type": "revolute",
     "parts": ["stabilizer_front", "track_base"],
     "origin": yc_pt(XF, "track_base", [STAB_X, STAB_Y, STAB_Z]),
     "axis": yc_dir(XF, "track_base", [0, 1, 0]), "limits": [0.0, 80.0]},
    {"name": "stabilizer_rear_deploy", "type": "revolute",
     "parts": ["stabilizer_rear", "track_base"],
     "origin": yc_pt(XF, "track_base", [-STAB_X, -STAB_Y, STAB_Z]),
     "axis": yc_dir(XF, "track_base", [0, 1, 0]), "limits": [0.0, 80.0]},
]
LINKS = []

# ── Animation: arrive → deploy → dig → lift → slew → dump → reset. The turret
#    slew (t8→10.5) happens ONLY while the arm is folded high (boom +55), so the
#    bucket never sweeps into the tracks or the rear hopper. ───────────────────
_SW = TURRET_SWING_DEG
_SD = STAB_DEPLOY
ANIMATION = {"fps": 30, "duration": 17.0, "tracks": [
    # outriggers stow, then deploy after arrival, retract at the end
    {"joint": "stabilizer_front_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [1.5, 0.0], [3.0, _SD], [14.0, _SD], [15.5, 0.0], [17.0, 0.0]]},
    {"joint": "stabilizer_rear_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [1.5, 0.0], [3.0, _SD], [14.0, _SD], [15.5, 0.0], [17.0, 0.0]]},
    # boom: travel raised, drop to dig, lift high for the carry/slew, settle
    {"joint": "boom_lift", "unit": "deg",
     "keys": [[0.0, 25.0], [1.5, 25.0], [3.0, 20.0], [6.0, -15.0], [7.0, -15.0],
              [8.0, 55.0], [13.0, 55.0], [15.5, 25.0], [17.0, 25.0]]},
    # stick: extend to reach ground, curl in to hold the load, reset
    {"joint": "stick_curl", "unit": "deg",
     "keys": [[0.0, -30.0], [1.5, -30.0], [3.0, -20.0], [6.0, 15.0], [8.0, -70.0],
              [13.0, -70.0], [15.5, -30.0], [17.0, -30.0]]},
    # bucket: curl closed to scoop, hold, then open over the hopper to dump
    {"joint": "bucket_curl", "unit": "deg",
     "keys": [[0.0, -20.0], [1.5, -20.0], [3.0, -45.0], [6.0, -120.0], [8.0, -120.0],
              [10.5, -120.0], [12.0, 15.0], [12.5, 15.0], [15.0, -20.0], [17.0, -20.0]]},
    # turret: flat yaw toward the rear hopper for the dump, then yaw back
    {"joint": "turret_swing", "unit": "deg",
     "keys": [[0.0, 0.0], [8.0, 0.0], [10.5, _SW], [13.0, _SW], [15.5, 0.0], [17.0, 0.0]]},
    # hopper: tilt open to release material as the bucket dumps, then close
    {"joint": "hopper_dump", "unit": "deg",
     "keys": [[0.0, 0.0], [11.5, 0.0], [12.5, 50.0], [13.0, 50.0], [15.0, 0.0], [17.0, 0.0]]},
]}

# Generated with YoCAD — https://yotown.com · 2026-07-19 · run 019f792a-a0c4-78c6-a87e-0132395b4342
"""Mars habitat delivery pod — an original deployable surface concept (not a
pressure- or human-rated design). A compact faceted octagonal shell rides on
four articulated landing legs with independent vertical leveling-slider feet.
Two side shells unfold down into floor platforms, a telescoping roof section
raises off the body top, and a front access ramp lowers to the surface — the
animation plays the whole landing-and-deploy sequence, then stows back to rest.

Parts are built in their own local frames and placed by the runtime connector
assembler (yc_assemble); every joint origin/axis is derived with yc_pt/yc_dir
so slider edits keep the pod attached to itself."""
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
BODY_ACROSS = 2400.0      # habitat shell flat-to-flat width [1800:3200]
BODY_HEIGHT = 1800.0      # habitat shell height [1200:2400]
LEG_LENGTH = 1400.0       # landing leg length [1000:1900]
LEG_DEPLOY_DEG = 72.0     # leg outward deploy angle [30:85]
LEVEL_TRAVEL = 300.0      # leveling foot vertical travel [100:400]
SHELL_UNFOLD_DEG = 90.0   # side-shell unfold angle [45:100]
ROOF_RAISE = 500.0        # telescoping roof rise [200:700]
RAMP_LOWER_DEG = 70.0     # access ramp lowering angle [30:85]

# ── Derived proportions (kept in sync with the parameters above) ─────────────
HALF = BODY_ACROSS / 2.0            # body half-width (side / ramp faces)
HZ = BODY_HEIGHT / 2.0             # body half-height (top face = roof seat)
LEG_OFF = (HALF - 200.0) / math.sqrt(2.0)  # leg mount on the octagon diagonal, kept 200 mm inside the corner vertex
LEG_Z = 200.0                    # hip height on the body
LEG_HALF = LEG_LENGTH / 2.0        # hip at +LEG_HALF, foot at -LEG_HALF (local Z)
LEG_R = 60.0                     # leg tube radius (120 mm dia)
FOOT_CONN = 150.0                # slider connector height inside the foot shaft
FOOT_SHAFT = 210.0               # foot shaft length
FOOT_PAD_R = 100.0               # foot pad radius (200 mm dia)
SHELL_W = 2200.0                 # side-shell width (local X)
SHELL_H = 1400.0                 # side-shell height/reach (local Z)
SHELL_T = 80.0                   # side-shell thickness (local Y)
ROOF_ACROSS = 2000.0             # roof cap flat-to-flat width
ROOF_H = 400.0                   # roof cap height
RAMP_W = 1200.0                  # ramp width (local X)
RAMP_DEPTH = 900.0               # ramp length (local Y)
RAMP_T = 60.0                    # ramp thickness (local Z)

# ── Colors (from the plan) ───────────────────────────────────────────────────
COL_BODY = cq.Color(0.82, 0.40, 0.24)
COL_LEG = cq.Color(0.35, 0.37, 0.40)
COL_FOOT = cq.Color(0.55, 0.57, 0.60)
COL_SHELL = cq.Color(0.86, 0.60, 0.30)
COL_ROOF = cq.Color(0.70, 0.72, 0.75)
COL_RAMP = cq.Color(0.50, 0.45, 0.42)


# ── Part builders — each in its OWN local frame, Z up ────────────────────────
def make_body():
    """Faceted octagonal shell centered at the origin, top-open roof pocket,
    front door opening for the ramp, and a mid belt ring."""
    b = (cq.Workplane("XY").polygon(8, BODY_ACROSS).extrude(BODY_HEIGHT)
         .translate((0.0, 0.0, -HZ)))
    # roof pocket in the top so the telescoping cap nests when stowed
    pocket = (cq.Workplane("XY").polygon(8, ROOF_ACROSS - 60.0)
              .extrude(ROOF_H + 10.0)
              .translate((0.0, 0.0, HZ - ROOF_H / 2.0 + 5.0)))
    b = b.cut(pocket)
    # front door opening (covered by the ramp) on the +Y face
    door = (cq.Workplane("XY").box(RAMP_W, 400.0, 1000.0)
            .translate((0.0, HALF, 180.0)))
    b = b.cut(door)
    # decorative belt ring at mid-height
    belt_out = cq.Workplane("XY").polygon(8, BODY_ACROSS + 24.0).extrude(70.0)
    belt_in = (cq.Workplane("XY").polygon(8, BODY_ACROSS - 24.0).extrude(90.0)
               .translate((0.0, 0.0, -10.0)))
    belt = belt_out.cut(belt_in).translate((0.0, 0.0, -35.0))
    b = b.union(belt)
    try:
        b = b.edges(">Z").chamfer(28.0)
    except Exception:
        pass
    return b


def make_leg():
    """Vertical leg tube along local Z: hip clevis at +LEG_HALF, foot cap at
    -LEG_HALF; both connector points lie on the solid."""
    shaft = (cq.Workplane("XY").circle(LEG_R).extrude(LEG_LENGTH)
             .translate((0.0, 0.0, -LEG_HALF)))
    hip = (cq.Workplane("XY").box(LEG_R * 2.6, LEG_R * 1.7, LEG_R * 1.8)
           .translate((0.0, 0.0, LEG_HALF - 50.0)))
    knee = cq.Workplane("XY").sphere(LEG_R * 1.05)
    foot_cap = (cq.Workplane("XY").circle(LEG_R * 1.35).extrude(46.0)
                .translate((0.0, 0.0, -LEG_HALF)))
    return shaft.union(hip).union(knee).union(foot_cap)


def make_foot():
    """Leveling foot: shaft up the local +Z (slider rides in it) and a wide
    ground pad below; slider connector at z = FOOT_CONN sits on the shaft."""
    shaft = cq.Workplane("XY").circle(46.0).extrude(FOOT_SHAFT)
    pad = (cq.Workplane("XY").circle(FOOT_PAD_R).extrude(52.0)
           .translate((0.0, 0.0, -52.0)))
    rim = (cq.Workplane("XY").circle(FOOT_PAD_R * 0.6).extrude(24.0)
           .translate((0.0, 0.0, -52.0)))
    return shaft.union(pad).union(rim)


def make_shell():
    """Unfolding side shell / floor platform: width along local X, reach along
    local Z (hinge edge at +SHELL_H/2), thin in local Y, with a floor rail."""
    panel = cq.Workplane("XY").box(SHELL_W, SHELL_T, SHELL_H)
    rail = (cq.Workplane("XY").box(SHELL_W, SHELL_T * 2.4, 90.0)
            .translate((0.0, 0.0, -SHELL_H / 2.0 + 45.0)))
    ribs = None
    for xr in [-SHELL_W / 2.0 + 300.0, 0.0, SHELL_W / 2.0 - 300.0]:
        rib = (cq.Workplane("XY").box(70.0, SHELL_T * 1.8, SHELL_H - 160.0)
               .translate((xr, 0.0, 0.0)))
        ribs = rib if ribs is None else ribs.union(rib)
    return panel.union(rail).union(ribs)


def make_roof():
    """Telescoping roof cap centered at the origin (bottom at -ROOF_H/2 = the
    connector seat), with a stepped cap and a slim mast."""
    r = (cq.Workplane("XY").polygon(8, ROOF_ACROSS).extrude(ROOF_H)
         .translate((0.0, 0.0, -ROOF_H / 2.0)))
    cap = (cq.Workplane("XY").polygon(8, ROOF_ACROSS * 0.68).extrude(130.0)
           .translate((0.0, 0.0, ROOF_H / 2.0)))
    mast = (cq.Workplane("XY").circle(28.0).extrude(320.0)
            .translate((0.0, 0.0, ROOF_H / 2.0 + 130.0)))
    roof = r.union(cap).union(mast)
    try:
        roof = roof.edges(">Z").chamfer(24.0)
    except Exception:
        pass
    return roof


def make_ramp():
    """Front access ramp: hinge edge at local -RAMP_DEPTH/2, treads along Y,
    side rails for a walkway read."""
    deck = cq.Workplane("XY").box(RAMP_W, RAMP_DEPTH, RAMP_T)
    rail_l = (cq.Workplane("XY").box(60.0, RAMP_DEPTH, 90.0)
              .translate((RAMP_W / 2.0 - 30.0, 0.0, 30.0)))
    rail_r = (cq.Workplane("XY").box(60.0, RAMP_DEPTH, 90.0)
              .translate((-RAMP_W / 2.0 + 30.0, 0.0, 30.0)))
    ramp = deck.union(rail_l).union(rail_r)
    treads = None
    n = 6
    for i in range(n):
        yt = -RAMP_DEPTH / 2.0 + (i + 0.5) * RAMP_DEPTH / n
        t = (cq.Workplane("XY").box(RAMP_W - 140.0, 40.0, 24.0)
             .translate((0.0, yt, RAMP_T / 2.0 + 12.0)))
        treads = t if treads is None else treads.union(t)
    return ramp.union(treads)


body = make_body()
leg = make_leg()
foot = make_foot()
shell = make_shell()
roof = make_roof()
ramp = make_ramp()

# ── Local-frame parts, connector frames, and mates (declared from the plan) ──
PARTS = {
    "habitat_body": body,
    "landing_leg_fl": leg, "landing_leg_fr": leg,
    "landing_leg_rl": leg, "landing_leg_rr": leg,
    "level_foot_fl": foot,
    "side_shell_left": shell, "side_shell_right": shell,
    "telescoping_roof": roof,
    "access_ramp": ramp,
}

# body-side leg-hip mount points (reused for the deploy joints)
BLC = {
    "fl": [-LEG_OFF, LEG_OFF, LEG_Z], "fr": [LEG_OFF, LEG_OFF, LEG_Z],
    "rl": [-LEG_OFF, -LEG_OFF, LEG_Z], "rr": [LEG_OFF, -LEG_OFF, LEG_Z],
}

CONNECTORS = {
    "c_body_leg_fl": {"part": "habitat_body", "origin": BLC["fl"], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_body_leg_fr": {"part": "habitat_body", "origin": BLC["fr"], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_body_leg_rl": {"part": "habitat_body", "origin": BLC["rl"], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_body_leg_rr": {"part": "habitat_body", "origin": BLC["rr"], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_leg_fl_hip": {"part": "landing_leg_fl", "origin": [0, 0, LEG_HALF], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_leg_fr_hip": {"part": "landing_leg_fr", "origin": [0, 0, LEG_HALF], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_leg_rl_hip": {"part": "landing_leg_rl", "origin": [0, 0, LEG_HALF], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_leg_rr_hip": {"part": "landing_leg_rr", "origin": [0, 0, LEG_HALF], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_leg_fl_foot": {"part": "landing_leg_fl", "origin": [0, 0, -LEG_HALF], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_foot_fl_slide": {"part": "level_foot_fl", "origin": [0, 0, FOOT_CONN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_body_shell_left": {"part": "habitat_body", "origin": [-HALF, 0, 0], "z": [0, 1, 0], "x": [0, 0, 1]},
    "c_shell_left_hinge": {"part": "side_shell_left", "origin": [0, 0, SHELL_H / 2.0], "z": [0, 1, 0], "x": [0, 0, 1]},
    "c_body_shell_right": {"part": "habitat_body", "origin": [HALF, 0, 0], "z": [0, 1, 0], "x": [0, 0, 1]},
    "c_shell_right_hinge": {"part": "side_shell_right", "origin": [0, 0, SHELL_H / 2.0], "z": [0, 1, 0], "x": [0, 0, 1]},
    "c_body_roof": {"part": "habitat_body", "origin": [0, 0, HZ], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_roof_base": {"part": "telescoping_roof", "origin": [0, 0, -ROOF_H / 2.0], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_body_ramp": {"part": "habitat_body", "origin": [0, HALF, -100], "z": [1, 0, 0], "x": [0, 1, 0]},
    "c_ramp_hinge": {"part": "access_ramp", "origin": [0, -RAMP_DEPTH / 2.0, 0], "z": [1, 0, 0], "x": [0, 1, 0]},
}

MATES = [
    ["c_body_leg_fl", "c_leg_fl_hip"],
    ["c_body_leg_fr", "c_leg_fr_hip"],
    ["c_body_leg_rl", "c_leg_rl_hip"],
    ["c_body_leg_rr", "c_leg_rr_hip"],
    ["c_leg_fl_foot", "c_foot_fl_slide"],
    ["c_body_shell_left", "c_shell_left_hinge"],
    ["c_body_shell_right", "c_shell_right_hinge"],
    ["c_body_roof", "c_roof_base"],
    ["c_body_ramp", "c_ramp_hinge"],
]

_asm = yc_assemble(PARTS, CONNECTORS, MATES, "habitat_body")
PLACED = _asm[0]
XF = _asm[1]

# ── Assembly (placed parts, plan colors) ─────────────────────────────────────
ADD_LIST = [
    ("habitat_body", COL_BODY),
    ("landing_leg_fl", COL_LEG), ("landing_leg_fr", COL_LEG),
    ("landing_leg_rl", COL_LEG), ("landing_leg_rr", COL_LEG),
    ("level_foot_fl", COL_FOOT),
    ("side_shell_left", COL_SHELL), ("side_shell_right", COL_SHELL),
    ("telescoping_roof", COL_ROOF),
    ("access_ramp", COL_RAMP),
]
result = cq.Assembly()
for _nm, _col in ADD_LIST:
    result = result.add(PLACED[_nm], name=_nm, color=_col)

# ── Joints — origins/axes derived via yc_pt / yc_dir ─────────────────────────
JOINTS = [
    {"name": "leg_fl_deploy", "type": "revolute", "parts": ["landing_leg_fl", "habitat_body"],
     "origin": yc_pt(XF, "habitat_body", BLC["fl"]), "axis": yc_dir(XF, "habitat_body", [0, 1, 0]),
     "limits": [0.0, LEG_DEPLOY_DEG]},
    {"name": "leg_fr_deploy", "type": "revolute", "parts": ["landing_leg_fr", "habitat_body"],
     "origin": yc_pt(XF, "habitat_body", BLC["fr"]), "axis": yc_dir(XF, "habitat_body", [0, -1, 0]),
     "limits": [0.0, LEG_DEPLOY_DEG]},
    {"name": "leg_rl_deploy", "type": "revolute", "parts": ["landing_leg_rl", "habitat_body"],
     "origin": yc_pt(XF, "habitat_body", BLC["rl"]), "axis": yc_dir(XF, "habitat_body", [0, 1, 0]),
     "limits": [0.0, LEG_DEPLOY_DEG]},
    {"name": "leg_rr_deploy", "type": "revolute", "parts": ["landing_leg_rr", "habitat_body"],
     "origin": yc_pt(XF, "habitat_body", BLC["rr"]), "axis": yc_dir(XF, "habitat_body", [0, -1, 0]),
     "limits": [0.0, LEG_DEPLOY_DEG]},
    {"name": "level_slider_fl", "type": "slider", "parts": ["level_foot_fl", "landing_leg_fl"],
     "origin": yc_pt(XF, "landing_leg_fl", [0, 0, -LEG_HALF]),
     "axis": yc_dir(XF, "landing_leg_fl", [0, 0, 1]), "limits": [0.0, LEVEL_TRAVEL]},
    {"name": "shell_left_unfold", "type": "revolute", "parts": ["side_shell_left", "habitat_body"],
     "origin": yc_pt(XF, "habitat_body", [-HALF, 0, 0]), "axis": yc_dir(XF, "habitat_body", [0, 1, 0]),
     "limits": [0.0, SHELL_UNFOLD_DEG]},
    {"name": "shell_right_unfold", "type": "revolute", "parts": ["side_shell_right", "habitat_body"],
     "origin": yc_pt(XF, "habitat_body", [HALF, 0, 0]), "axis": yc_dir(XF, "habitat_body", [0, 1, 0]),
     "limits": [-SHELL_UNFOLD_DEG, 0.0]},
    {"name": "roof_raise", "type": "slider", "parts": ["telescoping_roof", "habitat_body"],
     "origin": yc_pt(XF, "habitat_body", [0, 0, HZ]), "axis": yc_dir(XF, "habitat_body", [0, 0, 1]),
     "limits": [0.0, ROOF_RAISE]},
    {"name": "ramp_lower", "type": "revolute", "parts": ["access_ramp", "habitat_body"],
     "origin": yc_pt(XF, "habitat_body", [0, HALF, -100]), "axis": yc_dir(XF, "habitat_body", [1, 0, 0]),
     "limits": [0.0, RAMP_LOWER_DEG]},
]
LINKS = []

# ── Animation: land → deploy legs → level → unfold shells → raise roof →
#    lower ramp → hold the workspace → stow everything back to rest ───────────
D = LEG_DEPLOY_DEG
L = LEVEL_TRAVEL
S = SHELL_UNFOLD_DEG
R = ROOF_RAISE
P = RAMP_LOWER_DEG
ANIMATION = {"fps": 24, "duration": 13.0, "tracks": [
    {"joint": "leg_fl_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [0.5, 0.0], [2.0, D], [9.0, D], [11.5, 0.0], [13.0, 0.0]]},
    {"joint": "leg_fr_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [0.5, 0.0], [2.0, D], [9.0, D], [11.5, 0.0], [13.0, 0.0]]},
    {"joint": "leg_rl_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [0.5, 0.0], [2.0, D], [9.0, D], [11.5, 0.0], [13.0, 0.0]]},
    {"joint": "leg_rr_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [0.5, 0.0], [2.0, D], [9.0, D], [11.5, 0.0], [13.0, 0.0]]},
    {"joint": "level_slider_fl", "unit": "mm",
     "keys": [[0.0, 0.0], [2.0, 0.0], [3.0, L], [9.0, L], [11.0, 0.0], [13.0, 0.0]]},
    {"joint": "shell_left_unfold", "unit": "deg",
     "keys": [[0.0, 0.0], [3.0, 0.0], [4.5, S], [9.0, S], [10.5, 0.0], [13.0, 0.0]]},
    {"joint": "shell_right_unfold", "unit": "deg",
     "keys": [[0.0, 0.0], [3.0, 0.0], [4.5, -S], [9.0, -S], [10.5, 0.0], [13.0, 0.0]]},
    {"joint": "roof_raise", "unit": "mm",
     "keys": [[0.0, 0.0], [4.5, 0.0], [5.8, R], [9.0, R], [10.0, 0.0], [13.0, 0.0]]},
    {"joint": "ramp_lower", "unit": "deg",
     "keys": [[0.0, 0.0], [5.8, 0.0], [7.0, P], [9.0, P], [9.8, 0.0], [13.0, 0.0]]},
]}
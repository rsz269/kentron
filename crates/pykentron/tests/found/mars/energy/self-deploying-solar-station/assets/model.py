"""Compact Mars solar-power station — an unbranded industrial product concept.

Arrives as a closed rectangular lander module. From it deploy: four landing
feet that swing down to plant a wide stance, a telescoping central mast that
rises from the deck center, two inner solar wings that hinge from the mast tip,
and their photovoltaic panels that roll to a sun-tracking orientation. Every
part is built in its own local frame and placed by the deterministic assembler
from the supplied connector plan (yc_assemble); joints are derived from the
placed frames with yc_pt / yc_dir. This is a concept, not a validated power
system. Left/right wings and front/rear feet share single dimension constants
so edits stay symmetric.
"""
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
BODY_L = 1200.0          # module body length X [800:1600]
BODY_W = 800.0           # module body width Y [500:1100]
BODY_H = 500.0           # module body height Z [300:700]
MAST_H = 600.0           # collapsed mast height [400:900]
MAST_TRAVEL = 600.0      # mast rise travel [200:900]
WING_LEN = 700.0         # inner wing spar length [400:1000]
PANEL_L = 1200.0         # solar panel length [600:1600]
FOOT_DEPLOY_DEG = 90.0   # landing-foot deploy angle [30:90]
WING_FOLD_DEG = 55.0     # wing hinge demo angle [20:80]
PANEL_TRACK_DEG = 60.0   # panel sun-track sweep [20:90]

# ── Fixed proportions (shared so mirrored parts can't desync) ────────────────
MAST_SIDE = 180.0                     # mast square cross-section
MAST_WING_X = 70.0                    # wing-hinge offset from mast center
WING_W = 120.0                        # wing spar width
WING_T = 30.0                         # wing spar thickness
PANEL_W = 600.0                       # panel width Y
PANEL_T = 20.0                        # panel thickness
FOOT_LEN = 350.0                      # landing-leg length
FOOT_W = 120.0                        # landing-leg width
FOOT_T = 80.0                         # landing-leg thickness
FOOT_PAD = 160.0                      # landing-pad size
FOOT_HINGE_Z = -(BODY_H / 2.0 - 50.0) # foot hinge height on the body side

# ── Colors ───────────────────────────────────────────────────────────────────
COL_BODY = cq.Color(0.62, 0.63, 0.64)
COL_MAST = cq.Color(0.45, 0.46, 0.48)
COL_WING = cq.Color(0.55, 0.56, 0.58)
COL_PANEL = cq.Color(0.11, 0.16, 0.36)
COL_FOOT = cq.Color(0.30, 0.31, 0.33)


# ── Part builders (each in its OWN local frame, Z up) ────────────────────────
def build_body():
    """Closed carcass centered on origin; top at +BODY_H/2, bottom the ground
    footprint at -BODY_H/2. A low collar around the mast well adds an
    industrial deck detail."""
    b = cq.Workplane("XY").box(BODY_L, BODY_W, BODY_H)
    try:
        b = b.edges("|Z").fillet(24.0)
    except Exception:
        pass
    # mast well collar on the top deck
    collar = (cq.Workplane("XY").rect(MAST_SIDE + 60.0, MAST_SIDE + 60.0)
              .rect(MAST_SIDE + 12.0, MAST_SIDE + 12.0).extrude(28.0)
              .translate((0.0, 0.0, BODY_H / 2.0)))
    b = b.union(collar)
    # shallow equipment recess band along the sides
    band = (cq.Workplane("XY").box(BODY_L - 120.0, BODY_W + 4.0, 60.0)
            .translate((0.0, 0.0, BODY_H / 2.0 - 90.0)))
    try:
        b = b.cut(band.faces("<Y").workplane().rect(BODY_L - 120.0, 60.0)
                  .extrude(-8.0))
    except Exception:
        pass
    return b


def build_mast():
    """Telescoping mast: wide base stage plus a narrower upper stage for a
    two-segment look. Base at local z=0, top at z=MAST_H."""
    base = (cq.Workplane("XY").box(MAST_SIDE, MAST_SIDE, MAST_H)
            .translate((0.0, 0.0, MAST_H / 2.0)))
    upper = (cq.Workplane("XY").box(MAST_SIDE * 0.68, MAST_SIDE * 0.68, MAST_H * 0.55)
             .translate((0.0, 0.0, MAST_H * 0.78)))
    m = base.union(upper)
    # hinge lugs at the top where the two wings mount
    for sx in (1.0, -1.0):
        lug = (cq.Workplane("XY").box(60.0, WING_W + 24.0, 60.0)
               .translate((sx * (MAST_WING_X), 0.0, MAST_H - 20.0)))
        m = m.union(lug)
    try:
        m = m.edges("|Z").fillet(8.0)
    except Exception:
        pass
    return m


def build_wing():
    """Inner wing spar from root (local origin) out to the tip at +X."""
    w = (cq.Workplane("XY").box(WING_LEN, WING_W, WING_T)
         .translate((WING_LEN / 2.0, 0.0, 0.0)))
    # lightening pockets so it reads as a truss spar
    try:
        for i in range(3):
            xc = WING_LEN * (0.28 + 0.22 * i)
            w = w.cut(cq.Workplane("XY").box(90.0, WING_W - 44.0, WING_T + 4.0)
                      .translate((xc, 0.0, 0.0)))
    except Exception:
        pass
    return w


def build_panel():
    """Photovoltaic panel; root (local origin) at the wing tip, extends +X.
    A raised border frame and a few groove lines suggest cell strings."""
    p = (cq.Workplane("XY").box(PANEL_L, PANEL_W, PANEL_T)
         .translate((PANEL_L / 2.0, 0.0, 0.0)))
    try:
        # engrave a light cell grid on the top face
        for i in range(1, 4):
            yy = -PANEL_W / 2.0 + PANEL_W * i / 4.0
            p = p.cut(cq.Workplane("XY").box(PANEL_L - 30.0, 4.0, 3.0)
                      .translate((PANEL_L / 2.0, yy, PANEL_T / 2.0)))
        for j in range(1, 6):
            xx = PANEL_L * j / 6.0
            p = p.cut(cq.Workplane("XY").box(4.0, PANEL_W - 30.0, 3.0)
                      .translate((xx, 0.0, PANEL_T / 2.0)))
    except Exception:
        pass
    return p


def build_foot():
    """Landing leg: hinge at local origin, leg out along +X with a pad at the
    far end. Stowed at 0 (leg horizontal, retracted); deploying rotates it
    down to plant."""
    leg = (cq.Workplane("XY").box(FOOT_LEN, FOOT_W, FOOT_T)
           .translate((FOOT_LEN / 2.0, 0.0, 0.0)))
    pad = (cq.Workplane("XY").box(FOOT_PAD, FOOT_PAD * 1.3, 55.0)
           .translate((FOOT_LEN, 0.0, 0.0)))
    f = leg.union(pad)
    try:
        f = f.edges("|Y").fillet(12.0)
    except Exception:
        pass
    return f


# ── Local-frame solids ───────────────────────────────────────────────────────
PARTS = {
    "module_body": build_body(),
    "central_mast": build_mast(),
    "wing_inner_left": build_wing(),
    "wing_inner_right": build_wing(),
    "panel_left": build_panel(),
    "panel_right": build_panel(),
    "foot_front": build_foot(),
    "foot_rear": build_foot(),
}

# ── Connector frames (straight from the plan, expressed via parameters) ──────
CONNECTORS = {
    "c_body_mast": {"part": "module_body", "origin": [0.0, 0.0, BODY_H / 2.0],
                    "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_body_foot_front": {"part": "module_body", "origin": [BODY_L / 2.0, 0.0, FOOT_HINGE_Z],
                          "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_body_foot_rear": {"part": "module_body", "origin": [-BODY_L / 2.0, 0.0, FOOT_HINGE_Z],
                         "z": [0, 0, 1], "x": [-1, 0, 0]},
    "c_mast_base": {"part": "central_mast", "origin": [0.0, 0.0, 0.0],
                    "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_mast_wing_left": {"part": "central_mast", "origin": [MAST_WING_X, 0.0, MAST_H],
                         "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_mast_wing_right": {"part": "central_mast", "origin": [-MAST_WING_X, 0.0, MAST_H],
                          "z": [0, 0, 1], "x": [-1, 0, 0]},
    "c_wing_left_root": {"part": "wing_inner_left", "origin": [0.0, 0.0, 0.0],
                         "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_wing_left_tip": {"part": "wing_inner_left", "origin": [WING_LEN, 0.0, 0.0],
                        "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_wing_right_root": {"part": "wing_inner_right", "origin": [0.0, 0.0, 0.0],
                          "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_wing_right_tip": {"part": "wing_inner_right", "origin": [WING_LEN, 0.0, 0.0],
                         "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_panel_left_root": {"part": "panel_left", "origin": [0.0, 0.0, 0.0],
                          "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_panel_right_root": {"part": "panel_right", "origin": [0.0, 0.0, 0.0],
                           "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_foot_front_hinge": {"part": "foot_front", "origin": [0.0, 0.0, 0.0],
                           "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_foot_rear_hinge": {"part": "foot_rear", "origin": [0.0, 0.0, 0.0],
                          "z": [0, 0, 1], "x": [1, 0, 0]},
}

# ── Mates (child connector ↔ parent connector) ───────────────────────────────
MATES = [
    ["c_body_mast", "c_mast_base"],
    ["c_mast_wing_left", "c_wing_left_root"],
    ["c_mast_wing_right", "c_wing_right_root"],
    ["c_wing_left_tip", "c_panel_left_root"],
    ["c_wing_right_tip", "c_panel_right_root"],
    ["c_body_foot_front", "c_foot_front_hinge"],
    ["c_body_foot_rear", "c_foot_rear_hinge"],
]

# ── Deterministic placement (module_body is the fixed root/ground) ──────────
_asm = yc_assemble(PARTS, CONNECTORS, MATES, "module_body")
PLACED = _asm[0]
XF = _asm[1]

result = (cq.Assembly()
          .add(PLACED["module_body"], name="module_body", color=COL_BODY)
          .add(PLACED["central_mast"], name="central_mast", color=COL_MAST)
          .add(PLACED["wing_inner_left"], name="wing_inner_left", color=COL_WING)
          .add(PLACED["wing_inner_right"], name="wing_inner_right", color=COL_WING)
          .add(PLACED["panel_left"], name="panel_left", color=COL_PANEL)
          .add(PLACED["panel_right"], name="panel_right", color=COL_PANEL)
          .add(PLACED["foot_front"], name="foot_front", color=COL_FOOT)
          .add(PLACED["foot_rear"], name="foot_rear", color=COL_FOOT))

# ── Joints — every origin/axis derived from the placed frames ───────────────
# Mirrored parts get opposite axis signs so a single positive animation value
# drives left/right (and front/rear) symmetrically.
JOINTS = [
    {"name": "mast_rise", "type": "slider",
     "parts": ["central_mast", "module_body"],
     "origin": yc_pt(XF, "module_body", [0.0, 0.0, BODY_H / 2.0]),
     "axis": yc_dir(XF, "module_body", [0, 0, 1]),
     "limits": [0.0, MAST_TRAVEL]},

    {"name": "wing_left_unfold", "type": "revolute",
     "parts": ["wing_inner_left", "central_mast"],
     "origin": yc_pt(XF, "central_mast", [MAST_WING_X, 0.0, MAST_H]),
     "axis": yc_dir(XF, "central_mast", [0, 1, 0]),
     "limits": [0.0, 90.0]},

    {"name": "wing_right_unfold", "type": "revolute",
     "parts": ["wing_inner_right", "central_mast"],
     "origin": yc_pt(XF, "central_mast", [-MAST_WING_X, 0.0, MAST_H]),
     "axis": yc_dir(XF, "central_mast", [0, -1, 0]),
     "limits": [0.0, 90.0]},

    {"name": "panel_left_track", "type": "revolute",
     "parts": ["panel_left", "wing_inner_left"],
     "origin": yc_pt(XF, "wing_inner_left", [WING_LEN, 0.0, 0.0]),
     "axis": yc_dir(XF, "wing_inner_left", [1, 0, 0]),
     "limits": [-90.0, 90.0]},

    {"name": "panel_right_track", "type": "revolute",
     "parts": ["panel_right", "wing_inner_right"],
     "origin": yc_pt(XF, "wing_inner_right", [WING_LEN, 0.0, 0.0]),
     "axis": yc_dir(XF, "wing_inner_right", [1, 0, 0]),
     "limits": [-90.0, 90.0]},

    {"name": "foot_front_deploy", "type": "revolute",
     "parts": ["foot_front", "module_body"],
     "origin": yc_pt(XF, "module_body", [BODY_L / 2.0, 0.0, FOOT_HINGE_Z]),
     "axis": yc_dir(XF, "module_body", [0, 1, 0]),
     "limits": [0.0, 90.0]},

    {"name": "foot_rear_deploy", "type": "revolute",
     "parts": ["foot_rear", "module_body"],
     "origin": yc_pt(XF, "module_body", [-BODY_L / 2.0, 0.0, FOOT_HINGE_Z]),
     "axis": yc_dir(XF, "module_body", [0, -1, 0]),
     "limits": [0.0, 90.0]},
]
LINKS = []

# ── Deployment timeline ──────────────────────────────────────────────────────
# Staged so the folding sequence reads clearly and every joint returns to rest:
#   feet swing down to plant (1-3s) → mast telescopes up (3-5s) →
#   wings hinge open/closed while the mast is clear of the deck (5-11s) →
#   panels roll to sun-track (11-16s) → everything retracts back to rest.
# The wing hinge only moves while the mast is fully raised, so the long panels
# never sweep below the ground plane.
ANIMATION = {"fps": 24, "duration": 20.0, "tracks": [
    {"joint": "foot_front_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [1.0, 0.0], [3.0, FOOT_DEPLOY_DEG],
              [17.0, FOOT_DEPLOY_DEG], [19.0, 0.0], [20.0, 0.0]]},
    {"joint": "foot_rear_deploy", "unit": "deg",
     "keys": [[0.0, 0.0], [1.0, 0.0], [3.0, FOOT_DEPLOY_DEG],
              [17.0, FOOT_DEPLOY_DEG], [19.0, 0.0], [20.0, 0.0]]},

    {"joint": "mast_rise", "unit": "mm",
     "keys": [[0.0, 0.0], [3.0, 0.0], [5.0, MAST_TRAVEL],
              [17.0, MAST_TRAVEL], [19.0, 0.0], [20.0, 0.0]]},

    {"joint": "wing_left_unfold", "unit": "deg",
     "keys": [[0.0, 0.0], [5.0, 0.0], [7.0, WING_FOLD_DEG],
              [9.0, WING_FOLD_DEG], [11.0, 0.0], [20.0, 0.0]]},
    {"joint": "wing_right_unfold", "unit": "deg",
     "keys": [[0.0, 0.0], [5.0, 0.0], [7.0, WING_FOLD_DEG],
              [9.0, WING_FOLD_DEG], [11.0, 0.0], [20.0, 0.0]]},

    {"joint": "panel_left_track", "unit": "deg",
     "keys": [[0.0, 0.0], [11.0, 0.0], [13.0, PANEL_TRACK_DEG],
              [14.0, -PANEL_TRACK_DEG], [15.0, PANEL_TRACK_DEG],
              [16.0, 0.0], [20.0, 0.0]]},
    {"joint": "panel_right_track", "unit": "deg",
     "keys": [[0.0, 0.0], [11.0, 0.0], [13.0, -PANEL_TRACK_DEG],
              [14.0, PANEL_TRACK_DEG], [15.0, -PANEL_TRACK_DEG],
              [16.0, 0.0], [20.0, 0.0]]},
]}
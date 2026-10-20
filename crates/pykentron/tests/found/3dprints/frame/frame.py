"""Picture frame with plain solid sides, for a sheet of
glass, a picture and a backing board, in as few sections as fit the bed.

Run `uv run frame.py` from this folder to write STEP and STL files to out/.

The frame is modelled face down, which is also how it prints: z = 0 is the
front and z grows towards the back. It is centred on the origin with its long
side along x, and every section is modelled where it sits in the frame. The
contents go in from the back and are held by turn buttons on pegs.
"""

import math
from dataclasses import dataclass
from pathlib import Path

import cadquery as cq


@dataclass(frozen=True)
class Params:
    # Contents: glass, picture and backing board, as one stack. content_w is
    # the long side.
    content_w: float = 279.4  # 11 in
    content_h: float = 215.9  # 8.5 in
    stack: float = 6.5  # at the edge, where the buttons bear: measured; the folded easel is thicker, in the middle
    clr: float = 0.4  # per side, between the contents and the rabbet wall

    # Border, from the edge of the opening to the outside
    border: float = 21.0
    lip: float = 4.5  # how far the front lip overlaps the contents
    skin: float = 2.0  # thickness of the front lip
    wall: float = 1.8  # rabbet wall
    rim: float = 1.6  # outer rim
    depth: float = 13.5
    chamfer: float = 0.6  # 45 degrees, round the front edges inside and out
    corner_r: float = 1.0  # outside corners seen from the front

    pad: float = 20.0  # length of the seam fit-test pieces

    # Support for the lip while printing face up: a wall one line thick under
    # it, round the opening, peeled off before use, held by tabs at its foot
    sup_in: float = 0.8  # from the opening's edge, so its scar is hidden
    sup_w: float = 0.5
    sup_gap: float = 0.2  # under the lip, like a slicer support's top gap
    tab_w: float = 2.0
    tab_h: float = 0.4  # two layers
    tab_every: float = 50.0

    bed: float = 250.0  # largest part either way; the P1S bed is 256

    # Sliding dovetail at each seam, put in from the back. Its front is a 45
    # degree ramp that rests on a matching ramp in the socket.
    tongue_neck: float = 6.0
    tongue_tip: float = 10.0
    tongue_len: float = 7.0
    ramp_z: float = 1.0  # where the ramp meets the seam; the socket's thinnest front
    joint_clr: float = 0.0  # per side, when seated; the taper gives room going in
    joint_taper: float = 1.0  # degrees, narrower towards the front

    # Turn buttons: a flat bar on a peg in a pocket in the back of a pad,
    # swung over the back of the backing board
    slack: float = 0.0  # between the back of the stack and the buttons
    btn_t: float = 3.0
    btn_lap: float = 5.0  # how far a closed button reaches over the contents
    btn_pivot: float = 3.6  # peg, behind the face of the rabbet wall; at least btn_half, so it opens clear
    btn_off: float = 3.5  # peg to the hole that lets the button off over its head
    btn_rim: float = 1.2  # round the holes
    btn_chamfer: float = 0.5  # tall, under the edges, to ride onto a thick backing
    btn_lead: float = 0.3  # its inset: steeper than 48 degrees, so it prints on the bed
    btn_clr: float = 0.1  # radial, round the peg
    peg_r: float = 1.5
    peg_head: float = 0.8  # head overhang; its underside is 50 degrees, a margin over 45
    pocket_clr: float = 0.5
    pocket_drop: float = 0.4  # pocket floor ahead of the stack's back: it prints as a bridge, which may sag
    # Detent: a finger on the button's tail curves round inside the pocket
    # wall and presses a nub on it, which clicks into notches at closed and open
    det_w: float = 1.1  # finger, bending in the plane of the button
    det_h: float = 4.0  # finger and nub, above the 3 mm bar: stiffer at the same strain
    det_press: float = 0.7  # how far the nub reaches past the wall: 0.4 was lost in the print's tolerances
    nub_r: float = 1.0  # at least det_press, so the nub ramps out of a notch rather than hooking
    notch_r: float = 1.1

    # Keyholes for a picture nail or a small screw
    key_head: float = 6.0
    key_neck: float = 3.0
    key_throw: float = 5.0  # entry centre to slot end, towards the outside edge
    key_lip: float = 1.5
    key_room: float = 3.5  # under the lips, including their 45 degree undersides
    key_web: float = 1.2  # walls between a keyhole and the rabbet or a pocket


P = Params()


# --- Sizes --------------------------------------------------------------------


def rabbet(p=P):
    return p.content_w + 2 * p.clr, p.content_h + 2 * p.clr


def sight(p=P):
    return p.content_w - 2 * p.lip, p.content_h - 2 * p.lip


def outer(p=P):
    w, h = sight(p)
    return w + 2 * p.border, h + 2 * p.border


def band(p=P):
    return p.border - p.lip - p.clr - p.wall - p.rim


def stack_back(p=P):
    """Floor of the button pockets: the back of the stack, plus slack."""
    return p.skin + p.stack + p.slack


def pocket_floor(p=P):
    return stack_back(p) - p.pocket_drop


def btn_arm(p=P):
    """Peg to the tip of a button."""
    return p.btn_pivot + p.clr + p.btn_lap


def btn_hole(p=P):
    """Radius of the hole that passes over the peg's head."""
    return p.peg_r + p.peg_head + p.btn_clr


def btn_half(p=P):
    return btn_hole(p) + p.btn_rim


def tail_half(p=P):
    """Half width of the button's tail, narrower than the arm so the detent
    finger leaves it sooner and is longer."""
    return p.peg_r + p.btn_clr + p.btn_rim


def btn_tail(p=P):
    """Peg to the back end of a button, which reaches almost to the pocket
    wall. That stops it sliding back as far as btn_off, which would free it,
    and carries the detent finger."""
    return pocket_r(p) - 0.1


def pocket_r(p=P):
    return btn_arm(p) + p.pocket_clr


def btn_x(p=P):
    """Peg along the side from the keyhole, clear of it by key_web."""
    return p.key_head / 2 + p.key_web + pocket_r(p)


def notch_c(p=P):
    """Peg to the centre of a detent notch in the pocket wall."""
    return pocket_r(p) + p.det_press - p.nub_r + 0.5


def lock_span(p=P):
    """A lock's keyhole, pocket and notches along the side, from its keyhole:
    (start, end)."""
    return -(p.key_head / 2 + p.key_web), btn_x(p) + notch_c(p) + p.notch_r + p.key_web


def det_len(p=P):
    """Free length of the detent finger, from where it leaves the tail to the nub."""
    r = btn_tail(p) - p.det_w / 2
    return r * math.acos(tail_half(p) / r)


def key_entry(yr, p=P):
    """y of a keyhole's entry, its slot running out towards the edge."""
    return yr + p.wall + p.key_web + p.key_head / 2


def sections(p=P):
    """How many sections along x the frame is cut into to fit the bed."""
    w, h = outer(p)
    if p.content_h > p.content_w:
        raise ValueError("content_w must be the long side")
    if h > p.bed:
        raise ValueError(f"frame is {h:.1f} mm on its short side, more than the {p.bed} mm bed")
    n = 1
    while w / n + (p.tongue_len if n > 1 else 0) > p.bed:
        n += 1
    return n


def seams(p=P):
    w, _ = outer(p)
    n = sections(p)
    return [-w / 2 + k * w / n for k in range(1, n)]


def lock_spots(p=P):
    """x of the lock pads on the long sides: the middle of each section."""
    w, _ = outer(p)
    n = sections(p)
    return [-w / 2 + (k + 0.5) * w / n for k in range(n)]


# --- Solid-building helpers ---------------------------------------------------


def box(x0, x1, y0, y1, z0, z1):
    return cq.Workplane("XY").box(x1 - x0, y1 - y0, z1 - z0, centered=False).translate((x0, y0, z0))


def cyl(r, z0, z1, x=0.0, y=0.0):
    return cq.Workplane("XY").workplane(offset=z0).center(x, y).circle(r).extrude(z1 - z0)


def rprism(w, d, r, z0, z1):
    """Rounded-rectangle prism centred on the z axis."""
    s = cq.Workplane("XY").workplane(offset=z0).rect(w, d).extrude(z1 - z0)
    return s.edges("|Z").fillet(r) if r > 0 else s


def _rrect(w, d, r):
    s = cq.Sketch().rect(w, d)
    return s.vertices().fillet(r) if r > 0 else s


def flare(w, d, r, z0, h, grow):
    """Rounded rectangle (w, d, r) at z0, offset outward by `grow` at z0 + h."""
    a = _rrect(w, d, r).moved(cq.Location(cq.Vector(0, 0, z0)))
    b = _rrect(w + 2 * grow, d + 2 * grow, r + grow if r > 0 else 0).moved(cq.Location(cq.Vector(0, 0, z0 + h)))
    return cq.Workplane("XY").placeSketch(a, b).loft(ruled=True)


def grow(pts, d):
    """A convex polygon with every edge moved outward by d (inward if d < 0)."""
    if sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1])) < 0:
        pts = pts[::-1]  # counter-clockwise, so the right of each edge is outside
    lines = []
    for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]):
        length = math.hypot(x1 - x0, y1 - y0)
        nx, ny = (y1 - y0) / length, -(x1 - x0) / length
        lines.append(((x0 + nx * d, y0 + ny * d), (x1 - x0, y1 - y0)))
    out = []
    for ((ax, ay), (adx, ady)), ((bx, by), (bdx, bdy)) in zip(lines[-1:] + lines[:-1], lines):
        t = ((bx - ax) * bdy - (by - ay) * bdx) / (adx * bdy - ady * bdx)
        out.append((ax + adx * t, ay + ady * t))
    return out


def loft(pts0, z0, pts1, z1):
    """Ruled solid from polygon pts0 at z0 to polygon pts1 at z1."""
    w = cq.Workplane("XY").workplane(offset=z0).polyline(pts0).close()
    return w.workplane(offset=z1 - z0).polyline(pts1).close().loft(ruled=True)


def prism(pts, z0, z1):
    return loft(pts, z0, pts, z1)


def combo(*parts):
    return cq.Workplane(obj=cq.Compound.makeCompound([s for p in parts for s in p.solids().vals()]))


def cut_all(part, cutters):
    """Cut many solids in one boolean. They must not overlap one another."""
    return part.cut(combo(*cutters)) if cutters else part


def turn(part, deg, x=0.0, y=0.0):
    return part.rotate((x, y, 0), (x, y, 1), deg)


def face_up(part, p=P):
    """A part turned over to print front up: z = 0 is then the back."""
    return part.rotate((0, 0, 0), (1, 0, 0), 180).translate((0, 0, p.depth))


# Features are built on the top side (+y) and turned onto the others. Each side
# is (name, turn in degrees, distance of the rabbet wall from the centre).
def sides(p=P):
    rw, rh = rabbet(p)
    return [("top", 0, rh / 2), ("bottom", 180, rh / 2), ("right", -90, rw / 2), ("left", 90, rw / 2)]


def locks(p=P):
    """(turn, yr, position) of every lock: a keyhole at the position and its
    button to the right of it. One in the middle of each section on the long
    sides, one in the middle of each short side."""
    return [(deg, yr, c) for side, deg, yr in sides(p) for c in (lock_spots(p) if side in ("top", "bottom") else [0.0])]


# --- Keyholes and turn buttons -------------------------------------------------


def _stadium(c, y0, length, width, z):
    return cq.Workplane("XY").workplane(offset=z).center(c, y0 + length / 2).slot2D(length + width, width, 90)


def keyhole_cutter(c, yr, p=P):
    """Keyhole in the back of a pad, its slot running out towards the edge."""
    y0 = key_entry(yr, p)
    zl = p.depth - p.key_lip  # underside of the lips
    zs = zl - (p.key_head - p.key_neck) / 2  # where the 45 degree underside starts
    zf = zl - p.key_room  # floor
    th, hd, nk = p.key_throw, p.key_head, p.key_neck
    entry = cyl(hd / 2, zf, p.depth + 1, c, y0)
    neck = _stadium(c, y0, th, nk, zl - 0.01).extrude(p.depth + 1 - zl)
    # Room for the head under the lips, as one loft so there is no seam where
    # the straight sides meet the 45 degree undersides. Face up its floor is a
    # ceiling, so the loft starts from a 45 degree gabled roof.
    zr = zf - hd / 2 + 0.01  # ridge
    room = _stadium(c, y0, th, 0.02, zr).workplane(offset=zf - zr).slot2D(th + hd, hd, 90)
    room = room.workplane(offset=zs - zf).slot2D(th + hd, hd, 90)
    room = room.workplane(offset=zl - zs).slot2D(th + nk, nk, 90).loft(ruled=True)
    return entry.union(neck).union(room)


def peg(x, y, p=P):
    """The pivot for a button: a neck and a head whose top is flush with the
    back, so face up it prints from the bed, head first. The head's underside
    is 50 degrees."""
    zn = stack_back(p) + p.btn_t + 0.2  # 0.2 under the head
    rh, rise = p.peg_r + p.peg_head, p.peg_head * math.tan(math.radians(50))
    neck = cyl(p.peg_r, pocket_floor(p) - 0.01, zn + 0.01, x, y)
    # One loft for the underside and the top, so there is no seam between them
    rings = [cq.Sketch().circle(r).moved(cq.Location(cq.Vector(x, y, z))) for r, z in ((p.peg_r, zn), (rh, zn + rise), (rh, p.depth))]
    return neck.union(cq.Workplane("XY").placeSketch(*rings).loft(ruled=True))


def _bar(y0, y1, hw, p=P):
    """A rounded bar along y from y0 to y1, its lower edges chamfered."""
    ch, ld, length, yc = p.btn_chamfer, p.btn_lead, y1 - y0, (y0 + y1) / 2
    lead = cq.Workplane("XY").center(0, yc).slot2D(length - 2 * ld, 2 * (hw - ld), 90)
    lead = lead.workplane(offset=ch).slot2D(length, 2 * hw, 90).loft(ruled=True)
    bar = cq.Workplane("XY").workplane(offset=ch - 0.01).center(0, yc).slot2D(length, 2 * hw, 90).extrude(p.btn_t - ch + 0.01)
    return lead.union(bar)


def button(p=P, finger=True):
    """A turn button as fitted, closed on the top side: peg at the origin, arm
    towards -y, underside at z = 0. It goes on over the peg's head by the hole
    in its arm, then slides towards its tip so the peg sits at the end of the
    slot. To come off it has to slide back the other way, which the pocket
    wall stops unless it is turned with its arm pointing out from the opening.
    The detent finger runs a quarter turn round from the tail, its nub on +x,
    and stands det_h tall, above the bar."""
    hw, tw, arm, tail = btn_half(p), tail_half(p), btn_arm(p), btn_tail(p)
    # The arm runs from its tip to round the peg, the narrower tail from the peg back
    part = _bar(-arm, hw, hw, p).union(_bar(-tw, tail, tw, p))
    holes = _stadium(0, -p.btn_off, p.btn_off, 2 * (p.peg_r + p.btn_clr), -1).extrude(p.btn_t + 2)
    holes = holes.union(cyl(btn_hole(p), -1, p.btn_t + 1, 0, -p.btn_off))
    part = part.cut(holes)
    if finger:
        ro = btn_tail(p)
        ring = cyl(ro, 0, p.det_h).cut(cyl(ro - p.det_w, -1, p.det_h + 1))
        ring = ring.intersect(box(0, ro + 1, 0, ro + 1, -1, p.det_h + 1))
        nub = cyl(p.nub_r, 0, p.det_h, pocket_r(p) + p.det_press - p.nub_r)
        part = part.union(ring).union(nub)
    return part


def button_print(p=P):
    """The button as fitted prints on its underside, the finger standing up."""
    return button(p)


def button_at(deg, yr, c, theta=0.0, p=P, part=None):
    """A button fitted at a lock pad, turned theta degrees from closed:
    90 is open, 180 is the pose it comes off in."""
    part = part if part is not None else button(p)
    part = part.rotate((0, 0, 0), (0, 0, 1), theta).translate((c + btn_x(p), yr + p.btn_pivot, stack_back(p)))
    return turn(part, deg)


def buttons(p=P):
    """How many buttons the frame takes: one per lock pad."""
    return len(locks(p))


def lock_features(p=P):
    """(cutters, pegs) for every lock pad: its keyhole, and its button's
    pocket and peg."""
    cutters, pegs = [], []
    for deg, yr, c in locks(p):
        x, y = c + btn_x(p), yr + p.btn_pivot
        pocket = cyl(pocket_r(p), pocket_floor(p), p.depth + 1, x, y)
        for a in (0, 90):  # notches for the nub, closed and open
            a = math.radians(a)
            pocket = pocket.union(cyl(p.notch_r, pocket_floor(p), p.depth + 1, x + notch_c(p) * math.cos(a), y + notch_c(p) * math.sin(a)))
        cutters.append(turn(keyhole_cutter(c, yr, p).union(pocket), deg))
        pegs.append(turn(peg(x, y, p), deg))
    return cutters, pegs


def lip_support(p=P):
    """The support wall under the lip, in frame coordinates, with tabs at its
    foot joining it to the rabbet wall. The tabs keep clear of the seams and of
    the button pockets, which cut the rabbet wall away."""
    sw, sh = sight(p)
    a, b = p.sup_in, p.sup_in + p.sup_w
    z0 = p.skin + p.sup_gap
    wall = box(-sw / 2 - b, sw / 2 + b, -sh / 2 - b, sh / 2 + b, z0, p.depth)
    wall = wall.cut(box(-sw / 2 - a, sw / 2 + a, -sh / 2 - a, sh / 2 + a, z0 - 1, p.depth + 1))
    tabs = []
    for side, deg, yr in sides(p):
        along = sw / 2 if side in ("top", "bottom") else sh / 2
        pockets = [(c + btn_x(p) - pocket_r(p), c + btn_x(p) + pocket_r(p)) for d, _, c in locks(p) if d == deg]
        n = int((along - 5) // p.tab_every)
        for k in range(-n, n + 1):
            u = k * p.tab_every
            if side in ("top", "bottom") and any(abs(u - s) < 5 for s in seams(p)):
                continue
            if any(lo - p.tab_w < u < hi + p.tab_w for lo, hi in pockets):
                continue
            y0 = yr - p.clr - p.lip + b - 0.01  # the wall's outer face
            tabs.append(turn(box(u - p.tab_w / 2, u + p.tab_w / 2, y0, yr + 0.01, p.depth - p.tab_h, p.depth), deg))
    return wall.union(combo(*tabs))


# --- Frame ----------------------------------------------------------------------


def make_body(p=P):
    """The whole frame as one piece, before it is cut into sections."""
    w, h = outer(p)
    sw, sh = sight(p)
    rw, rh = rabbet(p)
    c, r = p.chamfer, p.corner_r
    # One ruled loft for the chamfer and the sides, so there is no seam between them
    rings = [_rrect(w - 2 * c, h - 2 * c, r - c), _rrect(w, h, r), _rrect(w, h, r)]
    rings = [s.moved(cq.Location(cq.Vector(0, 0, z))) for s, z in zip(rings, (0, c, p.depth))]
    body = cq.Workplane("XY").placeSketch(*rings).loft(ruled=True)
    body = body.cut(box(-rw / 2, rw / 2, -rh / 2, rh / 2, p.skin, p.depth + 1))
    body = body.cut(box(-sw / 2, sw / 2, -sh / 2, sh / 2, -1, p.skin + 0.1))
    body = body.cut(flare(sw + 2 * c, sh + 2 * c, 0, -0.01, c + 0.01, -c - 0.01))
    cutters, pegs = lock_features(p)
    body = cut_all(body, cutters)
    for pg in pegs:
        body = body.union(pg)
    return body


def tongue(xs, p=P, clr=0.0, extra=0.0):
    """The dovetail on the top side at the seam x = xs, pointing +x. With clr
    and extra it is the socket cutter: grown by clr and run out past the back."""
    _, rh = rabbet(p)
    yc = rh / 2 + (p.wall + band(p) + p.rim) / 2
    n, t = p.tongue_neck / 2, p.tongue_tip / 2
    pts = [(xs - 1, yc - n), (xs, yc - n), (xs + p.tongue_len, yc - t), (xs + p.tongue_len, yc + t), (xs, yc + n), (xs - 1, yc + n)]
    z0, z1 = p.ramp_z - 0.5, p.depth + extra
    taper = math.tan(math.radians(p.joint_taper))
    # The taper is measured from the seated front of the tongue
    flanks = loft(grow(pts, clr - taper * (p.ramp_z - z0)), z0, grow(pts, clr + taper * (z1 - p.ramp_z)), z1)
    ramp = [(xs - 2, p.ramp_z - 2), (xs + p.tongue_len + 2, p.ramp_z + p.tongue_len + 2), (xs + p.tongue_len + 2, z1 + 1), (xs - 2, z1 + 1)]
    ramp = cq.Workplane("XZ").polyline(ramp).close().extrude(20, both=True).translate((0, yc, 0))
    return flanks.intersect(ramp)


def tongues(xs, p=P, clr=0.0, extra=0.0):
    """Both dovetails at the seam x = xs, top and bottom side, pointing +x.
    They belong to the section on the left, so at every seam that section goes
    in from the back of the one on its right."""
    top = tongue(xs, p, clr, extra)
    return top.union(top.mirror("XZ"))


def make_sections(p=P, body=None, support=False):
    """The frame cut into sections at the seams, left to right; with support,
    as printed, with the lip's support wall."""
    body = body if body is not None else make_body(p)
    if support:
        body = body.union(lip_support(p))
    xs = seams(p)
    w, h = outer(p)
    lo, hi = -w, w
    out = []
    edges = [lo] + xs + [hi]
    for k, (x0, x1) in enumerate(zip(edges, edges[1:])):
        region = box(x0, x1, -h, h, -1, p.depth + 1)
        if k > 0:  # sockets for the left neighbour's tongues
            region = region.cut(tongues(x0, p, p.joint_clr, 1.0))
        if k < len(edges) - 2:
            region = region.union(tongues(x1, p))
        out.append(body.intersect(region))
    return out


def section_names(p=P):
    n = sections(p)
    if n == 1:
        return ["frame"]
    if n == 2:
        return ["frame_left", "frame_right"]
    return [f"frame_section_{k + 1}" for k in range(n)]


def fit_test_parts(p=P, parts=None):
    """Small pieces cut from the real sections: the side either side of one
    seam, and the short side's lock with its keyhole, button peg and a length
    of lip (with its support wall, if the sections have it)."""
    parts = parts if parts is not None else make_sections(p)
    w, h = outer(p)
    sw, _ = sight(p)
    k = p.pad / 2
    out = {}
    if len(parts) > 1:
        xs = seams(p)[0]
        crop = box(xs - k, xs + k, rabbet(p)[1] / 2, h / 2 + 1, -1, p.depth + 1)  # no lip, so no support
        out["test_seam_left"] = parts[0].intersect(crop)
        out["test_seam_right"] = parts[1].intersect(crop)
    _, b = lock_span(p)  # the right side's pad runs towards -y
    out["test_lock"] = parts[-1].intersect(box(sw / 2 - 1, w / 2 + 1, -(b + 1), k, -1, p.depth + 1))
    return out


# --- Checks and export ------------------------------------------------------------


def bridge_zones(p=P):
    """Where a part printed face up may have flat undersides, in frame
    coordinates: the lip, which bridges to its support wall, and the pocket
    floors."""
    rw, rh = rabbet(p)
    sw, sh = sight(p)
    lip = box(-rw / 2, rw / 2, -rh / 2, rh / 2, p.skin - 0.1, p.skin + 0.1).cut(box(-sw / 2, sw / 2, -sh / 2, sh / 2, 0, p.depth))
    zf = pocket_floor(p)
    floors = [turn(cyl(pocket_r(p) + 2 * p.notch_r, zf - 0.1, zf + 0.1, c + btn_x(p), yr + p.btn_pivot), deg) for deg, yr, c in locks(p)]
    return [lip] + floors


def overhang_area(part, max_deg=48.0, allow=()):
    """Area (mm^2) of downward faces flatter than `max_deg` from vertical,
    ignoring the face on the bed and any inside the `allow` solids, which are
    planned bridges. The margin over 45 degrees absorbs meshing error on the
    45 degree faces."""
    verts, tris = part.val().tessellate(0.02, 0.1)
    zmin = min(v.z for v in verts)
    limit = -math.sin(math.radians(max_deg))
    total = 0.0
    for i, j, k in tris:
        a, b, c = verts[i], verts[j], verts[k]
        n = (b - a).cross(c - a)
        if n.Length < 1e-12:
            continue
        if n.z / n.Length < limit and (a.z + b.z + c.z) / 3 > zmin + 0.05:
            m = (a + b + c) * (1 / 3)
            if any(s.val().isInside(m, 0.01) for s in allow):
                continue
            total += n.Length / 2
    return total


def grams(part):
    """Solid PLA weight; the print is lighter where infill replaces solid."""
    return part.val().Volume() / 1000 * 1.24


def main(out=Path(__file__).parent / "out", p=P):
    out.mkdir(exist_ok=True)
    body = make_body(p)
    parts = dict(zip(section_names(p), make_sections(p, body)))
    printed = dict(zip(section_names(p), make_sections(p, body, support=True)))
    asm = cq.Assembly()
    for name, part in parts.items():
        asm.add(part, name=name, color=cq.Color("gray20"))
    asm.export(str(out / "assembly.step"))
    parts = {n: face_up(s, p) for n, s in printed.items()} | {n: face_up(t, p) for n, t in fit_test_parts(p, list(printed.values())).items()}
    parts["button"] = button_print(p)
    for name, part in parts.items():
        bb = part.val().BoundingBox()
        cq.exporters.export(part, str(out / f"{name}.step"))
        cq.exporters.export(part, str(out / f"{name}.stl"), tolerance=0.01, angularTolerance=0.1)
        print(f"{name:16s} {bb.xlen:6.1f} x {bb.ylen:6.1f} x {bb.zlen:5.1f} mm  {grams(part):5.1f} g")
    w, h = outer(p)
    print(f"frame {w:.1f} x {h:.1f} mm, opening {sight(p)[0]:.1f} x {sight(p)[1]:.1f} mm, {grams(body):.0f} g")


if __name__ == "__main__":
    main()

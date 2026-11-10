#!/usr/bin/env python3
"""Flightscope Mk I enclosure - parametric CadQuery model.

Landscape handheld. Front to back: bezel, display module, mid-frame (retains the display),
Raspberry Pi 4 (component side toward the display), UPS HAT (D) with two 21700 cells, rear cover.
The FlightAware Pro Stick Plus lies along the top edge in a clip cradle on the rear cover.

Parts (all print without supports on a P1S, PETG, 0.2 mm):
  front_shell.stl   bezel + walls, print face-down          mid_frame.stl     flat plate, print flat
  rear_cover.stl    print inner-face-up (bosses + cradle)

Run:  python flightscope_case.py           -> writes ./out/*.stl and ./out/*.svg previews
Every dimension is a named constant below. Ones marked MEASURE need calipers on the real part.
"""
import cadquery as cq
from pathlib import Path

OUT = Path(__file__).resolve().parent / "out"; OUT.mkdir(exist_ok=True)

# ----------------------------------------------------------------------------- components (mm)
# Waveshare 4.3" DSI LCD (SKU 16249). Glass and PCB are the same outline.
DISPLAY_W        = 106.0   # measured
DISPLAY_H        = 67.7
DISPLAY_GLASS_T  = 6.5     # whole LCD module in front of the PCB (glass + panel + backlight)
DISPLAY_PCB_T    = 1.6
DISPLAY_REAR_H   = 4.0     # corner standoff tops: the mid-frame bears here
DISPLAY_HOLE_IN  = 3.0     # corner standoff centres from the PCB corners (x and y)
ACTIVE_W, ACTIVE_H = 95.0, 54.0          # 4.3" 800x480 active area
ACTIVE_LEFT, ACTIVE_TOP = 5.0, 4.0       # glass edge to the lit picture (right border is 6, bottom 9.7: the flex tail end)
DSI_X, DSI_Y, DSI_LEN, DSI_H = 1.5, 24.0, 20.0, 6.0   # ribbon socket on the back, from the edge it sits against / from the bottom; it stands 2 mm above the standoffs

# Raspberry Pi 4B
PI_W, PI_H, PI_T = 85.0, 56.0, 1.6
PI_HOLES = [(3.5, 3.5), (61.5, 3.5), (3.5, 52.5), (61.5, 52.5)]   # from the SD-card end, GPIO edge at y=56
PI_COMP_H = 16.0           # USB/Ethernet stacks, on the top side (faces the display)
PI_AIR    = 2.6            # mid-frame back face -> those stacks. Not a choice: the stack hangs from the cover on the studs, so this is what is left

# Waveshare UPS HAT (D), from its STEP model. Same outline and hole pattern as the Pi.
HAT_T          = 1.6
HAT_CELL_DEPTH = 22.1      # cell crowns behind the HAT PCB's back face (STEP: 23.7 below its top face). The holder is an open tray: only its end walls come up that far
HAT_STANDOFF_END = 25.0    # the kit's 25 mm studs: their female ends, behind the back face, 2.9 mm beyond the cells
HAT_PI_GAP     = 9.0       # standoff height between HAT top and Pi PCB bottom (pogo pins)
HAT_USBC_X     = (44.8, 54.3)  # along the GPIO edge, from the SD-card end: the middle (STEP solid 9.6 x 7.5 x 4.2 at the edge). The 9 x 11.7 part at 9.5-18.5 that overhangs the edge is not the port
HAT_USBC_Z     = (1.5, 7.0)    # above the HAT PCB top face
HAT_BTN_X      = (70.4, 76.5)  # power button on the same edge
HAT_BTN_Z      = (0.0, 4.0)

# FlightAware Pro Stick Plus
DONGLE_L, DONGLE_W, DONGLE_T = 70.0, 30.4, 13.8  # measured: body only
SMA_HOLE_D = 6.6                          # SMA bulkhead: plain round hole (no flat on this one)
SMA_JACK_D = 6.3                          # the stick's own jack

# Panel-mount USB-C extension (charging). The HAT's own receptacle stays inside; a right-angle
# male plug goes into it. The flange sits on a pad inside the rear cover and is screwed from OUTSIDE. Set these from the one you buy.
PANEL_USBC_FLANGE = (24.0, 14.0)          # MEASURE  flange outline (BTFO panel extension, two diagonal ears)
PANEL_USBC_HOLES  = [(-8.75, 0.0), (8.75, 0.0)]  # MEASURE  two holes either side of the port on the centreline, ~17.5 mm apart (dx along the flange, dy across)
PANEL_USBC_PAD    = 4.5                   # local wall thickness under the flange: the panel its nose passes through
PANEL_USBC_BODY   = (13.0, 8.0)           # cutout for the connector body through the wall
PANEL_USBC_ROT    = 90                    # 0: flange along x (slot across the cover's width, centre 14 mm up); 90: flange along y (slot up and down, centre 19 mm up)
PANEL_USBC_XY     = None                  # set below once OUT_W is known
PLUG_ROOM         = 10.0                  # rigid right-angle USB-C adapter on the HAT's port, its socket facing sideways

# ----------------------------------------------------------------------------- case
WALL     = 2.4
BEZEL_T  = 2.5
CLR      = 0.5                            # general clearance around boards
R_OUT    = 6.0                            # outer corner radius
# screws: M3, except the stack's four, which are M2.5: they go into
# the HAT kit's M2.5 studs, whose male ends pass through the Pi/HAT's 2.7 mm holes. So do the mid-frame's optional screws into the display's own posts.
INSERT_D = 4.0                            # M3 brass heat-set insert, 4.2-4.6 mm OD (generic M3 x 4.2, Ruthex, CNC Kitchen): 4.0 mm hole in PETG
SCREW_D  = 3.4                            # M3 clearance
M25_D    = 2.8                            # M2.5 clearance: the stack studs, the display's posts
PILLAR_R = 3.8                            # bezel pillars: 1.8 mm of PETG round the insert, as the rim bosses
# finish
EDGE_R_FRONT = 2.0                        # round-over on the bezel's outer edge (prints face down: the first two layers step in, fine in PETG)
EDGE_R_BACK  = 1.5                        # round-over on the rear cover's outer edge
CSK_D        = 6.4                        # countersink for M3 countersunk screws (ISO 7046 / 14581: head 5.5, 6.3 to the cone's sharp edge), so every screw on the back sits flush or a hair under
M25_CSK_D    = 5.6                        # the same for the stack's M2.5s (head 4.7, 5.5 to the sharp edge)
WINDOW_CHAMFER = 1.5                      # 45-degree bevel round the screen window: a finger reaches the picture's edges (the buttons there) past the 2.5 mm bezel
LIP_T, LIP_H, LIP_CLR = 1.2, 2.0, 0.3   # locating lip on the cover's inner face : 1.2 thick, 2 tall, 0.3 inside the shell's walls
LANYARD_W, LANYARD_OUT, LANYARD_SLOT = 28.0, 15.0, (16.0, 5.0)   # one lanyard loop on the left side edge of the cover 
SEAM_V       = 0.4                        # 45-degree chamfer on both sides of the shell/cover joint: a small V line that hides any step between two separately printed parts
DISP_GUIDE_H, DISP_GUIDE_T, DISP_GUIDE_L = 4.0, 1.6, 8.0   # L-shaped corner guides on the bezel: the display drops in square to the window (CLR all round) and can't slide
CELL_KEEP_GAP = 1.0                       # cover ribs across the two 21700s stop this far short of them: a drop on the screen can't lift a cell off its end contacts
LOGO_DEPTH, LOGO_SCALE = 0.6, 1.4         # the app's own plane silhouette (ui/splash.py PLANE) with its dotted trail, recessed into the top wall
LOGO_PLANE = [(0, -11), (2, -7), (2, -2), (11, 3), (11, 5), (2, 3), (2, 7), (5, 10), (5, 11), (0, 10),
              (-5, 11), (-5, 10), (-2, 7), (-2, 3), (-11, 5), (-11, 3), (-2, -2), (-2, -7)]   # screen coords: nose at -y
GRIP_R, GRIP_PITCH, GRIP_N = 1.0, 8.0, 9   # grip ribs on both side walls : half-round, 1 mm proud, on the vent bars' 8 mm pitch so each rib runs between two vent slots
LOGO_NAME, LOGO_NAME_SIZE = "", 9.0     # a name in the bottom-right corner of the back, in the app's font (Nunito Black), e.g. "SAM"; "" for none
_HERE = Path(__file__).resolve().parent
LOGO_FONT = next((str(f) for f in (_HERE / "fonts/Nunito-Black-flat.ttf",) if f.exists()), None)   # Nunito Black with its overlapping outlines merged (fontTools removeOverlaps): the stock TTF makes an invalid solid that wrecks the cut

# internal envelope
# The Pi's USB-A port feeds the stick through the measured right-angle USB-A extension: elbow in the Pi's top USB2 port
# (lead exits toward the cover), right-angle female block on the stick's plug with its lead leaving toward the bezel and
# looping back into the Pi's port - nothing beyond the block along the stick (USB_LINK = "cable"). A flat FPC ribbon (USB_LINK = "ribbon") makes the box ~12 mm narrower.
USB_ELBOW      = 22.0     # measured: right-angle USB-A cable, how far it reaches beyond the Pi's PCB edge
USB_ELBOW_W    = 17.0     # along the port row
USB_PORT_Y     = 18.5     # the elbow's low edge in the Pi's frame: in the top port of the middle (blue) USB3 pair. In the black USB2 pair (3.0), nearer the port edge, it sits right under the power switch once the cover is turned. The stick runs at USB2 speed in either
USB_FEMALE     = 9.9      # measured: the cable's right-angle female block beyond the stick's plug tip
USB_FEMALE_H   = 20.0     # MEASURE  the block: tall (across the stick's width; the lead leaves its bezel-side face)
USB_FEMALE_T   = 14.0     # MEASURE  and deep (across the stick's thickness); it covers the last 8 mm of the plug
USB_FPC_PLUG   = 6.0      # ribbon's bare A plug: beyond the Pi's port face before the ribbon can bend (MEASURE)
USB_SOCKET_L   = 14.5     # ribbon's A socket: length along the plug axis (MEASURE)
USB_SOCKET_W   = 13.3     # across the socket mouth (MEASURE)
USB_SOCKET_T   = 5.8      # socket height, ribbon included (MEASURE)
DONGLE_PLUG    = 13.0     # measured: the stick's USB-A plug beyond its body
SMA_JACK       = 10.0     # the stick's SMA jack beyond the other end
ANT_PLUG       = "bulkhead" # where the antenna lead's straight SMA plug screws onto the stick's jack:
                          #   "bulkhead" - an SMA female bulkhead in the TOP wall behind the power switch (round hole, nut
                          #               inside), antenna screwed on outside. A right-angle SMA adapter on the stick's jack turns the
                          #               pigtail's plug toward the screen; the pigtail runs along the band in front of the stick to it
                          #   "inside"  - inside the box; only the lead's neck passes out through the left wall, so the stick sits
                          #               ANT_PLUG_L in from the wall (the cradle sits mid-band rather than hard against one end)
                          #   "outside" - the stick is hard against the wall, its jack passes through, the plug screws on outside
ANT_PLUG_L     = 12.0     # MEASURE  inside variant: stick end face -> the lead's neck (knurled nut + body over the 10 mm jack)
ANT_LEAD_D     = 5.0      # MEASURE  the lead's neck/boot where it passes the wall (the wall slot is SMA_HOLE_D wide)
ANT_RA_X       = 15.0     # MEASURE  bulkhead: stick end face -> the right-angle adapter's leg axis (screwed on)
ANT_RA_LEG     = 16.0     # MEASURE  bulkhead: jack axis -> the adapter's female face (the leg points toward the screen)
ANT_PIG_PLUG   = 20.0     # MEASURE  bulkhead: the pigtail's straight SMA plug + boot beyond that face
SMA_BH_IN      = 12.0     # MEASURE  bulkhead: thread + nut + crimp inside the top wall before its cable can bend
SMA_END_MIN    = {"outside": 1.0, "bulkhead": ANT_RA_X + 4.5 + 0.5}.get(ANT_PLUG, 0.5 + ANT_PLUG_L + 0.5)   # least room between the left wall and the stick end face: the jack, the right-angle adapter's body (4.5 beyond its leg axis), or the antenna plug
POCKET_L       = 1.0 + USB_SOCKET_L + 0.8 + 2.0   # socket pocket at the stick's plug end: 1 mm to the socket mouth, socket, play, end wall
USB_LINK       = "cable"  # "cable" (the measured right-angle extension, default) or "ribbon" (flat FPC extension: 12 mm narrower box)
USB_CABLE_ROOM = 1.0      # cable variant: air beyond the female block (its lead leaves at a right angle, so no bend room along the stick)
PLUG_END_L     = POCKET_L if USB_LINK == "ribbon" else DONGLE_PLUG + USB_FEMALE + USB_CABLE_ROOM
BOSS_R         = 3.8      # cover screw bosses: 1.8 mm of PETG around the 4.0 mm M3 insert hole; moulded into two walls at each corner
GPS_BACK_T     = 3.5      # GPS board: CR1220 holder + solder stubs on the back
GPS_PATCH_T    = 4.5      # its 15x15x4 ceramic patch on the front
GPS_SLOT       = 0.5 + GPS_PATCH_T + 1.6 + GPS_BACK_T + 1.5   # the GPS stands on edge in the band between the bottom wall and the stack
RIBBON_ROOM    = 7.0      # right of the stack (SD end): the DSI ribbon's loop out of the Pi's DISPLAY socket and the SD card
IN_W_PRINTED   = 120.9    # minimum inside width (box 125.7 wide)
IN_W = max(DISPLAY_W + 2 * CLR, RIBBON_ROOM + PI_W + (2 + USB_FPC_PLUG if USB_LINK == "ribbon" else USB_ELBOW) + 3,   # ribbon at the SD end; ports + plug at the other
           SMA_END_MIN + DONGLE_L + PLUG_END_L + 0.5 + (2 * BOSS_R - 0.3 if USB_LINK == "ribbon" else 0), IN_W_PRINTED)   # top band at its tightest
SMA_END_L = SMA_END_MIN if ANT_PLUG == "outside" else max(SMA_END_MIN, min((IN_W - DONGLE_L) / 2, IN_W - DONGLE_L - PLUG_END_L - 0.5))   # jack through the wall: stick against it; otherwise centred in the box, as far as the plug end allows
IN_H = (7 + PLUG_ROOM) + PI_H + 4 + DONGLE_T + 3.0                  # plug room, stack, gap, dongle band, 0.5 air + 2 mm outer rail + 0.5 to the wall
IN_D = (BEZEL_T + DISPLAY_GLASS_T + DISPLAY_PCB_T + DISPLAY_REAR_H + 2.0    # mid-frame
        + PI_AIR + PI_COMP_H + PI_T + HAT_PI_GAP + HAT_T + HAT_STANDOFF_END + 0.6)
OUT_W, OUT_H, OUT_D = IN_W + 2 * WALL, IN_H + 2 * WALL, IN_D + WALL

# z positions measured from the FRONT face (z grows toward the back)
Z_GLASS   = BEZEL_T
Z_PCB     = Z_GLASS + DISPLAY_GLASS_T
Z_MID     = Z_PCB + DISPLAY_PCB_T + DISPLAY_REAR_H          # mid-frame front face
MID_T     = 2.0
Z_PI_TOP  = Z_MID + MID_T + PI_AIR + PI_COMP_H                  # Pi PCB front face (components point forward)
Z_HAT_TOP = Z_PI_TOP + PI_T + HAT_PI_GAP
Z_HAT_BOT = Z_HAT_TOP + HAT_T
Z_CELLS   = Z_HAT_BOT + HAT_CELL_DEPTH
Z_STANDOFF = Z_HAT_BOT + HAT_STANDOFF_END                    # rear cover bosses meet the stack here
Z_COVER   = OUT_D - WALL                                     # rear cover inner face

# xy: origin at the device's front-left-bottom outer corner; internal cavity starts at (WALL, WALL)
DISP_X = WALL + (IN_W - DISPLAY_W) / 2
DISP_Y = WALL + 7.0 + PLUG_ROOM                              # display above the plug zone, dongle above the stack
PILLARS = [(WALL + 3.5, WALL + 3.5), (OUT_W - WALL - 3.5, WALL + 3.5), (WALL + 3.5, OUT_H - WALL - 3.5), (OUT_W - WALL - 3.5, OUT_H - WALL - 3.5)]
# cover screw bosses: one in each corner. Short (10 mm), gussets receding into the corner beneath. Both top ones
# sit beside parts whose size is a guess (the female block's section, the antenna adapter): fit_test's usb_corner print tries the USB one
RIM_BOSSES = [(x, y) for y in (WALL + BOSS_R - 0.3, OUT_H - WALL - BOSS_R + 0.3) for x in (WALL + BOSS_R - 0.3, OUT_W - WALL - BOSS_R + 0.3)]
BOSS_H, GUSSET_H = 10.0, 9.0                                    # boss height; 45-degree gusset tapering into the corner beneath it
# The stack, in the cover's frame (see COVER_ROT): HAT ports DOWN into the plug room, the Pi's USB/Ethernet facing the LEFT wall, the SD end at
# the right. The Pi's own port edge (USB-C power, micro-HDMI, audio) is the bottom one, next to the HAT's ports, and its GPIO header runs along
# the top (the USB ports need the open space for the right-angle plug).
# pi_x/pi_y map a point in the Pi's own frame (x from the SD-card end, y from the edge opposite the GPIO) into the cover's frame.
STACK_X = OUT_W - WALL - RIBBON_ROOM - PI_W                  # ribbon room at the SD end (right); USB elbow room at the ports end (left)
STACK_Y = WALL + 7.0 + PLUG_ROOM                             # HAT ports face down; room below for the right-angle charging plug
def pi_x(x, w=0.0): return STACK_X + PI_W - x - w             # (w: a feature's extent along x, so a box's low corner lands right)
def pi_y(y, h=0.0): return STACK_Y + y
# The cover goes on turned 180 degrees in its own plane: everything on it (the stack, the stick in its cradle, the GPS, the
# charge socket) turns with it. The GPS patch then looks out through the TOP wall, at the sky, and the GPIO header runs along the bottom, away from
# the power switch. The four corner screws line up either way. Everything above is in the cover's own frame; on_cover() turns a solid from that
# frame into the box, cover_pt() a point. The front shell and what it carries (switch, SMA bulkhead, display, mid-frame) are in the box frame.
COVER_ROT = 180
def cover_pt(x, y): return (OUT_W - x, OUT_H - y) if COVER_ROT else (x, y)
def on_cover(w): return w.rotate((OUT_W / 2, OUT_H / 2, 0), (OUT_W / 2, OUT_H / 2, 1), COVER_ROT) if COVER_ROT else w
DONGLE_Y = STACK_Y + PI_H + 4.0                              # 2 mm rib + 2 mm air above the stack
DONGLE_X = WALL + SMA_END_L                                  # stick's SMA end: antenna plug between it and the left wall; USB plug + female end to the right
BTN_HEAD_OD, BTN_BODY_D, BTN_BODY_L, BTN_DEPTH = 18.0, 14.0, 17.0, 24.0   # MEASURE  switch: head ring outside, body behind the panel, body length, depth to the terminal tips
BTN_D, BTN_X, BTN_Z = 12.2, OUT_W - R_OUT - BTN_HEAD_OD / 2 - 2.0, 25.0   # latching 12 mm panel switch (two pins, no LED: across the HAT's slide switch, see case/README), right end of the top wall: its head ring sits on the flat, 2.0 mm short of the corner's round-over 
SMA_Y, SMA_Z = DONGLE_Y + DONGLE_T / 2, Z_COVER - DONGLE_W / 2   # the stick's jack axis: slot in the LEFT wall (open to the back edge) at the cradle's centre line ("inside"/"outside")
SMA_TOP_X, SMA_TOP_Z = BTN_X, BTN_Z + BTN_HEAD_OD / 2 + 3.0 + 6.0   # "bulkhead": top wall, right behind the switch, toward the cover: 3 mm between the switch head and the antenna's ~12 mm base
if PANEL_USBC_ROT == 90:                                     # as measured the flange is long along x; turn it to run along y
    PANEL_USBC_FLANGE = (PANEL_USBC_FLANGE[1], PANEL_USBC_FLANGE[0]); PANEL_USBC_BODY = (PANEL_USBC_BODY[1], PANEL_USBC_BODY[0])
    PANEL_USBC_HOLES = [(-dy, dx) for dx, dy in PANEL_USBC_HOLES]
    # rear cover, bottom-left corner, in the GPS's band at the other end of it: 0.5 mm
    # off the corner boss and the bottom wall, standing up the strip beside the Pi's ports (which are far below it). The charge lead runs along the bottom band to the adapter.
    PANEL_USBC_XY = (WALL + 2 * BOSS_R - 0.3 + 0.5 + (PANEL_USBC_FLANGE[0] + 8) / 2, WALL + 0.5 + (PANEL_USBC_FLANGE[1] + 8) / 2)
else:
    PANEL_USBC_XY = (WALL + 2 * BOSS_R + 1.0 + (PANEL_USBC_FLANGE[0] + 8) / 2, WALL + 4.0 + (PANEL_USBC_FLANGE[1] + 8) / 2)   # rear cover, bottom-left corner, below the stack

# Adafruit Ultimate GPS Breakout (PA1616D / v3): stands UPRIGHT on its short edge in the bottom-right
# corner, opposite the charge socket, long side along z, patch out through the bottom wall, header end up (toward the
# screen) so the leads go straight to the GPIO header on the Pi's bottom edge. The HAT's charge port and its adapter are
# in the middle of the band, clear of it. Both long edges slide down grooves in two full-height posts.
GPS_W, GPS_H, GPS_T = 25.5, 35.0, 1.6                          # PCB: short side (across, along x), long side (stands up, along z), thickness
GPS_Y     = WALL + 2.0                                          # patch face, 2.0 off the bottom wall (room for the patch module and the header leads)
GPS_PCB_Y = GPS_Y + GPS_PATCH_T                                 # PCB front face (patch side)
GPS_POST  = 1.5 + 0.3 + 2.0                                     # groove posts: 1.5 over the board's edge strip, 0.3 play, 2 mm of post beyond
GPS_X     = OUT_W - WALL - 2 * BOSS_R + 0.3 - 0.5 - (GPS_POST - 1.5) - GPS_W   # board's left edge: right post 0.5 mm clear of the corner boss
GPS_Z_BOT = Z_COVER - 0.3                                       # the short edge it stands on; header end at GPS_Z_BOT - GPS_H
GPS_HOLES = [(hx, hz) for hx in (3.0, GPS_W - 3.0) for hz in (3.0, GPS_H - 3.0)]   # not used: the board slides in   # one in each corner

print(f"outer {OUT_W:.1f} x {OUT_H:.1f} x {OUT_D:.1f} mm   (internal {IN_W:.1f} x {IN_H:.1f} x {IN_D:.1f})")

def rbox(w, h, d, r=R_OUT):
    return cq.Workplane("XY").box(w, h, d, centered=False).edges("|Z").fillet(r)

def stack_holes(z0, z1, d):
    """Cylinders at the Pi/HAT hole pattern, in case coordinates."""
    res = None
    for hx, hy in PI_HOLES:
        c = cq.Workplane("XY").workplane(offset=z0).center(pi_x(hx), pi_y(hy)).circle(d / 2).extrude(z1 - z0)
        res = c if res is None else res.union(c)
    return res

def rim_boss(bx, by):
    """A screw boss moulded INTO the wall(s): cylinder + full-height web to the side wall (and the end wall in a corner)
    + a 45-degree gusset below it that recedes into the wall, so on a face-down print it grows out of solid wall."""
    z0 = Z_COVER - BOSS_H
    left = bx < OUT_W / 2; top = by > OUT_H / 2
    xw = WALL - 0.5 if left else OUT_W - WALL + 0.5                      # wall faces (with a little embed)
    yw = OUT_H - WALL + 0.5 if top else WALL - 0.5
    b = cq.Workplane("XY").workplane(offset=z0).center(bx, by).circle(BOSS_R).extrude(BOSS_H)
    x0, x1 = (xw, bx) if left else (bx, xw)                              # web to the side wall
    b = b.union(cq.Workplane("XY").box(x1 - x0, 2 * BOSS_R, BOSS_H, centered=False).translate((x0, by - BOSS_R, z0)))
    if WALL + BOSS_R + 1 < by < OUT_H - WALL - BOSS_R - 1:               # side-wall boss: gusset from a sliver on the wall up to the footprint
        sx0, sx1 = (xw, xw + 0.8) if left else (xw - 0.8, xw)
        fx0, fx1 = (xw, bx + BOSS_R) if left else (bx - BOSS_R, xw)
        g = (cq.Workplane("XY").workplane(offset=z0 - GUSSET_H).center((sx0 + sx1) / 2, by).rect(sx1 - sx0, 2 * BOSS_R)
             .workplane(offset=GUSSET_H).center((fx0 + fx1) / 2 - (sx0 + sx1) / 2, 0).rect(fx1 - fx0, 2 * BOSS_R).loft())
        return b.union(g)
    y0, y1 = (by, yw) if top else (yw, by)                               # web to the end wall
    b = b.union(cq.Workplane("XY").box(2 * BOSS_R, y1 - y0, BOSS_H, centered=False).translate((bx - BOSS_R, y0, z0)))
    fx0, fx1 = (xw, bx + BOSS_R) if left else (bx - BOSS_R, xw)          # gusset: loft from a sliver in the corner up to the boss footprint
    fy0, fy1 = (by - BOSS_R, yw) if top else (yw, by + BOSS_R)
    sx0, sx1 = (xw, xw + 0.8) if left else (xw - 0.8, xw)
    sy0, sy1 = (yw - 0.8, yw) if top else (yw, yw + 0.8)
    g = (cq.Workplane("XY").workplane(offset=z0 - GUSSET_H).center((sx0 + sx1) / 2, (sy0 + sy1) / 2).rect(sx1 - sx0, sy1 - sy0)
         .workplane(offset=GUSSET_H).center((fx0 + fx1) / 2 - (sx0 + sx1) / 2, (fy0 + fy1) / 2 - (sy0 + sy1) / 2).rect(fx1 - fx0, fy1 - fy0).loft())
    return b.union(g)

def box2(x0, y0, x1, y1, z0, z1):
    """Box between two corners in any order."""
    return cq.Workplane("XY").box(abs(x1 - x0), abs(y1 - y0), z1 - z0, centered=False).translate((min(x0, x1), min(y0, y1), z0))

def countersink(x, y, d=SCREW_D, csk=CSK_D):
    """90-degree cone in the cover's outer face for a countersunk head, M3 unless told otherwise (prints fine face down: a 45-degree overhang)."""
    h = (csk - d) / 2
    return cq.Workplane("XY").add(cq.Solid.makeCone(d / 2, csk / 2 + 0.5, h + 0.5, cq.Vector(x, y, OUT_D - h), cq.Vector(0, 0, 1)))

def top_logo():
    """The app's plane with its dotted trail, recessed into the TOP wall, flying in line with the antenna and the power
    switch. Symmetric about its own axis."""
    sc = LOGO_SCALE; cz = BTN_Z                                          # on the switch's centre line
    cx = min(SMA_TOP_X - 6.0, BTN_X - BTN_HEAD_OD / 2) - 8.0 - 11 * sc   # nose 8 mm short of the antenna's ~12 mm base or the switch's head ring
    pts = [(cx - py * sc, cz + px * sc) for px, py in LOGO_PLANE]        # nose (-py) toward +x, wings along z
    wp = lambda: cq.Workplane("XZ").workplane(offset=-(OUT_H + 1))       # XZ normal is -Y: start outside the top face, cut inward
    g = wp().polyline(pts).close().extrude(1 + LOGO_DEPTH)
    tx = cx - 11 * sc - 2.5
    for i, r in enumerate((1.3, 1.15, 1.0, 0.85, 0.7)):                 # trail: shrinking dots, as on the splash screen
        g = g.union(wp().center(tx - i * 4.2, cz).circle(r).extrude(1 + LOGO_DEPTH))
    return g

def back_name():
    """LOGO_NAME in the bottom-right corner of the back (as seen from behind), in the app's font. XY text reads correctly
    from +Z, i.e. from behind. In the box frame, wherever the cover's own parts end up: 3 mm in from the corner boss's footprint, clear of the
    corner screw and of the stack's countersinks."""
    kw = {"fontPath": LOGO_FONT} if LOGO_FONT else {"kind": "bold"}
    hx = OUT_W - WALL - 2 * BOSS_R - 3.0                                  # right-aligned here
    return (cq.Workplane("XY").workplane(offset=OUT_D - LOGO_DEPTH).center(hx, WALL + BOSS_R + 5.5)
            .text(LOGO_NAME, LOGO_NAME_SIZE, LOGO_DEPTH + 1, combine=False, halign="right", valign="center", **kw))

# ----------------------------------------------------------------------------- front shell
def front_shell():
    s = (rbox(OUT_W, OUT_H, Z_COVER).faces("<Z").edges().fillet(EDGE_R_FRONT)       # walls run from the bezel to the cover plane; rounded front edge
         .faces(">Z").edges().chamfer(SEAM_V))                                        # half of the V at the joint with the cover
    cavity = rbox(IN_W, IN_H, Z_COVER, R_OUT - WALL).translate((WALL, WALL, BEZEL_T))
    s = s.cut(cavity)
    # display corner guides: locate the display sideways, centred in the window
    t, L, z0, z1 = DISP_GUIDE_T, DISP_GUIDE_L, BEZEL_T - 0.1, BEZEL_T + DISP_GUIDE_H
    for cx, sx in ((DISP_X - CLR, -1), (DISP_X + DISPLAY_W + CLR, 1)):
        for cy, sy in ((DISP_Y - CLR, -1), (DISP_Y + DISPLAY_H + CLR, 1)):
            s = s.union(box2(cx - sx * L, cy, cx + sx * t, cy + sy * t, z0, z1)).union(box2(cx, cy - sy * L, cx + sx * t, cy + sy * t, z0, z1))
    # bezel window
    wx = DISP_X + ACTIVE_LEFT - 1.0
    wy = DISP_Y + (DISPLAY_H - ACTIVE_TOP - ACTIVE_H) - 1.0
    s = s.cut(cq.Workplane("XY").box(ACTIVE_W + 2, ACTIVE_H + 2, BEZEL_T + 2, centered=False).translate((wx, wy, -1)))
    cw = WINDOW_CHAMFER                                                  # 45-degree chamfer round the window's front edge
    s = s.cut(cq.Workplane("XY").workplane(offset=-0.01).center(wx + (ACTIVE_W + 2) / 2, wy + (ACTIVE_H + 2) / 2)
              .rect(ACTIVE_W + 2 + 2 * cw + 0.02, ACTIVE_H + 2 + 2 * cw + 0.02).workplane(offset=cw + 0.01).rect(ACTIVE_W + 2, ACTIVE_H + 2).loft())
    # four corner pillars from the bezel plane: the mid-frame screws to these and clamps the display against the bezel lip
    for (px, py) in PILLARS:
        s = s.union(cq.Workplane("XY").workplane(offset=BEZEL_T - 0.1).center(px, py).circle(PILLAR_R).extrude(Z_MID - BEZEL_T + 0.1))
        s = s.cut(cq.Workplane("XY").workplane(offset=Z_MID - 7).center(px, py).circle(INSERT_D / 2).extrude(8))
    # rim bosses with heat-set inserts for the rear cover (4), one in each corner
    for (bx, by) in RIM_BOSSES:
        s = s.union(rim_boss(bx, by))
        s = s.cut(cq.Workplane("XY").workplane(offset=Z_COVER - 8).center(bx, by).circle(INSERT_D / 2).extrude(9))
    # the antenna lead (or, ANT_PLUG = "outside", the stick's own jack) passes out through the left wall: a round-bottomed
    # slot open to the back edge, so it drops in as the cover goes on; the cover plate closes the slot
    if ANT_PLUG == "bulkhead":                                   # bulkhead: closed round hole in the top wall beside the switch
        s = s.cut(cq.Workplane("XZ").workplane(offset=-(OUT_H + 1)).center(SMA_TOP_X, SMA_TOP_Z).circle(SMA_HOLE_D / 2).extrude(WALL + 2))
    else:                                                        # lead (or jack) drops into a slot open to the back edge
        s = s.cut(cq.Workplane("YZ").workplane(offset=-1).center(SMA_Y, SMA_Z).circle(SMA_HOLE_D / 2).extrude(WALL + 2))
        s = s.cut(cq.Workplane("XY").box(WALL + 2, SMA_HOLE_D, Z_COVER - SMA_Z + 1, centered=False).translate((-1, SMA_Y - SMA_HOLE_D / 2, SMA_Z)))
    # power switch beside it: a latching 12 mm panel button, wired to the HAT's slide-switch pads
    s = s.cut(cq.Workplane("XZ").workplane(offset=-(OUT_H + 1)).center(BTN_X, BTN_Z).circle(BTN_D / 2).extrude(WALL + 2))   # plain round hole
    # side vents at the Pi's depth
    for side_x in (WALL / 2, OUT_W - WALL / 2):
        for i in range(6):
            vy = WALL + 12 + i * 8
            s = s.cut(cq.Workplane("YZ").workplane(offset=side_x - WALL / 2 - 1).center(vy + 2.25, Z_MID + 4 + 7)
                      .slot2D(14, 4.5, angle=90).extrude(WALL + 2))           # round-ended slots
    s = s.cut(top_logo())
    # grip ribs: horizontal (front-to-back) half-round ribs on both side walls, where the hands go, so it can't slide
    # through small fingers. Centred on the bars between the vent slots (vy + 4.5 + 1.75), 45-degree cone ends (no overhang when
    # printed face down). They make the shell 2 mm wider overall; the cover and the box inside are unchanged.
    z0, z1 = EDGE_R_FRONT + 4.0, Z_COVER - 3.0
    for i in range(GRIP_N):
        gy = WALL + 12 + 4.5 + 1.75 - GRIP_PITCH + i * GRIP_PITCH
        for gx in (0.0, OUT_W):
            rib = (cq.Solid.makeCylinder(GRIP_R, z1 - z0, cq.Vector(gx, gy, z0), cq.Vector(0, 0, 1))
                   .fuse(cq.Solid.makeCone(0, GRIP_R, GRIP_R, cq.Vector(gx, gy, z0 - GRIP_R), cq.Vector(0, 0, 1)))
                   .fuse(cq.Solid.makeCone(GRIP_R, 0, GRIP_R, cq.Vector(gx, gy, z1), cq.Vector(0, 0, 1))))
            s = s.union(cq.Workplane("XY").add(rib))
    return s

# ----------------------------------------------------------------------------- mid-frame
def mid_frame():
    """Ring that clamps the display: screwed to the four bezel pillars. Side rails are opened up so the
    DSI ribbon and the display's own connector are free; the display's corners sit under the ring."""
    m = rbox(IN_W - 0.6, IN_H - 0.6, MID_T, R_OUT - WALL).translate((WALL + 0.3, WALL + 0.3, Z_MID))
    inner = cq.Workplane("XY").box(DISPLAY_W - 4, DISPLAY_H - 2 * (DISPLAY_HOLE_IN + 5), MID_T + 2, centered=False)
    m = m.cut(inner.translate((DISP_X + 2, DISP_Y + DISPLAY_HOLE_IN + 5, Z_MID - 1)))
    # the DSI socket (6 mm) stands 2 mm above the standoffs the frame bears on: open both side rails over it,
    # so it clears whichever way round the display goes in
    for sx in (DISP_X + DSI_X - 2.0, DISP_X + DISPLAY_W - DSI_X - DSI_LEN - 2.0):
        m = m.cut(cq.Workplane("XY").box(DSI_LEN + 4.0, DSI_LEN + 6.0, MID_T + 2, centered=False).translate((sx, DISP_Y + DSI_Y - 3.0, Z_MID - 1)))
    for (px, py) in PILLARS:
        m = m.cut(cq.Workplane("XY").workplane(offset=Z_MID - 1).center(px, py).circle(SCREW_D / 2).extrude(MID_T + 2))
    for hx in (DISP_X + DISPLAY_HOLE_IN, DISP_X + DISPLAY_W - DISPLAY_HOLE_IN):          # optional: into the display's own standoffs
        for hy in (DISP_Y + DISPLAY_HOLE_IN, DISP_Y + DISPLAY_H - DISPLAY_HOLE_IN):
            m = m.cut(cq.Workplane("XY").workplane(offset=Z_MID - 1).center(hx, hy).circle(M25_D / 2).extrude(MID_T + 2))
    return m

# ----------------------------------------------------------------------------- rear cover
def rear_cover():
    c = rbox(OUT_W, OUT_H, WALL).faces(">Z").edges().fillet(EDGE_R_BACK).faces("<Z").edges().chamfer(SEAM_V).translate((0, 0, Z_COVER))   # rounded back edge; the other half of the V at the joint
    # screw holes matching the shell's rim bosses, countersunk from the back (flush heads)
    for (bx, by) in RIM_BOSSES:
        c = c.cut(cq.Workplane("XY").workplane(offset=Z_COVER - 1).center(bx, by).circle(SCREW_D / 2).extrude(WALL + 2))
        c = c.cut(countersink(bx, by))
    # bosses that meet the stack's standoff ends; M2.5 screws from the back into the HAT standoffs
    c = c.union(stack_holes(Z_STANDOFF, Z_COVER, 7.0))                 # short bosses that meet the standoff ends
    c = c.cut(stack_holes(Z_STANDOFF - 1, OUT_D + 1, M25_D))
    for hx, hy in PI_HOLES: c = c.cut(countersink(pi_x(hx), pi_y(hy), M25_D, M25_CSK_D))
    # cell keepers: two 2 mm ribs across both cells, CELL_KEEP_GAP short of their crowns. The HAT's holder is an open tray that holds
    # each cell only by its end contacts, so without them a drop on the screen could lift a cell 3.5 mm, off the contacts. Along the
    # cells at a quarter of the way in from each end: clear of the holder's end walls and the stack bosses
    for kx in (24.5, 60.5):
        c = c.union(box2(pi_x(kx) - 1.0, pi_y(10.0), pi_x(kx) + 1.0, pi_y(46.0), Z_CELLS + CELL_KEEP_GAP, Z_COVER + 0.1))
    # dongle cradle: the stick stands on its long edge along the top wall. A solid outer rail (backed by the top wall once
    # the cover is on), three snap hooks on the inner side that latch over its front edge, a low stiffening rail between
    # them, and an end saddle at each end so it cannot slide along - the SMA jack and the USB plug pass through notches.
    zd = Z_COVER - DONGLE_W / 2                                          # mid-height of the stick: SMA jack / plug axis
    ry = DONGLE_Y + DONGLE_T + 0.5                                       # outer rail, 0.5 mm off the stick
    rx0 = DONGLE_X - 0.3 - 3.0 if ANT_PLUG != "outside" else R_OUT       # rails start at the SMA end stop (or at the wall fillet when the stick is against the wall)
    c = c.union(cq.Workplane("XY").box(DONGLE_X + DONGLE_L + 3.5 - rx0, 2.0, 12, centered=False).translate((rx0, ry, Z_COVER - 12)))
    c = c.union(cq.Workplane("XY").box(DONGLE_X + DONGLE_L + 3.5 - rx0, 2.4, 6, centered=False).translate((rx0, DONGLE_Y - 2.4, Z_COVER - 6)))
    hl = DONGLE_W + 0.3                                                  # lip underside: 0.3 mm in front of the stick's edge
    hook = (cq.Workplane("YZ").polyline([(0, 0), (2.4, 0), (2.4, -hl), (3.6, -hl), (3.6, -hl - 0.7), (2.4, -hl - 1.5), (0, -hl - 1.5)]).close()
            .extrude(12).translate((0, DONGLE_Y - 2.4, Z_COVER)))      # 2.4 mm fin, 1.2 mm lip, entry chamfer
    for hx in (DONGLE_X + 14, DONGLE_X + DONGLE_L / 2, DONGLE_X + DONGLE_L - 14):
        c = c.union(hook.translate((hx - 6, 0, 0)))
    px0 = DONGLE_X + DONGLE_L + 0.3; cy = DONGLE_Y + DONGLE_T / 2
    if ANT_PLUG != "outside":
        # SMA end: two 3 mm pillars the stick's end face bears on, either side of a 9.5 mm gap the jack and the antenna
        # plug's knurled nut sit in (the nut comes right up to the stick's face, so nothing may stand in front of the jack)
        sx0 = DONGLE_X - 0.3 - 3.0
        stop = cq.Workplane("XY").box(3.0, DONGLE_T + 4.9, DONGLE_W / 2, centered=False).translate((sx0, DONGLE_Y - 2.4, zd))
        stop = stop.cut(cq.Workplane("XY").box(5.0, DONGLE_T - 2.4, DONGLE_W / 2 + 6, centered=False).translate((sx0 - 1, DONGLE_Y + 1.2, zd - 5)))   # same gap as the USB end
        c = c.union(stop)
    # (ANT_PLUG = "outside": the stick's end face bears on the shell's left wall, its jack goes out through the slot - no stop)
    if USB_LINK == "ribbon":
        # plug end: a U pocket (open toward the front) the ribbon's USB-A socket drops into; the stick's plug then mates into
        # it and the block is the stick's other end stop. End wall slotted along both faces so the ribbon can leave either way.
        ct = USB_SOCKET_T + 0.8; cw = USB_SOCKET_W + 1.1
        blk = cq.Workplane("XY").box(POCKET_L, DONGLE_T + 4.9, 24.0, centered=False).translate((px0, DONGLE_Y - 2.4, Z_COVER - 24.0))
        blk = blk.cut(cq.Workplane("XY").box(1.0 + USB_SOCKET_L + 0.8 + 1.0, ct, cw + 4, centered=False).translate((px0 - 1.0, cy - ct / 2, zd - cw / 2 - 4)))
        for sy in (cy - ct / 2, cy + ct / 2 - 1.6):                      # ribbon slots through the end wall
            blk = blk.cut(cq.Workplane("XY").box(4.0, 1.6, cw + 4, centered=False).translate((px0 + POCKET_L - 3.0, sy, zd - cw / 2 - 4)))
    else:
        # plug end, cable variant: a 3 mm stop with a slot the stick's USB-A plug tongue passes through; the cable's
        # right-angle female block sits beyond it, its lead leaving toward the bezel
        blk = cq.Workplane("XY").box(3.0, DONGLE_T + 4.9, DONGLE_W / 2, centered=False).translate((px0, DONGLE_Y - 2.4, zd))
        # gap nearly the stick's full thickness (the plug is off-centre toward one face), cut down to the cover plate, so this
        # end is two pillars exactly like the SMA end
        blk = blk.cut(cq.Workplane("XY").box(5, DONGLE_T - 2.4, DONGLE_W / 2 + 6, centered=False).translate((px0 - 1, DONGLE_Y + 1.2, zd - 5)))   # through to the cover plate: two pillars, identical to the SMA end
    c = c.union(blk)
    # panel-mount USB-C charging socket: flange on a pad inside the cover, snug body cutout, and two countersunk clearance holes so
    # its screws go in from the back face into the flange's ears
    fw, fh = PANEL_USBC_FLANGE; bw, bh = PANEL_USBC_BODY; px, py = PANEL_USBC_XY
    c = c.union(cq.Workplane("XY").box(fw + 8, fh + 8, PANEL_USBC_PAD - WALL, centered=False).translate((px - fw / 2 - 4, py - fh / 2 - 4, Z_COVER - (PANEL_USBC_PAD - WALL))))
    c = c.cut(cq.Workplane("XY").box(bw, bh, PANEL_USBC_PAD + 2, centered=False).translate((px - bw / 2, py - bh / 2, Z_COVER - (PANEL_USBC_PAD - WALL) - 1)))
    for dx, dz in PANEL_USBC_HOLES:
        c = c.cut(cq.Workplane("XY").workplane(offset=Z_COVER - (PANEL_USBC_PAD - WALL) - 1).center(px + dx, py + dz).circle(SCREW_D / 2).extrude(PANEL_USBC_PAD + 2))
        c = c.cut(countersink(px + dx, py + dz))
    # GPS stand: the board stands upright on its short edge, patch out through the bottom wall. Its two long edges slide
    # down grooves (1.5 mm over the edge strip on the back, 0.5 on the component face, 0.2 mm play each side) in two full-height posts standing
    # on the cover plate (from 0.3 off the shell's bottom wall). Nothing touches the component faces; the top (header) end is open for the leads.
    gpy = GPS_PCB_Y; zt = Z_COVER - GPS_H - 0.5
    for ex, side in ((GPS_X, -1), (GPS_X + GPS_W, 1)):
        x0 = ex - 1.5 if side > 0 else ex + 1.5 - GPS_POST
        post = cq.Workplane("XY").box(GPS_POST, gpy + GPS_T + 0.2 + 1.6 - (WALL + 0.3), Z_COVER - zt, centered=False).translate((x0, WALL + 0.3, zt))
        gx = ex - 1.5 if side > 0 else ex - 0.3
        post = post.cut(cq.Workplane("XY").box(1.5 + 0.3, GPS_T + 0.4, Z_COVER - zt + 1, centered=False).translate((gx, gpy - 0.2, zt - 1)))
        # component (patch) face: the lip only covers 0.5 mm of the edge - the board has resistors/LED right up to it
        rx = ex + 0.5 if side < 0 else ex - 1.5 - 0.1
        post = post.cut(cq.Workplane("XY").box(1.0 + 0.1, gpy - 0.2 - (WALL + 0.3) + 0.1, Z_COVER - zt + 1, centered=False).translate((rx, WALL + 0.3 - 0.1, zt - 1)))
        c = c.union(post)
    # (no vent slots in the cover: it is the stack's mounting plate and the handle; the side walls vent)
    c = on_cover(c)
    # in the box frame, so they stay put whichever way round the cover's insides go: the locating lip, the lanyard loop and LOGO_NAME.
    # Lip: a rim 1.2 mm thick standing 2 mm off the inner face, LIP_CLR inside the shell's walls, broken at the four corner bosses, with a
    # 0.8 mm 45-degree lead-in on its outer edge, so the cover drops into place before the screws go in. Prints standing on the plate.
    zt, zb, k = Z_COVER + 0.1, Z_COVER - LIP_H, 0.8
    a = WALL + 2 * BOSS_R + 0.5                                          # clear of the corner bosses and their webs
    for yw, sg in ((WALL + LIP_CLR, 1), (OUT_H - WALL - LIP_CLR, -1)):   # bottom and top walls: profile in (y, z), run along x
        pr = [(yw, zt), (yw, zb + k), (yw + sg * k, zb), (yw + sg * LIP_T, zb), (yw + sg * LIP_T, zt)]
        c = c.union(cq.Workplane("YZ").workplane(offset=a).polyline(pr).close().extrude(OUT_W - 2 * a))
    for xw, sg in ((WALL + LIP_CLR, 1), (OUT_W - WALL - LIP_CLR, -1)):   # side walls: profile in (x, z), run along y (XZ's normal is -Y)
        pr = [(xw, zt), (xw, zb + k), (xw + sg * k, zb), (xw + sg * LIP_T, zb), (xw + sg * LIP_T, zt)]
        c = c.union(cq.Workplane("XZ").workplane(offset=-(OUT_H - a)).polyline(pr).close().extrude(OUT_H - 2 * a))
    # lanyard: one large loop off the middle of the left side edge, in the cover's plane (prints flat, loaded along its layers), 1.5 mm
    # deeper beside the shell's wall where the grip ribs have ended, with a 16 x 5 slot for a strap or a clip
    ly0 = OUT_H / 2 - LANYARD_W / 2
    tab = cq.Workplane("XY").box(LANYARD_OUT + 1, LANYARD_W, WALL, centered=False).translate((-LANYARD_OUT, ly0, Z_COVER)).edges("|Z").edges("<X").fillet(6)
    tab = tab.union(cq.Workplane("XY").box(LANYARD_OUT - 0.3, LANYARD_W, 1.6, centered=False).translate((-LANYARD_OUT, ly0, Z_COVER - 1.5)).edges("|Z").edges("<X").fillet(6))
    tab = tab.cut(cq.Workplane("XY").workplane(offset=Z_COVER - 3).center(-LANYARD_OUT + 4.0 + LANYARD_SLOT[1] / 2, OUT_H / 2).slot2D(LANYARD_SLOT[0], LANYARD_SLOT[1], 90).extrude(WALL + 6))
    c = c.union(tab)
    if LOGO_NAME: c = c.cut(back_name())
    return c

parts = {"front_shell": front_shell(), "mid_frame": mid_frame(), "rear_cover": rear_cover()}
for name, part in parts.items():
    cq.exporters.export(part, str(OUT / f"{name}.stl"), tolerance=0.05, angularTolerance=0.1)
    bb = part.val().BoundingBox()
    print(f"{name:15s} {bb.xlen:6.1f} x {bb.ylen:6.1f} x {bb.zlen:5.1f} mm")
# previews
asm = parts["front_shell"].union(parts["mid_frame"]).union(parts["rear_cover"])
for name, obj, dirn in (("preview_front", parts["front_shell"], (0.6, -0.8, -0.5)), ("preview_rear", parts["rear_cover"], (0.6, -0.8, 0.5)),
                        ("preview_mid", parts["mid_frame"], (0.3, -0.5, -1)), ("preview_assembly", asm, (0.7, -0.9, -0.6))):
    cq.exporters.export(obj, str(OUT / f"{name}.svg"), opt={"width": 900, "height": 600, "projectionDir": dirn, "showHidden": False, "strokeWidth": 0.6})
print("done ->", OUT)

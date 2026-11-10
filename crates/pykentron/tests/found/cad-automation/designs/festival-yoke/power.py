"""Power pod, its lid, the LED strip and a diffuser cover cut to clear the modules.
The pod sits at the back of the collar to balance the cup and phone at the front; one pouch cell and a charge/boost
board feed the LED strip (two halves from the pod to the front gap) and the fan.
Run: uv run --no-project --with cadquery power.py"""
import math

import cadquery as cq

import build
import layout
import module

POD_ANGLE = layout.angle("pod")
BATT = (60, 40, 10)      # pouch cell length (along the rail), width, thickness; ~4000 mAh 3.7 V, protected. Measure yours.
BOARD = (32, 22, 6)      # USB-C charge + 5 V boost board with protection
WALL = 2.4
POD_X = 14               # tray centre, outward from the rail line (puts the lead slots over the LED channel)

STRIP_WIDTH, STRIP_T = build.strip_width(), 2   # 8 mm where the channel allows, narrower on small sizes
STRIP_WATT_PER_M = 10    # 5 V COB strip at full white; many are 8-12 W/m
FAN_WATT = 1.5
CELL_WH = 3.7 * 4.0      # 14.8 Wh nominal
EFFICIENCY = 0.85        # boost converter and cable losses

MODULE_ANGLES = layout.angles()   # modules on the collar in the active kit, for the diffuser cut-outs

_inner = (BATT[1] + 3, BATT[0] + 6, BATT[2] + 1 + BOARD[2] + 3)         # x (across), y (along the rail), z
_outer = (_inner[0] + 2 * WALL, _inner[1] + 2 * WALL)


def _pod_local(angle: float) -> tuple[cq.Workplane, cq.Workplane]:
    _, _, _, _, z_arm = module.z_levels(angle)
    iw, il, ih = _inner
    ow, ol = _outer
    tray = cq.Workplane("XY").box(ow, ol, ih + WALL, centered=(True, True, False)).translate((POD_X, 0, z_arm))
    body = module.clip_local(angle).union(tray)
    body = body.cut(cq.Workplane("XY").box(iw, il, ih + 1, centered=(True, True, False)).translate((POD_X, 0, z_arm + WALL)))
    # corner screw bosses for the lid (M2 self-tapping)
    bx, by = iw / 2 - 3.2, il / 2 - 3.2
    for sx in (-1, 1):
        for sy in (-1, 1):
            post = cq.Workplane("XY").workplane(offset=z_arm + WALL).center(POD_X + sx * bx, sy * by).circle(3.2).extrude(ih)
            body = body.union(post).cut(
                cq.Workplane("XY").workplane(offset=z_arm + ih - 6).center(POD_X + sx * bx, sy * by).circle(0.85).extrude(9))
    # ports in the +y end wall at board height: USB-C, power switch, charge-LED window
    zp = z_arm + WALL + BATT[2] + 1 + BOARD[2] / 2
    for dx, w, h in ((0, 9.6, 3.6), (14, 7.2, 3.2)):
        body = body.cut(cq.Workplane("XY").box(w, WALL + 2, h, centered=(True, False, True)).translate((POD_X + dx, il / 2 - 0.5, zp)))
    body = body.cut(cq.Workplane("XZ").workplane(offset=-(il / 2 - 0.5)).center(POD_X - 14, zp).circle(1.6).extrude(-(WALL + 2)))
    # lead slots through the floor, over the LED channel, one per strip half
    for sy in (-1, 1):
        body = body.cut(cq.Workplane("XY").box(3.4, 8, WALL + 2, centered=(True, True, False)).translate((POD_X + 8, sy * 24, z_arm - 1)))
    # lid
    top = z_arm + WALL + ih
    lid = cq.Workplane("XY").box(ow, ol, WALL, centered=(True, True, False)).translate((POD_X, 0, top))
    rim = cq.Workplane("XY").box(iw - 0.6, il - 0.6, 1.5, centered=(True, True, False)).translate((POD_X, 0, top - 1.5))
    rim = rim.cut(cq.Workplane("XY").box(iw - 3.6, il - 3.6, 1.6, centered=(True, True, False)).translate((POD_X, 0, top - 1.5)))
    lid = lid.union(rim)
    for sx in (-1, 1):
        for sy in (-1, 1):
            lid = lid.cut(cq.Workplane("XY").workplane(offset=top - 2).center(POD_X + sx * bx, sy * by).circle(1.3).extrude(WALL + 4))
    for dy in (-14, 0, 14):  # vents
        lid = lid.cut(cq.Workplane("XY").box(2.4, 12, WALL + 4, centered=(True, True, False)).translate((POD_X, dy, top - 2)))
    return body, lid


def pod_module(angle: float = POD_ANGLE) -> cq.Workplane:
    return module.mount(_pod_local(angle)[0], angle)


def pod_lid(angle: float = POD_ANGLE) -> cq.Workplane:
    return module.mount(_pod_local(angle)[1], angle)


def keepout(angle: float) -> cq.Workplane:
    """Footprint of a module's bar and clip, used to open a hole in the diffuser cover where a module sits."""
    import dock
    f = dock.frame(angle)
    box = cq.Workplane("XY").box(30, dock.BAR_LENGTH + 8, 120, centered=(False, True, True)).translate((-14, 0, -30))
    return dock._place(box, f)


def fitted_diffuser() -> cq.Workplane:
    d = build.diffuser()
    for angle in MODULE_ANGLES.values():
        d = d.cut(keepout(angle))
    return d


def strip() -> cq.Workplane:
    """The LED strip: a ribbon lying on the skin in the channel (8 mm wide in plan), cut at the front gap."""
    u = (build.R["position"] + build.LIP_U) / 2
    plan = build.ring(u, STRIP_WIDTH / 2 - 0.6, STRIP_WIDTH / 2 + 0.6, -80, 160)   # nudged 0.6 mm outward: the cover's inner wall is not a uniform offset round an eccentric ellipse
    above_skin = build.cone(0.4 + STRIP_T).cut(build.cone(0.4))
    ribbon = build.safe_cut(plan.intersect(above_skin), build.gap())
    bars = build.rail()
    return ribbon.cut(bars) if bars is not None else ribbon   # the strip rides over the bar flanges


def strip_length() -> float:
    """Length of one half: from the pod at the back round to the front gap, along the channel centre line."""
    a, b, _ = build.surface((build.R["position"] + build.LIP_U) / 2)
    end = math.acos(min(1.0, (build.Y["front_gap"] / 2 + 2) / a))
    t0, t1, n = -math.pi / 2, end, 2000
    return sum(math.hypot(a * math.sin(t0 + (t1 - t0) * (i + .5) / n), b * math.cos(t0 + (t1 - t0) * (i + .5) / n)) for i in range(n)) * (t1 - t0) / n


def power_notes() -> list[str]:
    half = strip_length()
    total_m = 2 * half / 1000
    full = total_m * STRIP_WATT_PER_M + FAN_WATT
    dim = total_m * STRIP_WATT_PER_M * 0.3 + FAN_WATT * 0.6
    usable = CELL_WH * EFFICIENCY
    return [
        f"LED strip: two halves of {half:.0f} mm from the pod to the front gap, {2 * half:.0f} mm in all, {STRIP_WIDTH} mm wide",
        f"estimated draw: {full:.1f} W strip at full white plus fan, {dim:.1f} W at 30% strip and fan low (assumes {STRIP_WATT_PER_M} W/m strip, {FAN_WATT} W fan)",
        f"estimated runtime on {CELL_WH:.1f} Wh at {EFFICIENCY:.0%} efficiency: about {usable / full:.1f} h at full, {usable / dim:.1f} h dimmed",
    ]


def overlap(a: cq.Workplane, b: cq.Workplane) -> float:
    return build.overlap_volume(a, b)


def main():
    import kit
    if not layout.has_leds():
        print(f"kit {layout.KIT} has no power pod, so no strip, cover or pod to build")
        return
    build.OUT.mkdir(parents=True, exist_ok=True)
    pod, lid, dif, led = pod_module(), pod_lid(), fitted_diffuser(), strip()
    cq.exporters.export(pod, str(build.OUT / "module_pod.stl"))
    cq.exporters.export(lid, str(build.OUT / "module_pod_lid.stl"))
    module.report("module_pod", pod)
    print(f"  pod {_outer[0]:.0f} x {_outer[1]:.0f} mm footprint, holds a {BATT[0]}x{BATT[1]}x{BATT[2]} mm pouch cell and a {BOARD[0]}x{BOARD[1]}x{BOARD[2]} mm board; lid on 4 M2 self-tapping screws, vented")
    for q, (sx, sy) in build.QUARTERS.items():
        cq.exporters.export(build.quarter(dif, sx, sy), str(build.OUT / f"diffuser_fitted_{q}.stl"))
    # everything must clear everything
    mods = kit.modules()
    yoke = build.yoke()
    worst = {}
    for name, m in mods.items():
        worst[f"{name} vs fitted diffuser"] = overlap(m, dif)
        worst[f"{name} vs strip"] = overlap(m, led)
        worst[f"{name} vs yoke"] = overlap(m, yoke)
    worst["lid vs fitted diffuser"] = overlap(lid, dif)
    worst["strip vs yoke"] = overlap(led, yoke)
    worst["strip vs fitted diffuser"] = overlap(led, dif)
    for k, v in worst.items():
        if v > 0.5:
            print(f"WARN {k}: {v:.0f} mm3")
    others = list(mods)
    for i, a in enumerate(others):
        for b in others[i + 1:]:
            v = overlap(mods[a], mods[b])
            if v > 0.5:
                print(f"WARN {a} and {b} overlap {v:.0f} mm3")
    print(f"kit {layout.KIT} clash checks: {sum(1 for v in worst.values() if v <= 0.5)}/{len(worst)} clean")
    for line in power_notes():
        print(line)


if __name__ == "__main__":
    main()

import os

import cadquery as cq
from cadquery import Location as Loc, Vector as Vec

from cq_warehouse.thread import IsoThread

from util import get_bottom, get_center, knurl


class Capsule:
    innerDiameter: float
    outerDiameter: float
    wallThickness: float
    height: float
    capHeight: float
    dent: float
    turns: int
    threadPitch: float
    threadHeight: float
    threadClearance: float
    knurled: bool

    def __init__(
        self,
        innerDiameter=15.0,
        wallThickness=1.5,
        height=60.0,
        dent=4.0,
        turns=6,
        knurled=True,
    ) -> None:
        self.innerDiameter = innerDiameter
        self.wallThickness = wallThickness
        self.outerDiameter = self.innerDiameter + 2 * self.wallThickness
        self.height = height
        self.capHeight = height / 6
        self.dent = dent
        self.turns = turns
        self.threadHeight = 5.0
        self.threadClearance = 0.6
        self.threadPitch = self.threadHeight / self.turns
        self.knurled = knurled

        assert self.height > self.threadHeight
        assert self.wallThickness >= 1.5

    def make(self):
        body = (
            cq.Workplane("XY")
            .circle(self.outerDiameter / 2)
            .circle(self.innerDiameter / 2)
            .extrude(self.height)
        )

        cover = get_bottom(
            self.dent + self.wallThickness, self.outerDiameter, self.wallThickness
        )
        cover = cover.cut(get_bottom(self.dent, self.innerDiameter))

        handle = self.make_handle()
        handle = handle.translate((0, 0, self.height + self.wallThickness))

        body = (
            body.union(cover.mirror("XY"))
            .translate((0, 0, self.wallThickness))
            .union(handle)
        )
        body = body.union(cover.translate((0, 0, self.height)))

        if self.knurled:
            body = body.translate((0, 0, self.dent))
            body = knurl(
                body,
                self.height + 4 * self.wallThickness,
                self.outerDiameter / 2,
                90,
                self.wallThickness / 8,
                180,
                20,
            )
            body = body.translate((0, 0, -self.dent))

        offset = -(self.dent / 2 + self.height - self.capHeight - self.threadHeight)
        cap = body.faces("<Z").workplane(offset).split(keepBottom=True)
        body = body.faces("<Z").workplane(offset).split(keepTop=True)

        md = self.outerDiameter - 2 * (self.wallThickness / 3)

        innerThread = IsoThread(
            major_diameter=md - (self.threadClearance / 2),
            pitch=self.threadPitch,
            length=self.threadHeight,
            external=True,
            end_finishes=("fade", "fade"),
            hand="right",
        )

        innerThread = (
            cq.Workplane("XY", origin=(0, 0, -1))
            .circle(innerThread.min_radius)
            .circle(self.innerDiameter / 2)
            .extrude(self.threadHeight + 1)
            .union(innerThread)
        )

        outerThread = IsoThread(
            major_diameter=md + (self.threadClearance / 2),
            pitch=self.threadPitch,
            length=self.threadHeight - self.threadPitch,
            external=False,
            end_finishes=("fade", "fade"),
            hand="right",
        )
        outerThread = outerThread.cq_object.translate((0, 0, self.threadPitch / 2))

        cap = cap.cut(
            cq.Workplane("XY")
            .cylinder(
                2 * self.threadHeight, (self.outerDiameter / 2) - self.wallThickness / 3
            )
            .translate((0, 0, self.height - self.capHeight - self.threadHeight))
        )

        cap = cap.union(
            outerThread.translate(
                (0, 0, self.height - self.capHeight - self.threadHeight)
            )
        )

        body = body.union(
            innerThread.translate(
                (0, 0, self.height - self.capHeight - self.threadHeight)
            )
        )

        return body, cap

    def make_handle(self):
        base = self.outerDiameter if self.outerDiameter < 15.0 else 15.0
        th = base / 10
        d = (base / 2) - th
        dw = 0.5
        r = 2
        path = (
            cq.Workplane("XZ")
            .moveTo(-d + dw, 0)
            .lineTo(-d + dw, d - r + 2)
            .radiusArc((-(d - r) + dw, d + 2), r)
            .lineTo(d - r - dw, d + 2)
            .radiusArc((d - dw, d - r + 2), r)
            .lineTo(d - dw, 0)
            .close()
        )

        handle = cq.Workplane("YZ").ellipse(base / 8, th).sweep(path, isFrenet=True)

        r = (
            -get_center(
                (0, self.dent + self.wallThickness),
                (self.outerDiameter / 2, self.wallThickness),
            )
            + self.dent
            + self.wallThickness
        )

        handle = handle.cut(
            cq.Workplane("XY").sphere(r).translate((0, 0, -r + self.dent)), tol=1e-2
        )

        return handle


cfg = Capsule(knurled=False, height=68, innerDiameter=15)
objs = cfg.make()

parts = []
z_offsets = [0, cfg.height - cfg.capHeight - cfg.threadHeight]

for i, (obj, c) in enumerate(zip(objs, ["black", "blue"])):
    parts.append(obj.translate(((i + 1) * (cfg.outerDiameter + 20), 0, -z_offsets[i])))
    show_object(obj, options={"alpha": 0.5, "color": c})

show_object(parts)

cq.exporters.export(
    objs[0].union(objs[1].translate((40, 0, -z_offsets[1]))),
    os.path.join("output", "capsule.stl"),
)

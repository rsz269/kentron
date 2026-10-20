import os
from typing import Tuple

import math
import cadquery as cq
from sympy import Point, Line, Polygon, Segment


def cut_connector(sides, side, shape, x_offset, y_offset, depth):
    return (
        sides.faces(side)
        .workplane(centerOption="CenterOfBoundBox")
        .move(x_offset, y_offset)
        .placeSketch(shape)
        .cutBlind(-depth)
    )


def latch(t, h, l, w, a=0.5):
    s = (
        cq.Sketch()
        .polygon([(-t / 2, h / 2), (t / 2, h / 2), (0, -h / 2), (-t / 2, -h / 2)])
        .polygon(
            [(-t / 2, -h / 2), ((-t / 2) - a, (-h / 2) + l), (-t / 2, (-h / 2) + l)]
        )
        .clean()
    )
    return cq.Workplane("YZ").workplane(offset=-w / 2).placeSketch(s).extrude(w)


def slanted_slots(w, h, n, angle, ratio=1.0):
    rect = Polygon(
        *map(
            Point, [(-w / 2, -h / 2), (w / 2, -h / 2), (w / 2, h / 2), (-w / 2, h / 2)]
        )
    )

    if not angle:
        angle = math.degrees(math.atan(w / h))

    dy = h / (n + 1)
    dx = w / (n + 1)

    pts = [(-(w / 2) + dx * (i + 1), (h / 2) - dy * (i + 1)) for i in range(n)]
    slot_height = math.hypot(dx, dy) / 2

    m = cq.Sketch()
    for p in pts:
        line = Line(Point(p), slope=math.tan(math.radians(angle)))
        ints = rect.intersection(line)
        m.push([tuple(Segment(*ints).midpoint)]).slot(
            ints[0].distance(ints[1]), slot_height * ratio, angle
        )

    return m


def hex_mesh(n, m, r=2, mesh_th=2):

    m = 2 * m - 1
    x = r * math.cos(math.radians(30))
    xb = x + mesh_th / 2
    rb = xb * math.cos(math.radians(30))
    dx = rb * 1.5
    dy = rb * math.cos(math.radians(30))

    w = rb * (3 * n - 1)
    h = m * dy

    mesh_points = []

    for y in range(m):
        for x in range(n - 1 if (y % 2 == 1) else n):
            x_s = 0 if (y % 2) == 0 else dx
            mesh_points.append(
                (
                    (x * 2 * dx) + x_s - w / 2 + rb,
                    (y * dy) - h / 2 + dy / 2,
                )
            )

    return cq.Sketch().push(mesh_points).regularPolygon(r, 6, angle=90)


def ear(a, b, diam=3, r=5):
    r1 = r + diam / 2
    s = (
        cq.Sketch()
        .arc((a, -r1), (a + r1, 0), (a, r1))
        .segment((0.0, -b / 2), (0, b / 2))
        .hull()
        .reset()
        .push([(a, 0)])
        .circle(diam / 2, mode="s")
    )
    return s


class Enclosure:
    innerWidth: float
    innerDepth: float
    innerHeight: float
    wallThickness: float
    pcbHolesWidth: float
    pcbHolesDepth: float
    pcbThickness: float
    pcbHolesDiameter: float
    pcbOffset: Tuple[float, float]
    pcbStandoffsHeight: float
    rubberFeetDiam: float
    latchWidth: float
    numLatches: float

    def __init__(
        self,
        innerWidth=45,
        innerDepth=39.99,
        innerHeight=18,
        wallThickness=1.5,
        pcbHolesWidth=40,
        pcbHolesDepth=34.99,
        pcbThickness=1.5,
        pcbHolesDiameter=2.5,
        pcbOffset=(0, 0),
        pcbStandoffsHeight=3.0,
        rubberFeetDiam=10,
        latchWidth=6,
        numLatches=2,
    ) -> None:
        assert innerHeight > (pcbStandoffsHeight + pcbThickness)
        assert wallThickness >= 1.0

        self.innerWidth = innerWidth
        self.innerDepth = innerDepth
        self.innerHeight = innerHeight
        self.wallThickness = wallThickness
        self.pcbHolesWidth = pcbHolesWidth
        self.pcbHolesDepth = pcbHolesDepth
        self.pcbThickness = pcbThickness
        self.pcbHolesDiameter = pcbHolesDiameter
        self.pcbOffset = pcbOffset
        self.pcbStandoffsHeight = pcbStandoffsHeight
        self.rubberFeetDiam = rubberFeetDiam
        self.latchWidth = latchWidth
        self.numLatches = numLatches
        self.connectors = []

    def make(self):
        base_outline = (
            cq.Sketch()
            .rect(self.innerDepth, self.innerWidth)
            .vertices()
            .fillet(2)
            .reset()
        )
        inside_outline = base_outline.copy().wires().offset((0.5)).clean()
        outside_outline = (
            inside_outline.copy().wires().offset((self.wallThickness)).clean()
        )

        pcb = (
            cq.Workplane("XY")
            .workplane(offset=-self.pcbThickness / 2, centerOption="CenterOfBoundBox")
            .move(self.pcbOffset[1], self.pcbOffset[0])
            .sketch()
            .rect(
                self.pcbHolesDepth + self.pcbHolesDiameter + 3,
                self.pcbHolesWidth + self.pcbHolesDiameter + 3,
            )
            .vertices()
            .fillet(2.5)
            .finalize()
            .extrude(self.pcbThickness)
            .faces(">Z")
            .rect(self.pcbHolesDepth, self.pcbHolesWidth, forConstruction=True)
            .vertices()
            .hole(self.pcbHolesDiameter)
        )

        sides = (
            cq.Workplane("XY")
            .placeSketch(outside_outline)
            .extrude(self.innerHeight)
            .faces(">Z")
            .placeSketch(inside_outline)
            .cutThruAll()
        )

        bbox = sides.val().BoundingBox()
        BOX_H = bbox.xmax - bbox.xmin
        BOX_W = bbox.ymax - bbox.ymin

        topbottom = (
            cq.Workplane("XY").placeSketch(outside_outline).extrude(self.wallThickness)
        )

        top = (
            topbottom.faces(">Z")
            .fillet(self.wallThickness * 0.75)
            .faces("<Z")
            .placeSketch(
                inside_outline.copy()
                .rect(self.innerDepth - 2, self.innerWidth - 2, mode="s")
                .vertices()
                .fillet(2)
            )
            .extrude(-2)
            .faces(">X[2]")
            .edges(">Z")
            .fillet(0.5)
        )

        pcb_standoffs = (
            cq.Workplane("XY")
            .workplane(centerOption="CenterOfBoundBox")
            .move(self.pcbOffset[1], self.pcbOffset[0])
            .rect(self.pcbHolesDepth, self.pcbHolesWidth, forConstruction=True)
            .vertices()
            .cylinder(self.pcbStandoffsHeight, (self.pcbHolesDiameter + 3) / 2)
            .faces(">Z")
            .workplane()
            .rect(self.pcbHolesDepth, self.pcbHolesWidth, forConstruction=True)
            .vertices()
            .cboreHole(
                self.pcbHolesDiameter + 0.7, self.pcbHolesDiameter + 1.0, 1, 3
            )  # M2.5 threaded inserts
        )

        bottom = (
            topbottom.faces("<Z")
            .fillet(self.wallThickness * 0.75)
            .faces(">Z")
            .union(
                pcb_standoffs.translate(
                    (0, 0, self.wallThickness + self.pcbStandoffsHeight / 2)
                )
            )
        )

        if self.rubberFeetDiam > 0.0:
            bottom = (
                bottom.faces("<Z")
                .workplane()
                .rect(
                    BOX_H - self.rubberFeetDiam - 2 * self.wallThickness - 2,
                    BOX_W - self.rubberFeetDiam - 2 * self.wallThickness - 2,
                    forConstruction=True,
                )
                .vertices()
                .hole(self.rubberFeetDiam, 0.5)
            )

        a = BOX_H / (self.numLatches)
        for i in range(self.numLatches):
            # cut latch holes
            sides = cut_connector(
                sides,
                ">Y",
                cq.Sketch().rect(self.latchWidth, 2),
                ((BOX_H - a) / 2) - (a * i),
                -self.innerHeight / 2 + self.innerHeight - 4 + 1 - 2,
                self.wallThickness,
            )
            sides = cut_connector(
                sides,
                "<Y",
                cq.Sketch().rect(self.latchWidth, 2),
                ((BOX_H - a) / 2) - (a * i),
                -self.innerHeight / 2 + self.innerHeight - 4 + 1 - 2,
                self.wallThickness,
            )

        for c in self.connectors:
            sides = cut_connector(
                sides,
                c[0],
                c[1],
                c[2],
                c[3] - self.innerHeight / 2,
                c[4] if c[4] else self.wallThickness,
            )

        top = top.translate((0, 0, self.innerHeight + self.wallThickness))

        for i in range(self.numLatches):
            latch_1 = latch(1, 4, 2, self.latchWidth, self.wallThickness / 2).translate(
                (
                    ((BOX_H - a) / 2) - (a * i),
                    -self.innerWidth / 2,
                    self.wallThickness + self.innerHeight - 2 - 2,
                )
            )
            latch_2 = latch_1.rotate((0, 0, 0), (0, 0, 1), 180)
            top = top.union(latch_1).union(latch_2)

        sides = sides.translate((0, 0, self.wallThickness))
        bottom = sides.union(bottom)

        return top, bottom, pcb

    def add_connectors(self, connectors):
        self.connectors = connectors

    def flip_top(self, top):
        d1, d2 = (
            (0, self.innerWidth)
            if self.innerWidth < self.innerDepth
            else (self.innerDepth, 0)
        )
        return top.rotateAboutCenter((1, 0, 0), 180).translate(
            (
                1.2 * d1,
                1.2 * d2,
                -self.innerHeight - self.wallThickness + self.wallThickness / 2,
            )
        )


enclosure = Enclosure()

connectors = [
    (  # USB cable widening
        "<X",
        cq.Sketch().slot(10.5 - 5, 5),
        -((enclosure.innerWidth / 2) - 10),
        enclosure.pcbStandoffsHeight + enclosure.pcbThickness + (3.5 / 2),
        1,
    ),
    (  # USB
        "<X",
        cq.Sketch().slot(9.0 - 3.5, 3.5),
        -((enclosure.innerWidth / 2) - 10),
        enclosure.pcbStandoffsHeight + enclosure.pcbThickness + (3.5 / 2),
        None,
    ),
    (  # solar connector
        "<X",
        cq.Sketch().rect(7, 4),
        (enclosure.innerWidth / 2) - 19 + 3.5,
        enclosure.pcbStandoffsHeight + enclosure.pcbThickness + (4 / 2),
        None,
    ),
    (  #  LEDs hole
        "<X",
        cq.Sketch().slot(8.5, 1),
        -2,
        enclosure.pcbStandoffsHeight + enclosure.pcbThickness + (1.0 / 2),
        None,
    ),
    (  # battery connector
        ">X",
        cq.Sketch().rect(13, 7 - enclosure.pcbThickness),
        (enclosure.innerWidth / 2) - (19 - 3.5),
        enclosure.pcbStandoffsHeight
        + enclosure.pcbThickness
        + ((7 - enclosure.pcbThickness) / 2),
        None,
    ),
    (  # 5V connector
        ">Y",
        cq.Sketch().rect(7, 4),
        (enclosure.innerDepth / 2) - 16 + 3.5,
        enclosure.pcbStandoffsHeight + enclosure.pcbThickness + (4 / 2),
        None,
    ),
    (  # battery screw connector
        ">Y",
        cq.Sketch().rect(7, 4),
        -(enclosure.innerDepth / 2) + 13 - 3.5,
        enclosure.pcbStandoffsHeight + enclosure.pcbThickness + (4 / 2),
        None,
    ),
]

enclosure.add_connectors(connectors)
top, bottom, pcb = enclosure.make()

ear_outline = ear(6, 20, 3, 3)
mount = (
    cq.Workplane("XY")
    .workplane(centerOption="CenterOfBoundBox")
    .placeSketch(ear_outline)
    .extrude(3)
    .rotate((0, 0, 0), (0, 0, 1), 90)
    .translate((0, enclosure.innerWidth / 2 + 0.5, 0))
)
mount = mount.union(mount.rotate((0, 0, 0), (0, 0, 1), 180))
bottom = bottom.union(mount)

# mesh = slanted_slots(
#     enclosure.innerDepth - 8, enclosure.innerWidth / 2 - 8, 16, None, 0.5
# )
mesh = hex_mesh(6, 5, 1.5, 2)

mesh_cut = (
    cq.Workplane("XY")
    .workplane(offset=enclosure.wallThickness + enclosure.innerHeight)
    .placeSketch(mesh)
    .extrude(enclosure.wallThickness)
)

top = top.cut(mesh_cut.translate((0, (enclosure.innerWidth - 6) / 4, 0))).cut(
    mesh_cut.translate((0, (-enclosure.innerWidth + 6) / 4, 0))
)

top = enclosure.flip_top(top)

show_object(bottom, options={"alpha": 0.5, "color": (0, 0, 1)})
show_object(
    pcb.translate(
        (
            0,
            0,
            enclosure.wallThickness
            + enclosure.pcbStandoffsHeight
            + (enclosure.pcbThickness / 2),
        )
    ),
    options={"alpha": 0.7},
)
show_object(top, options={"alpha": 0.7, "color": (0, 0, 255)})

cq.exporters.export(
    bottom.union(top), os.path.join("output", "{}.stl".format("enclosure"))
)

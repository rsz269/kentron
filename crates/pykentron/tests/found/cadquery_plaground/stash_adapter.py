import math

from matplotlib import colors

import cadquery as cq
from cadquery import Location as Loc, Vector as Vec

from cq_warehouse.thread import IsoThread

from util import get_bottom, knurl


def sharp_ring(diameter, width, angle):
    radius = diameter / 2
    x = math.tan(math.radians(angle / 2)) * width
    return (
        cq.Workplane("XZ", origin=(radius, 0, 0))
        .polyline([(0, x), (-width, 0), (0, -x)])
        .close()
        .revolve(360, (-radius, 0, 0), (-radius, -1, 0))
    )


class Adapter:
    outerDiameter: float
    wallThickness: float
    height: float
    dent: float
    containerWallThickness: float
    pitch: float
    threadClearance: float
    knurled: bool

    def __init__(
        self,
        outerDiameter=62.0,
        wallThickness=2.5,
        height=25.0,
        dent=11,
        containerWallThickness=0.2,
        threadClearance=0.8,
        knurled=True,
    ):
        self.outerDiameter = outerDiameter
        self.wallThickness = wallThickness
        self.height = height
        self.dent = dent
        self.threadDiameter = self.outerDiameter - 1.5 * self.wallThickness
        self.dent_radius = ((self.outerDiameter**2) + (4 * self.dent**2)) / (
            8 * self.dent
        )
        self.containerWallThickness = containerWallThickness
        self.pitch = 6.0
        self.threadClearance = threadClearance
        self.knurled = knurled

    def make(self):
        container_side = (
            cq.Workplane("XY")
            .circle(self.outerDiameter / 2)
            .circle((self.outerDiameter / 2) + self.containerWallThickness)
            .extrude(2.0 * self.height)
        )

        bottom = get_bottom(
            self.dent + self.containerWallThickness,
            self.outerDiameter,
            self.containerWallThickness,
        )

        container_bottom = bottom.cut(
            get_bottom(
                self.dent,
                self.outerDiameter - self.containerWallThickness,
            )
        )

        container = (
            container_side.union(container_bottom)
            .faces("<Z")
            .edges()
            .fillet(self.containerWallThickness / 2 - 0.05)
        )

        thread = IsoThread(
            major_diameter=self.threadDiameter - (self.threadClearance / 2),
            pitch=self.pitch,
            length=self.height,
            external=True,
            end_finishes=("fade", "fade"),
            hand="right",
        )

        self.internalDent = self.dent_radius - (
            0.5 * math.sqrt((4 * self.dent_radius**2) - (4 * thread.min_radius**2))
        )

        core = (
            cq.Workplane("XY")
            .circle(thread.min_radius)
            .extrude(self.height + 3)
            .faces(">Z")
            .hole(2 * thread.min_radius - 4, self.height + 3 - self.internalDent - 2)
        )

        cap_internal = core.union(thread)
        cap_internal = cap_internal.cut(
            get_bottom(self.internalDent, 2 * thread.min_radius)
        )

        thread = IsoThread(
            major_diameter=self.threadDiameter + (self.threadClearance / 2),
            pitch=self.pitch,
            length=self.height,
            external=False,
            end_finishes=("fade", "fade"),
            hand="right",
        )

        self.assembly_offsets = [
            0,
            adapter.dent
            - adapter.internalDent
            + adapter.pitch / 2
            + adapter.containerWallThickness,
            adapter.dent + adapter.containerWallThickness - adapter.internalDent,
        ]

        cap_external = (
            cq.Workplane("XY")
            .workplane(offset=-self.assembly_offsets[1])
            .circle(self.outerDiameter / 2)
            .circle(self.threadDiameter / 2)
            .extrude(self.height + self.assembly_offsets[1] + 10)
            .faces(">Z")
            .wires(cq.selectors.RadiusNthSelector(0))
            .chamfer(8, self.wallThickness - 1)
            .cut(bottom.translate((0, 0, -self.assembly_offsets[1])))
        )

        cap_external = cap_external.union(thread)

        if self.knurled:
            cap_external = knurl(
                cap_external,
                self.height,
                self.outerDiameter / 2,
                90,
                self.wallThickness / 4,
                40,
                36,
            )

        return [container, cap_external, cap_internal]


adapter = Adapter()
elements = adapter.make()

split = []
parts = []

for i, obj in enumerate(elements):
    split.append(
        obj.faces(">Z")
        .workplane()
        .transformed(rotate=cq.Vector(90, 0, 0))
        .split(keepBottom=True)
        .translate((0, 0, adapter.assembly_offsets[i]))
    )

cap_external = elements[1].translate((0, 0, adapter.assembly_offsets[1]))
cap_internal = elements[2]

for i, obj in enumerate([cap_external, cap_internal]):
    parts.append(obj.translate(((i + 1) * (adapter.outerDiameter + 20), 0, 0)))

for s, c in zip(split, map(lambda cs: colors.to_rgba(cs), ["blue", "green", "yellow"])):
    show_object(s, options={"alpha": 0.5, "color": c})

show_object(parts)

cq.exporters.export(
    cap_external.union(cap_internal.translate((100, 0, 0))), "output/stash.stl"
)

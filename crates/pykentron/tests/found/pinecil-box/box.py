import cadquery as cq
from types import SimpleNamespace as d
import math
import itertools

dims = d()
dims.tip = d()
dims.tip.total_length = 105
dims.tip.welding_tip = d()
dims.tip.welding_tip.length = 13.5  # for the skinny one, rest are 11
dims.tip.welding_tip.diam = 4.6
dims.tip.cone_base_to_welding_tip = 55
dims.tip.heater_diam = 6.5
dims.tip.cone = d()
dims.tip.cone.od = 11
dims.tip.cone.id = 5
dims.tip.cone.sloped_length = 3
dims.tip.cone.length_from_cone_to_body = 3.5
dims.tip.contacts = d()
dims.tip.contacts.length = 34
dims.body = d()
dims.body.width = 14.5
dims.body.height = 17.4
dims.body.length = 103.8
dims.body.neg_height = dims.tip.cone.od / 2 + 2.2  # distance from the centre axis of the tip to the bottom of the body

dims.clear_radius = 0.5

dims.slot_radius = 17

tip_cutout = (
    cq.Workplane()
    .workplane(offset=-5)
    .circle(dims.tip.welding_tip.diam / 2)
    .workplane(offset=5)
    .circle(dims.tip.welding_tip.diam / 2)
    .workplane(offset=(dims.tip.welding_tip.length))
    .circle(dims.tip.heater_diam / 2 + dims.clear_radius)
    .workplane(offset=(dims.tip.cone_base_to_welding_tip))
    .circle(dims.tip.heater_diam / 2 + dims.clear_radius)
    .workplane(offset=(dims.tip.cone.sloped_length))
    .circle(dims.tip.cone.od / 2 + dims.clear_radius)
    .loft(ruled=True)
    .tag("just_cone")
    .faces(">Z")
    .workplane()
    .circle(dims.tip.cone.od / 2 + dims.clear_radius)
    .extrude(dims.tip.contacts.length)
)

iron_cutout = (
    tip_cutout
    ._getTagged("just_cone")
    .faces(">Z")
    .workplane()
    .center(0, -dims.body.neg_height - dims.clear_radius)
    .rect(
        dims.body.width + dims.clear_radius * 2,
        dims.body.height + dims.clear_radius * 2,
        centered=(True, False)
    )
    .extrude(dims.body.length)
)

all_cutouts = (
    cq.Assembly()
    .add(iron_cutout, name="iron")
)
for idx in range(4):
    name = f"tip{idx}"
    all_cutouts.add(tip_cutout, name=name)
    all_cutouts.constrain(
        "iron",
        cq.Face.makePlane(1, 1, dir=(0, 0, 1)),
        name,
        cq.Face.makePlane(1, 1, dir=(0, 0, 1)),
        "Axis",
        # param=math.pi / 2
    )
    all_cutouts.constrain(
        "iron",
        cq.Vertex.makeVertex(
            math.sin(math.pi * idx / 2) * dims.slot_radius,
            math.cos(math.pi * idx / 2) * dims.slot_radius,
            0
        ),
        name,
        tip_cutout.faces(">Z").val(),
        "Point"
    )

all_cutouts.solve()

acbb = all_cutouts.toCompound().BoundingBox()
main_body_rad = (max(
    acbb.xmax,
    acbb.ymax,
    -acbb.xmin,
    -acbb.ymin
) + 3)  # 3mm wall thickness

main_body = (
    cq.Workplane()
    .circle(main_body_rad)
    .extrude(acbb.zmax)
    .cut(all_cutouts.toCompound())
    .faces("<Z")
    .workplane(offset=(-(
        dims.tip.total_length + 4
    )))
    .tag("tip_cutout")
)
for ang in [x * 90 for x in range(4)]:
    main_body = (
        main_body
        .workplaneFromTagged("tip_cutout")
        .transformed(rotate=(0, 0, ang))
        .transformed(offset=(0, dims.slot_radius, 0))
        .rect(dims.tip.heater_diam + 1, main_body_rad, centered=(True, False))
        .cutBlind(dims.tip.welding_tip.length + 4)
    )

dims.top_cap = d()
dims.top_cap.length = 20
dims.top_cap.wallthick = 2
icbb = iron_cutout.val().BoundingBox()
iron_rad = 0
for xval, yval in itertools.product([icbb.xmax, icbb.xmin], [icbb.ymax, icbb.ymin]):
    iron_rad = max(iron_rad, math.sqrt(xval ** 2 + yval ** 2))

dims.top_cap.id = 2 * (iron_rad + 3)
dims.top_cap.od = dims.top_cap.id + dims.top_cap.wallthick * 2
main_body = (main_body.cut(
    cq.Workplane("XZ", origin=main_body.faces(">Z").val().Center().toTuple())
    .move(dims.top_cap.id / 2, 0)
    .vLine(-dims.top_cap.length)
    .hLine(dims.top_cap.wallthick)
    .lineTo(main_body_rad, -60)
    .hLine(10)
    .vLineTo(0)
    .close()
    .revolve(axisStart=(0, 0), axisEnd=(0, 1))
))

top_cap = (
    cq.Workplane()
    .circle(dims.top_cap.od / 2)
    .extrude(-dims.top_cap.wallthick)
    .faces(">Z")
    .workplane()
    .tag("cap_underside")
    .circle(dims.top_cap.od / 2)
    .extrude(dims.top_cap.length)
    .tag("solid_body")
)
top_cap = top_cap.cut(
    top_cap
    .workplaneFromTagged("cap_underside")
    .circle(dims.top_cap.id / 2)
    .faces(">Z", tag="solid_body")
    .workplane()
    .circle(dims.top_cap.id / 2 - 0.1)
    .loft(ruled=True, combine=False)
)
top_cap = (
    top_cap
    .faces("<Z")
    .workplane()
    .center(0, dims.top_cap.id / 2 - 3)
    .tag("cut_here")
    .circle(2)
    .cutThruAll()
    .workplaneFromTagged("cut_here")
    .rect(2, 10, centered=(True, False))
    .cutThruAll()
)

main_body.val().exportStl("body.stl")
top_cap.val().exportStl("top_cap.stl")

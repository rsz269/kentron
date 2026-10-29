import cadquery as cq

body = cq.Workplane("XY").workplane(offset=10).circle(30).circle(15).extrude(30)
body = body.edges().fillet(2)
for i in range(12):
    tooth = (cq.Workplane("XY").box(5.5, 8, 20, centered=(False, True, False))
             .translate((27, 0, 0)).rotate((0, 0, 0), (0, 0, 1), 30 * i))
    body = body.union(tooth)
solid = body

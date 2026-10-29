import cadquery as cq

body = cq.Workplane("XY").workplane(offset=10).circle(30).extrude(30)
body = body.faces(">Z").edges().fillet(2)
for i in range(12):
    tooth = (cq.Workplane("XY").box(6.5, 8, 20, centered=(False, True, False))
             .translate((26, 0, 0)).rotate((0, 0, 0), (0, 0, 1), 30 * i))
    body = body.union(tooth)
body = body.cut(cq.Workplane("XY").circle(15).extrude(40))
body = body.faces(">Z").edges(cq.selectors.RadiusNthSelector(0)).chamfer(2)
solid = body

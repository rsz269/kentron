import cadquery as cq
# 1.50 x 1.125 x 0.1875 plate with "Unity" cut through
plate = cq.Workplane("XY").box(1.5, 1.125, 0.1875)
text = (cq.Workplane("XY").workplane(offset=-0.2)
        .text("Unity", 0.42, 0.6, font="Arial", kind="bold", halign="center", valign="center"))
solid = plate.cut(text)

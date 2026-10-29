import cadquery as cq
L, W, T = 1.5, 1.125, 0.1875
plate = cq.Workplane("XY").box(L, W, T, centered=(True, True, False))
txt = (cq.Workplane("XY").workplane(offset=-0.05)
       .text("Unity", 0.45, T + 0.1, kind="bold", halign="center", valign="center", combine=False))
bb = txt.val().BoundingBox()
txt = txt.translate((-(bb.xmin + bb.xmax) / 2, -(bb.ymin + bb.ymax) / 2 - 0.03, 0))
solid = plate.cut(txt)

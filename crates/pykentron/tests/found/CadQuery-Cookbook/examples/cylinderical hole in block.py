#a flat square plate with a hole in the middle, sides 20mm, thickness 5mm
import cadquery as cq
# <parameters>
thickness = 5
width = 20
# </parameters>

result = (
    cq.Workplane("front")
    .box(width, width, thickness)
    .faces(">Z")
    .hole(thickness)
)
cq.exporters.export(result, 'result.stl')
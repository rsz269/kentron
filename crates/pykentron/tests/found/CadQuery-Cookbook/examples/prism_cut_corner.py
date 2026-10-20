#A rectangular plate with length 3mm, height 2mm and thickness of 0.5mm 
#The face parallel to the XY plane was chosen and the intersection of the X-axis and Y-axis was taken to determine the vertics
#At the center of the mass of the object a circle was extruded and then subtracted from the main object
#The circle has a diameter of 1mm
import cadquery as cq

# <parameters>
length = 3
height = 2
thickness = 0.5
diameter = 1
# </parameters>

# Create a basic prism and cut out a corner
result = (
    cq.Workplane("front")
    .box(length, height, thickness)
    .faces(">Z")
    .vertices("<XY")
    .workplane(centerOption="CenterOfMass")
    .circle(diameter)
    .cutThruAll()
)
cq.exporters.export(result, 'result.stl')
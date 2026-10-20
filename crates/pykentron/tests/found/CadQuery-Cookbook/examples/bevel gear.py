##A bevel Gear with number of teeth 11, cone angle 45 degrees, a face width of 3mm and a bore diameter of 4mm
import cadquery as cq
from cq_gears import BevelGear, BevelGearPair

# <parameters>
module = 1.0
teeth_number = 11
cone_angle = 45.0
face_width = 3.0
bore_d = 4.0
# </parameters>

bevel_gear = BevelGear(module=module, teeth_number=teeth_number,
                       cone_angle=cone_angle, face_width=face_width, bore_d=bore_d)

result = cq.Workplane('XY').gear(bevel_gear)
cq.exporters.export(result, 'result.stl')
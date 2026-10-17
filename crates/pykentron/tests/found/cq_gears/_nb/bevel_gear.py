import cadquery as cq
pass  # from jupyter_cadquery import show, set_defaults

set_defaults(axes=False, axes0=True, grid=(True, False, False),
             default_color='lightgrey', tools=False)

from cq_gears import BevelGear, BevelGearPair

bevel_gear = BevelGear(module=1.0, teeth_number=11,
                       cone_angle=45.0, face_width=3.0, bore_d=4.0)

wp = cq.Workplane('XY').gear(bevel_gear)
# show(wp)

bevel_gear = BevelGear(module=1.0, teeth_number=19, cone_angle=45.0,
                       face_width=7.0, helix_angle=30.0, bore_d=6.0)

wp = cq.Workplane('XY').gear(bevel_gear)
# show(wp)

bevel_gears = BevelGearPair(module=1.0, gear_teeth=40, pinion_teeth=20,
                            face_width=6.0, helix_angle=42.0, axis_angle=90.0,
                            bore_d=22.0, pinion_build_args={'bore_d': 6.0})

wp = cq.Workplane('XY').gear(bevel_gears)
# show(wp)


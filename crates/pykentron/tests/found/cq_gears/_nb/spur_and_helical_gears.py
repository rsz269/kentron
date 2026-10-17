import cadquery as cq
pass  # from jupyter_cadquery import show, set_defaults

set_defaults(axes=False, axes0=True, grid=(True, True, True),
             default_color='lightgrey', tools=False)

from cq_gears import SpurGear, HerringboneGear

spur_gear = SpurGear(module=1.0, teeth_number=19, width=5.0, bore_d=5.0)

wp = cq.Workplane('XY').gear(spur_gear)
# show(wp)

spur_gear = SpurGear(module=1.0, teeth_number=13, width=5.0, bore_d=5.0)

wp = (cq.Workplane('XY')
        # Pushing an array of 4 points with spacing equal to the gear's pitch diameter 
      .rarray(xSpacing=spur_gear.r0 * 2.0,
              ySpacing=1.0, xCount=4, yCount=1, center=False)
        # Create 4 meshing gears
      .gear(spur_gear)
      
      .moveTo(spur_gear.r0 * 2 * 4, 0.0)
        # Create an additional gear with the same profile but with different
        # bore diameter and hub
      .addGear(spur_gear, bore_d=3.0, hub_d=8.0, hub_length=4.0))

# show(wp)

# Create two meshing helical gears with different tooth count and face width

helical_gear1 = SpurGear(module=1.0, teeth_number=29, width=16.0,
                         # the sign of helix_angle determines direction:
                         # positive=clockwise, negative=counterclockwise
                         helix_angle=40.0,
                         bore_d=20.0)

helical_gear2 = SpurGear(module=1.0, teeth_number=17, width=8.0,
                         helix_angle=-40.0,
                         bore_d=10.0)

wp = (cq.Workplane('XY')
      .gear(helical_gear1)
      .moveTo(helical_gear1.r0 + helical_gear2.r0, 0.0)
      .addGear(helical_gear2))

# show(wp)

# Herringbone helical gear variation

hb_gear = HerringboneGear(module=1.0, teeth_number=42, width=10.0,
                          helix_angle=-20.0, bore_d=28.0)

wp = cq.Workplane('XY').gear(hb_gear)
# show(wp)

# Build a gear with some missing teeth
beaten_spur_gear = SpurGear(module=1.0, teeth_number=20, width=5.0,
                            pressure_angle=20.0, bore_d=5.0,
                            missing_teeth=(0, 10))

wp = cq.Workplane('XY').gear(beaten_spur_gear)
# show(wp)

helical_gear = SpurGear(module=1.0, teeth_number=60, width=8.0,
                        pressure_angle=20.0, helix_angle=45.0,
                        bore_d=10.0, hub_d=16.0, hub_length=10.0,
                        recess_d=52.0, recess=3.0, n_spokes=5,
                        spoke_width=6.0, spoke_fillet=4.0,
                        spokes_id=23.0, spokes_od=48.0)

wp = cq.Workplane('XY').gear(helical_gear)
# show(wp)


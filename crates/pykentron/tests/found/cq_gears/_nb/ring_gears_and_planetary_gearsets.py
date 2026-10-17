import cadquery as cq
pass  # from jupyter_cadquery import show, set_defaults

set_defaults(axes=False, axes0=True, grid=(True, True, True),
             default_color='lightgrey', tools=False)

from cq_gears import (RingGear, HerringboneRingGear,
                      PlanetaryGearset, HerringbonePlanetaryGearset)

straight_ring = RingGear(module=1.0, teeth_number=32, width=8.0, rim_width=3.0)

helical_ring = RingGear(module=1.0, teeth_number=32, width=8.0,
                        rim_width=3.0, helix_angle=30.0)

hb_ring = HerringboneRingGear(module=1.0, teeth_number=32, width=8.0,
                              rim_width=3.0, helix_angle=30.0)


wp = (cq.Workplane('XY')
      .gear(straight_ring)
      .moveTo(50.0, 0.0)
      .addGear(helical_ring)
      .moveTo(100.0, 0.0)
      .addGear(hb_ring))

# show(wp)

gearset = PlanetaryGearset(module=1.0,
                           sun_teeth_number=12, planet_teeth_number=18,
                           width=10.0, rim_width=3.0, n_planets=3,
                           bore_d=6.0)

wp = cq.Workplane('XY').gear(gearset)
# show(wp)

gearset = HerringbonePlanetaryGearset(module=1.0,
                                      sun_teeth_number=18,
                                      planet_teeth_number=18,
                                      width=10.0, rim_width=3.0, n_planets=4,
                                      helix_angle=30.0,
                            # Set backlash and clearance for 3d-prinatability
                                      backlash=0.3, clearance=0.2, bore_d=6.0)

wp = cq.Workplane('XY').gear(gearset)
# show(wp)

# Uncomment the following line to export as .stl
# wp.val().exportStl('./hb_planetary_gearset.stl')

# Internal gearset
gearset1 = HerringbonePlanetaryGearset(module=2.0,
                                       sun_teeth_number=13,
                                       planet_teeth_number=7,
                                       width=10.0, rim_width=1.0, n_planets=5,
                                       helix_angle=30.0,
                            # Set backlash and clearance for 3d-prinatability
                                       backlash=0.3, clearance=0.2,
                                       bore_d=4.0, 
                            # build args could be set for sun/planets separately
                                       sun_build_args={'bore_d': 15.0})

# External gearset
gearset2 = HerringbonePlanetaryGearset(module=2.0,
                                       sun_teeth_number=36,
                                       planet_teeth_number=8,
                                       width=10.0, rim_width=4.0, n_planets=11,
                                       helix_angle=30.0,
                                       backlash=0.3, clearance=0.2,
                                       bore_d=6.0,
                                       sun_build_args={'bore_d': 
                                              gearset1.ring.rim_r * 2.0 - 1.0})

wp = cq.Workplane('XY').gear(gearset1)      
wp2 = cq.Workplane('XY').gear(gearset2)

# Combining them together may take some significant amount of time
wp = wp.add(wp2).combine()

# show(wp)

# Uncomment the following line to export as .stl
# wp.val().exportStl('./double_planetary_gearset.stl')


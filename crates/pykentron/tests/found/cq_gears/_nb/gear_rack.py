import cadquery as cq
pass  # from jupyter_cadquery import show, set_defaults

set_defaults(axes=False, axes0=True, grid=(True, True, True),
             default_color='lightgrey', tools=False)

from cq_gears import SpurGear, HerringboneGear, RackGear, HerringboneRackGear
from math import pi

# Straight rack and pinion
rack = RackGear(module=1.0, length=150.0, width=10.0, height=20.0)
pinion = SpurGear(module=1.0, teeth_number=34, width=10.0, backlash=0.1)


# Helical rack and pinion
hc_rack = RackGear(module=1.0, length=150.0, width=10.0, height=20.0,
                   helix_angle=30.0)
hc_pinion = SpurGear(module=1.0, teeth_number=34, width=10.0, backlash=0.1,
                     helix_angle=-30.0)

# Herringbone rack and pinion
hb_rack = HerringboneRackGear(module=1.0, length=150.0, width=10.0, height=20.0,
                              helix_angle=30.0)
hb_pinion = HerringboneGear(module=1.0, teeth_number=34, width=10.0, backlash=0.1,
                            helix_angle=-30.0)


# Build and arrange different racks and pinions in one space
spacing = 70.0

wp = (cq.Workplane('XY')
      .addGear(rack)
      .moveTo(pi * pinion.r0 * 2.0 / pinion.z * 10, pinion.r0)
      .addGear(pinion, bore_d=12.0))

wp = (wp.moveTo(0.0, spacing)
      .addGear(hc_rack)
      .moveTo(pi * hc_pinion.r0 * 2.0 / hc_pinion.z * 23,
              hc_pinion.r0 + spacing)
      .addGear(hc_pinion, bore_d=12.0))

wp = (wp.moveTo(0.0, spacing * 2.0)
      .addGear(hb_rack)
      .moveTo(pi * hb_pinion.r0 * 2.0 / hb_pinion.z * 38,
              hb_pinion.r0 + spacing * 2.0)
      .addGear(hb_pinion, bore_d=12.0))

# show(wp)


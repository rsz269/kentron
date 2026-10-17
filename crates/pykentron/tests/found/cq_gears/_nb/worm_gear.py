import cadquery as cq
pass  # from jupyter_cadquery import show, set_defaults

set_defaults(axes=False, axes0=True, grid=(True, True, True),
             default_color='lightgrey', tools=False)

from cq_gears import Worm

worm = Worm(module=1.0, lead_angle=7.0, n_threads=2, length=60.0)
worm2 = Worm(module=1.0, lead_angle=-7.0, n_threads=2, length=60.0)

wp = (cq.Workplane('XY')
      .gear(worm, bore_d=8.0)
      .moveTo(0.0, 40.0)
      .addGear(worm2, bore_d=8.0))

# show(wp)


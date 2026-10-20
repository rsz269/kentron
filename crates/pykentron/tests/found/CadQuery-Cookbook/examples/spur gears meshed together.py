import cadquery as cq
from cq_gears import (SpurGear, HerringboneGear, RackGear, HerringboneRackGear,
                      PlanetaryGearset, HerringbonePlanetaryGearset,
                      BevelGearPair, Worm, CrossedGearPair, HyperbolicGearPair)


L = lambda *args: cq.Location(cq.Vector(*args))
Lrz = lambda loc, r: cq.Location(cq.Vector(*loc), cq.Vector(0.0, 0.0, 1.0), r)
C = lambda *args: cq.Color(*args)

def animate_pinion_and_gear(pinion, gear):
    ratio = pinion.z / gear.z

    asm = (cq.Assembly(name='gears')
           .add(pinion.build(), name='pinion', color=C('goldenrod'),
                loc=L(0.0, 0.0, 0.0))
           .add(gear.build(), name='gear', color=C('lightsteelblue'),
                loc=L(gear.r0 + pinion.r0, 0.0, 0.0)))

    
    return asm.toCompound()

pinion = SpurGear(module=2.0, teeth_number=17, width=10.0, helix_angle=0.0,
                  bore_d=10.0, chamfer=0.5)
gear = SpurGear(module=2.0, teeth_number=29, width=10.0, helix_angle=0.0,
                bore_d=10.0, chamfer=0.5, n_spokes=4, spokes_id=20.0,
                spokes_od=44.0, spoke_width=6.0, spoke_fillet=4.0)

result = animate_pinion_and_gear(pinion, gear)
cq.exporters.export(result, 'result.stl')
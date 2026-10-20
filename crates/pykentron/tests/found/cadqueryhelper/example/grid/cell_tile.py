import cadquery as cq
from cadqueryhelper.grid import cell_tile

points:list[tuple[float,float]] = [
    (0,0),
    (0,25),
    (25,25),
    (25,0)
]

def plain(length,width,height):
    log('plain')
    return cq.Workplane("XY").box(length,width,height)

ex_cell = cell_tile(
    points=points,
    height=3,
    tile_method = plain
)

#show_object(ex_cell)
cq.exporters.export(ex_cell,'stl/grid_cell_tile.stl')
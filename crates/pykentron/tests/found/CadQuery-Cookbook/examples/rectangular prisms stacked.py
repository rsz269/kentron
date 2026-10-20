
#Create 3D models of boxes with slits cut into them, demonstrating the difference between using `combine=False` and `combine=True` when performing the cut operation.
#
#### Box Construction
#1. **Shape**: A rectangular box.
#2. **Parameters**: Width, depth, height.
#
#### Slit Construction
#1. **Shape**: A rectangular box representing the slit.
#2. **Positioning**: Translate the slit box to be centered on the top face of the main box.
#3. **Parameters**: Width (same as main box), depth, height, combine flag.
#
#### Assembly
#1. Create an assembly containing the original box and two variations of the box with slits.
#2. Position the boxes side by side.
import cadquery as cq

# <parameters>
box_width = 5
box_depth = 10
box_height = 10
slit_depth = 3
slit_height = 7
# </parameters>

def make_box(width, depth, height):
    return cq.Workplane().box(width, depth, height)

def make_slit(box, depth, height, combine=False):
    return (
        box.faces(">Z")
        .box(box_width, depth, height, combine=combine)
        .translate((0, 0, -height / 2))
    )

box = make_box(box_width, box_depth, box_height)

# Create slit using combine=False
slit1 = make_slit(box, slit_depth, slit_height, combine=False)
body1 = box.cut(slit1)

# Create slit using combine=True
slit2 = make_slit(box, slit_depth, slit_height, combine=True)
body2 = box.cut(slit2)

# Create an assembly to showcase the different results
result = (
    cq.Assembly()
    .add(box, loc=cq.Location((0, 0, 0)))
    .add(body1, loc=cq.Location((box_width * 2, 0, 0)))
    .add(slit1, loc=cq.Location((box_width * 2, 0, 0)))
    .add(body2, loc=cq.Location((box_width * 4, 0, 0)))
    .add(slit2, loc=cq.Location((box_width * 4, 0, 0)))
)

filename = "render.stl"
result.save(filename, exportType="STL")
cq.exporters.export(result.toCompound(), 'result.stl')
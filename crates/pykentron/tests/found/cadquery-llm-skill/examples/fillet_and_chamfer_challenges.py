# Example: Shows the common failure modes of both fillets and chamfers

import cadquery as cq

# Parameters
box_width = 30
box_depth = 20
box_height = 20

# Base box that can have fillets or chamfers applied
result = cq.Workplane("XY").box(box_width, box_depth, box_height)

# This fillet works, but if it is increased at all, it will fail.
# The 0.0001 being subtracted is to keep the fillets from touching.
# Once they do, there will be a "BRep_API: command not done" error
# The two selected edges are box_depth apart, so a radius of box_depth / 2
# is where the fillets just touch - any larger and they overlap, causing
# the error.
result = result.faces(">Z").edges("|X").fillet(box_depth / 2.0 - 0.0001)

# The following line throws the error if uncommented
# result = result.faces(">Z").edges("|X").fillet(box_depth / 2.0)

# Resetting the box
result = cq.Workplane("XY").box(box_width, box_depth, box_height)

# A similar problem happens with chamfers
# Again, 0.0001 is subtracted so that the chamfers do not touch.
# If they do touch, there will be a "BRep_API: command not done" error
result = result.faces(">Z").edges("|X").chamfer(box_depth / 2.0 - 0.0001)

# The following lines throws the BRep_API error if uncommented
# result = result.faces(">Z").edges("|X").chamfer(box_depth / 2.0)

# The following call can be used if running in CQ-editor
# show_object(result)

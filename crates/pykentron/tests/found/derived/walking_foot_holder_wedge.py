# Derived script: not one of the hunted corpus's 500 (see derived.tsv and ../README.md).
#
# Derived from: https://github.com/kmeisthax/gridfinity-cadquery at 637a5c4,
#   examples/walking_foot_holder.py, MIT License, Copyright ©2022-2024 kmeisthax and
#   ©2022-2024 Zach Freedman (the license follows below). It imports the project's own
#   gridfinity.py (run it with that project's directory on the path).
#
# Change (2026-10-27): the original's asymmetric chamfer of the tool's front bottom edge,
# .chamfer(H / 4, D / 3), puts its first length on the edge's first face in the kernel's face
# order (OCCT's: the bottom; kentron's: the front), so the original's result depends on face
# order and is left out of the corpus counts. Here the wedge OCCT's chamfer cuts is cut
# explicitly (OCCT's result is unchanged); the rest of the model (fillets of 40 and 1 on the
# tool, its cut from a gridfinity block, the four fillets at the top) stays a comparison.
#
# Copyright ©2022-2024 kmeisthax.
#
# Copyright ©2022-2024 Zach Freedman.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the “Software”), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Holder for the walking foot that came with my Brother CS7000X"""

import cadquery as cq
import gridfinity

walking_foot_block = cq.Workplane("XY")\
    .gridfinity_block(1, 1, 3)\
    .gridfinity_block_stack(1, 1)\
    .gridfinity_block_lip(1, 1)

FUDGE_FACTOR = 0.5

WALKING_FOOT_WIDTH = 19.62 + FUDGE_FACTOR
WALKING_FOOT_HEIGHT = 27.35 + FUDGE_FACTOR
WALKING_FOOT_DEPTH = 42.56 + FUDGE_FACTOR

walking_foot = cq.Workplane("XY")\
    .rect(WALKING_FOOT_WIDTH, WALKING_FOOT_HEIGHT)\
    .extrude(WALKING_FOOT_DEPTH)\
    .edges(">Y")\
    .edges("<Z")\
    .fillet(WALKING_FOOT_HEIGHT / 2.5)

# (Derived: the original chamfers this edge with .chamfer(WALKING_FOOT_HEIGHT / 4,
# WALKING_FOOT_DEPTH / 3), which measures the first length on the edge's first face in the
# kernel's face order. OCCT's first is the bottom: the wedge it cuts is put here explicitly.)
chamfer_wedge = cq.Workplane("YZ", origin=(-WALKING_FOOT_WIDTH / 2, 0, 0))\
    .polyline([(-WALKING_FOOT_HEIGHT / 2, 0),
               (-WALKING_FOOT_HEIGHT / 2 + WALKING_FOOT_HEIGHT / 4, 0),
               (-WALKING_FOOT_HEIGHT / 2, WALKING_FOOT_DEPTH / 3)])\
    .close()\
    .extrude(WALKING_FOOT_WIDTH)

walking_foot = walking_foot.cut(chamfer_wedge)\
    .edges("<Y")\
    .edges("<Z")\
    .fillet(40)\
    .edges("<Z")\
    .edges("<Y")\
    .fillet(1)

walking_foot_block = walking_foot_block.cut(walking_foot)\
    .edges(cq.NearestToPointSelector([4,0,gridfinity.block_top_surface(3)]))\
    .fillet(1)\
    .edges(cq.NearestToPointSelector([0,4,gridfinity.block_top_surface(3)]))\
    .fillet(1)\
    .edges(cq.NearestToPointSelector([-4,0,gridfinity.block_top_surface(3)]))\
    .fillet(1)\
    .edges(cq.NearestToPointSelector([0,-4,gridfinity.block_top_surface(3)]))\
    .fillet(1)

del walking_foot
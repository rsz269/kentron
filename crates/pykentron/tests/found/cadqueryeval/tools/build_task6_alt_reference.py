"""Build reference/task6_alt.stl, the alternate reading of task6.

The hole sits 5mm from a long edge, centered along the length. The plate uses
the same frame as reference/task6.stl: x 0..40, y 0..20, z 0..5.

Run with:
    uv run --no-project --python 3.12 --with cadquery \\
        python tools/build_task6_alt_reference.py
"""

import cadquery as cq

OUTPUT = "src/cadqueryeval/data/reference/task6_alt.stl"

result = (
    cq.Workplane("XY")
    .box(40, 20, 5, centered=False)
    .faces(">Z")
    .workplane(origin=(0, 0, 5))
    .center(20, 5)
    .cboreHole(4, 8, 2)
)
cq.exporters.export(result, OUTPUT, tolerance=0.01, angularTolerance=0.1)

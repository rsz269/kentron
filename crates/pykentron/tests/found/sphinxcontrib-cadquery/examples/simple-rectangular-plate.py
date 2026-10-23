"""Simple rectangular plate."""

import cadquery

result = cadquery.Workplane().box(2, 2, 0.5)

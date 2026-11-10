import cadquery as cq
from typing import Sequence, Tuple, Literal

PlaneName = Literal["XY", "XZ", "YZ"]








def create_profile_from_straight_connections(
    l: Sequence[Tuple[float, float]],
    plane: PlaneName = "XY",
    closed: bool = False
) -> cq.Workplane:
    s = cq.Sketch()
    for i in range(len(l) - 1):
        s = s.segment(l[i], l[i + 1])
    if closed:
        s = s.close().assemble()
    wp = cq.Workplane(plane).placeSketch(s)
    wp._plane_name = plane  # store it manually  # type: ignore[attr-defined]
    return wp




















# def create_profile_from_straight_connections(
#     l: Sequence[Tuple[float,float]],
#     plane: PlaneName = "XY",
#     closed: bool = False
#     ) -> cq.Workplane:
#     """
#     l : list -> list of points to be connected through straight line
#     plane : str -> plane for the workplane
#     closed: False -> by default, do not connect the first and last points
#     """
#     wp = cq.Workplane(plane).polyline(l)
#     return wp.close() if closed else wp

# def create_profile_from_straight_connections(
#     l: Sequence[Tuple[float, float]],
#     plane: PlaneName = "XY",
# ) -> cq.Workplane:
#     """
#     l : list of (x, y) points connected by straight lines (auto-closed)
#     plane : workplane plane
#     Returns a Workplane with a closed face, ready for extrude/revolve.
#     """
#     sketch = cq.Sketch().polygon(l)
#     return cq.Workplane(plane).placeSketch(sketch)



# def create_profile_from_straight_connections(
#     l: Sequence[Tuple[float, float]],
#     plane: PlaneName = "XY",
#     closed: bool = False
# ) -> cq.Sketch:
#     """
#     l : list of points to be connected through straight line
#     plane : kept for API compatibility; Sketch itself is plane-agnostic
#     closed: False -> by default, do not connect the first and last points
#     """
#     s = cq.Sketch().segment(l)
#     return s.close() if closed else s



# def create_profile_from_straight_connections(
#     l: Sequence[Tuple[float, float]],
#     plane: PlaneName = "XY",
#     closed: bool = False
# ) -> cq.Sketch:
#     s = cq.Sketch()
#     for i in range(len(l) - 1):
#         s = s.segment(l[i], l[i + 1])
#     return s.close() if closed else s



# FUNCTIONA CON 12345
# def create_profile_from_straight_connections(
#     l: Sequence[Tuple[float, float]],
#     plane: PlaneName = "XY",
#     closed: bool = False
# ) -> cq.Workplane:
#     wp = cq.Workplane(plane).polyline(l)
#     return wp.close() if closed else wp

# def create_profile_from_straight_connections(
#     l: Sequence[Tuple[float, float]],
#     plane: PlaneName = "XY",
#     closed: bool = False
# ) -> cq.Sketch:
#     s = cq.Sketch()
#     for i in range(len(l) - 1):
#         s = s.segment(l[i], l[i + 1])
#     return s.close() if closed else s

"""
def create_profile_from_straight_connections(
    l: Sequence[Tuple[float, float]],
    plane: PlaneName = "XY",
    closed: bool = False
) -> Tuple[cq.Sketch, PlaneName]:
    s = cq.Sketch()
    for i in range(len(l) - 1):
        s = s.segment(l[i], l[i + 1])
    sketch = s.close() if closed else s
    return sketch, plane
"""





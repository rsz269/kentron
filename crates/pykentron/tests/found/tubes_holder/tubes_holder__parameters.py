import logging
import math as m
from importlib import reload

from typeguard import typechecked

log = logging.getLogger(__name__)

# ============
#  Tolerances
# ============


@typechecked
class Tol:
    """
    Tolerances constants for the model.
    All dimensions are in mm if not specified.
    """

    # General tolerance
    e = 0.1
    # Fit tolerance
    fit_e = 0.2


# ============
#  Dimensions
# ============


@typechecked
class TubeHolder:
    """
    Dimensions of the tube holder.
    All dimensions are in mm if not specified.
    """

    height = 60.0
    tube_height = 100.0
    tube_cap_diameter = 19.0
    tube_diameter = 15.0
    tube_hole_diameter = tube_diameter + 2.1
    tube_hole_radius = tube_hole_diameter / 2
    spacing_between_holes = 20.0
    nb_holes = 6
    width = 45
    length = (
        tube_hole_diameter * (nb_holes + 2) + (nb_holes - 1) * spacing_between_holes
    )
    wall_thickness = 3.0
    fillet_holes = wall_thickness / 2.0
    middle_height = 10 - wall_thickness / 2.0
    fillet_foot = 0.5

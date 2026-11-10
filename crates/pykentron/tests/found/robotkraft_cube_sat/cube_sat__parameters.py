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
    tol = 0.1


# ============
#  Dimensions
# ============


@typechecked
class CubeSat:
    """
    Dimensions of the cube sat.
    All dimensions are in mm if not specified.
    """

    pin_smooth_height = 1.5
    pin_smooth_disc_ratio = 0.2  # 20% the size of the disc above the flat cut pin to create the tip of the tooth like structure that will guide the piece to be assembled correctly by the robot. This is the size of the radius of the disc relative to the diameter of the bounding box of the face of the top of the pin when it was cut.

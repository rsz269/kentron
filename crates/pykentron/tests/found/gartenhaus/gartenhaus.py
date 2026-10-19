import json
import cadquery as cq
import math

# Constants
SPACING = 200
SAFETY_SPACING = 15
DEFAULT_WINDOW_WIDTH = 1095
DEFAULT_WINDOW_HEIGHT = 1525
DEFAULT_DOOR_WIDTH = 1000
DEFAULT_DOOR_HEIGHT = 2000
DEFAULT_JOIST_SPACING = 400
DEFAULT_FLOOR_LENGTH = 2000
DEFAULT_BEAM_LENGTH = 2000
DEFAULT_BEAM_WIDTH = 95
DEFAULT_BEAM_HEIGHT = 95
DEFAULT_BEAM_NOTCH_DEPTH = 30
DEFAULT_BEAM_NOTCH_WIDTH = 40
DEFAULT_BEAM_NOTCH_HEIGHT = 40


class CadQueryObject(cq.Workplane):
    def __init__(self, inPlane="XY", origin=(0, 0, 0), obj=None, spacing=(0, 0, 0)):
        super().__init__(inPlane=inPlane, origin=origin, obj=obj)
        self.spacing = spacing

    def add_spacing(self):
        obj = self.translate(self.spacing)
        return obj


class CadQueryObjectList(list):
    def add_spacing(self):
        return [obj.add_spacing() for obj in self]

    def union(self):
        result = self[0]
        for obj in self[1:]:
            result = result.union(obj)
        return result


def get_roof_construction(
    beam_length: float,
    beam_width: float,
    beam_height: float,
    cut_objects: CadQueryObjectList,
) -> CadQueryObjectList:
    distance = 485
    # Calculate the original x-position
    x_pos = ((beam_length) / 2) - (
        (beam_length - (DEFAULT_FLOOR_LENGTH + DEFAULT_BEAM_WIDTH)) / 2
    )

    # Adjust the x-position for rotation
    x_pos_after_rot = (
        x_pos - (x_pos - beam_length / 2) * (1 - math.cos(math.radians(8))) + 1
    )
    beam = (
        CadQueryObject("XY")
        .box(beam_length, beam_width, beam_height)
        .rotate((0, 0, 0), (0, 1, 0), 8)
        .translate(
            (
                x_pos_after_rot,
                1212.5,
                2300,
            )
        )
    )
    beams = [beam.translate((0, -distance * i, 0)) for i in range(6)]
    for cut_obj in cut_objects:
        for idx, beam in enumerate(beams):
            beams[idx] = beam.cut(cut_obj.union())
            beams[idx].tag(f"roof_beam_{idx}")
    for beam in beams:
        beam.spacing = (0, 0, SPACING * 4.9)
    return CadQueryObjectList(beams)


def get_front(
    floor_length: float,
    beam_length: float,
    beam_width: float,
    beam_height: float,
    beam_notch_depth: float,
    beam_notch_width: float,
    beam_notch_height: float,
    window_width: float = DEFAULT_WINDOW_WIDTH,
    window_height: float = DEFAULT_WINDOW_HEIGHT,
) -> CadQueryObjectList:
    beam = (
        get_join_beam(
            beam_length,
            beam_width,
            beam_height,
            beam_notch_depth,
            beam_notch_width,
            beam_notch_height,
        )
        .rotate((0, 0, 0), (0, 0, 1), -90)
        .translate(
            (
                floor_length / 2 + beam_width / 2,
                floor_length / 2 + beam_width / 2,
                +floor_length + beam_width,
            )
        )
    )
    window_beam_right = (
        get_join_beam(
            beam_length,
            beam_width,
            beam_height,
            beam_notch_depth,
            beam_notch_width,
            beam_notch_height,
        )
        .rotate((0, 0, 0), (1, 0, 0), -90)
        .translate(
            (
                floor_length / 2 - window_width / 2 - SAFETY_SPACING / 2,
                floor_length / 2 + beam_width / 2,
                beam_length / 2 + beam_width / 2,
            )
        )
    )
    window_beam_left = window_beam_right.translate(
        (beam_width + SAFETY_SPACING + window_width, 0, 0)
    )
    window_beam_bottom = (
        get_join_beam(
            window_width + SAFETY_SPACING,
            beam_width,
            beam_height,
            beam_notch_depth,
            beam_notch_width,
            beam_notch_height,
        )
        .rotate((0, 0, 0), (0, 0, 1), -90)
        .translate(
            (
                floor_length / 2 + beam_width / 2,
                floor_length / 2 + beam_width / 2,
                beam_length + beam_width - window_height - SAFETY_SPACING,
            )
        )
    )
    beam = beam.cut(window_beam_right).cut(window_beam_left)
    window_beam_right = window_beam_right.cut(window_beam_bottom)
    window_beam_left = window_beam_left.cut(window_beam_bottom)
    beam = beam.tag("front_beam")
    window_beam_right = window_beam_right.tag("window_beam_right")
    window_beam_left = window_beam_left.tag("window_beam_left")
    window_beam_bottom = window_beam_bottom.tag("window_beam_bottom")
    beam.spacing = (0, 0, SPACING * 2.9)
    window_beam_right.spacing = (-SPACING / 2, 0, SPACING)
    window_beam_left.spacing = (SPACING / 2, 0, SPACING)

    return CadQueryObjectList(
        [beam, window_beam_left, window_beam_right, window_beam_bottom]
    )


def get_cut(
    width: float,
    height: float,
    notch_depth: float,
    notch_width: float,
    notch_height: float,
) -> cq.Workplane:
    cut = (
        cq.Workplane("top")
        .rect(width, height)
        .rect(notch_width, notch_height)
        .extrude(notch_depth)
    )
    return cut


def get_join_beam(
    length: float,
    width: float,
    height: float,
    notch_depth: float,
    notch_width: float,
    notch_height: float,
) -> CadQueryObject:
    beam = CadQueryObject("YX").box(length + 2 * notch_depth, width, height)
    cut = get_cut(width, height, notch_depth, notch_width, notch_height)
    return beam.cut(cut.translate((0, length / 2, 0))).cut(
        cut.translate((0, -length / 2 - notch_depth, 0))
    )


def get_floor_beams(
    length: float,
    width: float,
    height: float,
    beam_length: float,
    beam_width: float,
    beam_height: float,
    beam_notch_depth: float,
    beam_notch_width: float,
    beam_notch_height: float,
    joist_spacing: float = DEFAULT_JOIST_SPACING,
) -> CadQueryObjectList:
    beam_count = int(width / joist_spacing) + 1
    spacing = width / (beam_count - 1) if beam_count > 1 else 0
    spacing = spacing - beam_width / (beam_count - 1)
    beams = CadQueryObjectList(
        [
            get_join_beam(
                beam_length,
                beam_width,
                beam_height,
                beam_notch_depth,
                beam_notch_width,
                beam_notch_height,
            ).translate((i * spacing, 0, 0))
            for i in range(beam_count)
        ]
    )
    return beams


def get_floor_construction(
    length: float,
    width: float,
    height: float,
    beam_length: float,
    beam_width: float,
    beam_height: float,
    beam_notch_depth: float,
    beam_notch_width: float,
    beam_notch_height: float,
    cut_objects: CadQueryObjectList,
    joist_spacing: float = DEFAULT_JOIST_SPACING,
) -> CadQueryObjectList:
    floor_beams = get_floor_beams(
        length,
        width + beam_width * 2,
        height,
        beam_length,
        beam_width,
        beam_height,
        beam_notch_depth,
        beam_notch_width,
        beam_notch_height,
        joist_spacing,
    )
    beam = (
        CadQueryObject("XY", spacing=(0, SPACING * 2, 0))
        .box(length + beam_width * 2, beam_width, beam_height)
        .translate((length / 2 + beam_width / 2, length / 2 + beam_height / 2, 0))
        .tag("floor_beam_main")
    )
    for cut_obj in cut_objects:
        cut_obj = cut_obj.union()
        for idx, obj in enumerate(floor_beams):
            obj = obj.cut(cut_obj)
            floor_beams[idx] = obj
            floor_beams[idx].tag(f"floor_beam_{idx}")
        beam = beam.cut(cut_obj.union())
    beam_2 = beam.translate((0, -width - beam_height, 0)).tag("floor_beam_main_2")
    beam_2.spacing = (0, -SPACING, 0)
    beam = beam.cut(floor_beams.union())
    beam_2 = beam_2.cut(floor_beams.union())
    beam.spacing = (0, SPACING * 2, 0)
    beam_2.spacing = (0, -SPACING, 0)
    return CadQueryObjectList(floor_beams + [beam, beam_2])


def get_right_construction(
    floor_length: float,
    beam_length: float,
    beam_width: float,
    beam_height: float,
    beam_notch_depth: float,
    beam_notch_width: float,
    beam_notch_height: float,
    cut_objects: CadQueryObjectList,
    door_width: float = DEFAULT_DOOR_WIDTH,
    door_height: float = DEFAULT_DOOR_HEIGHT,
) -> CadQueryObjectList:
    beam = get_join_beam(
        beam_length,
        beam_width,
        beam_height,
        beam_notch_depth,
        beam_notch_width,
        beam_notch_height,
    ).rotate((0, 0, 0), (1, 0, 0), -90)
    beam_top = (
        CadQueryObject("YX", spacing=(0, 0, SPACING * 5))
        .box(beam_length + 200, beam_width, beam_height)
        .translate((0, 0, beam_length + beam_width))
        .tag("right_beam_top")
    )
    beam_right = beam.translate(
        (0, floor_length / 2 + beam_width / 2, beam_length / 2 + beam_width / 2)
    ).tag("right_beam_right")
    beam_right.spacing = (0, 0, SPACING)
    beam_left = beam.translate(
        (0, -floor_length / 2 - beam_width / 2, beam_length / 2 + beam_width / 2)
    ).tag("right_beam_left")
    beam_left.spacing = (0, 0, SPACING)
    beam_door_right = beam_left.translate((0, 200, 0)).tag("right_beam_door_right")
    beam_door_left = beam_door_right.translate(
        (0, door_width + 2 * SAFETY_SPACING + beam_width, 0)
    ).tag("right_beam_door_left")
    beam_door_top = (
        get_join_beam(
            door_width + 2 * SAFETY_SPACING,
            beam_width,
            beam_height,
            beam_notch_depth,
            beam_notch_width,
            beam_notch_height,
        )
        .translate(
            (
                0,
                -floor_length / 2 + door_width / 2 + 200 + SAFETY_SPACING,
                door_height + beam_width + SAFETY_SPACING,
            )
        )
        .tag("right_beam_door_top")
    )
    beam_door_left = beam_door_left.cut(beam_door_top)
    beam_door_right = beam_door_right.cut(beam_door_top)
    beams = beam_right.union(beam_left).union(beam_door_left).union(beam_door_right)
    beam_top = beam_top.cut(beams)
    for cut_obj in cut_objects:
        beam_right = beam_right.cut(cut_obj.union())
        beam_left = beam_left.cut(cut_obj.union())
    beam_top.spacing = (0, 0, SPACING * 2.8)
    beam_door_left.spacing = (0, SPACING, SPACING)
    beam_door_right.spacing = (0, 0, SPACING)
    beam_door_top.spacing = (0, SPACING / 2, SPACING)
    beam_right.spacing = (0, 0, SPACING)
    beam_left.spacing = (0, 0, SPACING)
    return CadQueryObjectList(
        [
            beam_right,
            beam_left,
            beam_door_left,
            beam_door_right,
            beam_top,
            beam_door_top,
        ]
    )


def get_left_construction(
    floor_length: float,
    beam_length: float,
    beam_width: float,
    beam_height: float,
    beam_notch_depth: float,
    beam_notch_width: float,
    beam_notch_height: float,
    cut_objects: CadQueryObjectList,
) -> CadQueryObjectList:
    beam = (
        get_join_beam(
            beam_length,
            beam_width,
            beam_height,
            beam_notch_depth,
            beam_notch_width,
            beam_notch_height,
        )
        .rotate((0, 0, 0), (1, 0, 0), -90)
        .translate(
            (
                floor_length + beam_width,
                -floor_length / 2 - beam_width / 2,
                beam_length / 2 + beam_width / 2,
            )
        )
        .tag("left_beam_main")
    )
    beam_count = int((floor_length + beam_width * 2) / 500) + 1
    beams = [beam]
    for idx in range(1, beam_count):
        beams.append(
            beam.translate((0, idx * 500 + beam_width, 0)).tag(f"left_beam_{idx}")
        )
    beam_top = (
        CadQueryObject("YX", spacing=(0, 0, SPACING * 5))
        .box(floor_length + 200 + 300, beam_width, beam_height)
        .translate((floor_length + beam_width, 0, beam_length + beam_width))
        .tag("left_beam_top")
    )
    beam_top = beam_top.cut(CadQueryObjectList(beams).union())
    for cut_obj in cut_objects:
        beam_top = beam_top.cut(cut_obj.union())
    for beam in beams:
        beam.spacing = (0, 0, SPACING)
    beam_top.spacing = (0, 0, SPACING * 2)
    return CadQueryObjectList(beams + [beam_top])


def standardize_beam_orientation(beam: CadQueryObject) -> CadQueryObject:
    # Get the bounding box of the beam
    bounding_box = beam.findSolid().BoundingBox()

    # Determine the orientation based on the bounding box dimensions
    if bounding_box.xlen > bounding_box.ylen and bounding_box.xlen > bounding_box.zlen:
        # If the x dimension is the largest, no rotation is needed
        pass
    elif (
        bounding_box.ylen > bounding_box.xlen and bounding_box.ylen > bounding_box.zlen
    ):
        # If the y dimension is the largest, rotate around the z-axis to align along the x-axis
        beam = beam.rotate((0, 0, 0), (0, 0, 1), 90)
    elif (
        bounding_box.zlen > bounding_box.xlen and bounding_box.zlen > bounding_box.ylen
    ):
        # If the z dimension is the largest, rotate around the y-axis to align along the x-axis
        beam = beam.rotate((0, 0, 0), (0, 1, 0), -90)
    return beam


def arrange_beams_in_row(
    beams: CadQueryObjectList, spacing: float = 200
) -> CadQueryObjectList:
    arranged_beams = CadQueryObjectList()
    for i, beam in enumerate(beams):
        standardized_beam = standardize_beam_orientation(beam)
        bbox = standardized_beam.findSolid().BoundingBox()
        standardized_beam = standardized_beam.translate(
            (-bbox.center.x, -bbox.center.y, -bbox.center.z)
        )
        standardized_beam = standardized_beam.translate((0, i * spacing, 0))
        bbox = standardized_beam.findSolid().BoundingBox()
        zdiff = round(bbox.zmax, 2) - round(bbox.zmin, 2)
        if zdiff > 100:  # re rotate if the beam is too high
            standardized_beam = standardized_beam.rotate((0, 0, 0), (0, 1, 0), -8)
        arranged_beams.append(standardized_beam)

    return arranged_beams


def get_inventory(parts: CadQueryObjectList) -> dict:
    inventory = {}
    for part in parts:
        bbox = part.findSolid().BoundingBox()
        dimensions = (round(bbox.xlen, 2), round(bbox.ylen, 2), round(bbox.zlen, 2))
        dim_str = f"{dimensions[0]}x{dimensions[1]}x{dimensions[2]}"
        inventory[dim_str] = inventory.get(dim_str, 0) + 1
    return inventory


def main():
    front_obj = get_front(
        floor_length=DEFAULT_FLOOR_LENGTH,
        beam_length=DEFAULT_BEAM_LENGTH,
        beam_width=DEFAULT_BEAM_WIDTH,
        beam_height=DEFAULT_BEAM_HEIGHT,
        beam_notch_depth=DEFAULT_BEAM_NOTCH_DEPTH,
        beam_notch_width=DEFAULT_BEAM_NOTCH_WIDTH,
        beam_notch_height=DEFAULT_BEAM_NOTCH_HEIGHT,
    )

    back_obj = get_front(
        floor_length=DEFAULT_FLOOR_LENGTH,
        beam_length=DEFAULT_BEAM_LENGTH,
        beam_width=DEFAULT_BEAM_WIDTH,
        beam_height=DEFAULT_BEAM_HEIGHT,
        beam_notch_depth=DEFAULT_BEAM_NOTCH_DEPTH,
        beam_notch_width=DEFAULT_BEAM_NOTCH_WIDTH,
        beam_notch_height=DEFAULT_BEAM_NOTCH_HEIGHT,
    )
    for idx, obj in enumerate(back_obj):
        back_obj[idx] = obj.translate(
            (0, -DEFAULT_FLOOR_LENGTH - DEFAULT_BEAM_WIDTH, 0)
        )
        back_obj[idx].spacing = {
            0: (0, -SPACING * 2, SPACING * 2.9),
            1: (SPACING, -SPACING * 2, SPACING),
            2: (-SPACING / 2, -SPACING * 2, SPACING),
            3: (SPACING / 2, -SPACING * 2, SPACING),
        }[idx]

    right_obj = get_right_construction(
        floor_length=DEFAULT_FLOOR_LENGTH,
        beam_length=2300,
        beam_width=DEFAULT_BEAM_WIDTH,
        beam_height=DEFAULT_BEAM_HEIGHT,
        beam_notch_depth=DEFAULT_BEAM_NOTCH_DEPTH,
        beam_notch_width=DEFAULT_BEAM_NOTCH_WIDTH,
        beam_notch_height=DEFAULT_BEAM_NOTCH_HEIGHT,
        cut_objects=[front_obj, back_obj],
    )
    left_obj = get_left_construction(
        floor_length=DEFAULT_FLOOR_LENGTH,
        beam_length=DEFAULT_BEAM_LENGTH,
        beam_width=DEFAULT_BEAM_WIDTH,
        beam_height=DEFAULT_BEAM_HEIGHT,
        beam_notch_depth=DEFAULT_BEAM_NOTCH_DEPTH,
        beam_notch_width=DEFAULT_BEAM_NOTCH_WIDTH,
        beam_notch_height=DEFAULT_BEAM_NOTCH_HEIGHT,
        cut_objects=[front_obj, back_obj],
    )

    floor_obj = get_floor_construction(
        length=DEFAULT_FLOOR_LENGTH,
        width=DEFAULT_FLOOR_LENGTH,
        height=DEFAULT_BEAM_HEIGHT,
        beam_length=DEFAULT_BEAM_LENGTH,
        beam_width=DEFAULT_BEAM_WIDTH,
        beam_height=DEFAULT_BEAM_HEIGHT,
        beam_notch_depth=DEFAULT_BEAM_NOTCH_DEPTH,
        beam_notch_width=DEFAULT_BEAM_NOTCH_WIDTH,
        beam_notch_height=DEFAULT_BEAM_NOTCH_HEIGHT,
        cut_objects=[right_obj, left_obj, front_obj],
    )
    roof_obj = get_roof_construction(
        beam_length=3000,
        beam_width=75,
        beam_height=75,
        cut_objects=[right_obj, left_obj, front_obj, back_obj],
    )

    all_objects = CadQueryObjectList(
        right_obj + floor_obj + left_obj + front_obj + back_obj + roof_obj
    )

    right_obj_spaced = right_obj.add_spacing()
    floor_obj_spaced = floor_obj.add_spacing()
    left_obj_spaced = left_obj.add_spacing()
    front_obj_spaced = front_obj.add_spacing()
    back_obj_spaced = back_obj.add_spacing()
    roof_obj_spaced = roof_obj.add_spacing()

    exploded_garden_house = CadQueryObjectList(
        right_obj_spaced
        + floor_obj_spaced
        + left_obj_spaced
        + front_obj_spaced
        + back_obj_spaced
        + roof_obj_spaced
    ).union()
    garden_house = CadQueryObjectList(
        right_obj + floor_obj + left_obj + front_obj + back_obj + roof_obj
    ).union()
    parts = arrange_beams_in_row(all_objects, spacing=DEFAULT_BEAM_WIDTH + 100)
    inventory = get_inventory(parts)
    with open("gartenhaus_inventory.json", "w") as f:
        json.dump(inventory, f, indent=4)
    cq.exporters.export(exploded_garden_house, "gartenhaus_exploded.step")
    cq.exporters.export(garden_house, "gartenhaus.step")
    cq.exporters.export(parts.union(), "gartenhaus_parts.step")

    show_object(garden_house)


main()

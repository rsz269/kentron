import math

import cadquery as cq


def get_center(p1, p2):
    p1x, p1y = p1
    p2x, p2y = p2

    mx = (p1x + p2x) / 2
    my = (p1y + p2y) / 2

    a = -(p1x - p2x) / (p1y - p2y)
    b = my - a * mx

    return b


def get_bottom(dent, diameter, height=0.0):
    r = -get_center((0, dent), (diameter / 2, height)) + dent
    result = (
        cq.Workplane("XY")
        .sphere(r)
        .translate((0, 0, -r + dent))
        .intersect(
            cq.Workplane("XY")
            .cylinder(2 * dent, (diameter / 2))
            .translate((0, 0, dent))
        )
    )

    return result


def knurl(obj, height, radius, cut_angle, cut_depth, angle, n, point_radius=None):
    length = 2 * cut_depth
    x1 = length * math.cos(math.radians(cut_angle / 2))
    y1 = length * math.sin(math.radians(cut_angle / 2))

    def doprofile(loc):
        s = (
            cq.Sketch()
            .push([(-cut_depth, 0)])
            .segment((0, 0), (x1, y1))
            .segment((x1, -y1))
            .close()
            .assemble()
            .reset()
        )
        if point_radius:
            s = s.vertices("<X").fillet(point_radius).reset()
        return s.moved(loc)

    s1 = cq.Sketch().parray(radius, 0, 360, n).each(doprofile)

    w = cq.Workplane()
    w1 = w.placeSketch(s1).twistExtrude(height, angle)
    w2 = w.placeSketch(s1).twistExtrude(height, -angle)
    tool = w1.union(w2, clean=True, tol=1e-2)

    return obj.cut(tool)


def cut_in_half(obj):
    return (
        obj.faces("<Z")
        .workplane()
        .transformed(rotate=cq.Vector(-90, 0, 0))
        .split(keepBottom=True)
    )

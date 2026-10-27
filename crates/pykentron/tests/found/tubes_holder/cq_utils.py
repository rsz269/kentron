import math as m
import cadquery as cq
from cadquery.selectors import *
from cadquery import Workplane, Sketch, Vector, Location

"""
from OCC.Core.gp import (gp_Mat, gp_Trsf, gp_GTrsf)
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_GTransform,
                                     BRepBuilderAPI_Transform)

def extruder( face ):
    f2ext = face.wires()
    f2ext.ctx.pendingWires = f2ext.vals()
    return f2ext

def scaleXYZ( obj, sX, sY, sZ ):
    xform = gp_GTrsf()
    xform.SetVectorialPart(gp_Mat(
    sX, 0, 0,
    0, sY, 0,
    0, 0, sZ,
    ))
    
    brep = BRepBuilderAPI_GTransform(obj.wrapped, xform, False)
    brep.Build()
    fshape = brep.Shape()
    return cq.Shape.cast(fshape)
"""


def argcmp_f(l, cmp=lambda x, y: x > y, f=lambda x: x):
    i = 0
    v = None
    item = None
    argmax = None
    # Test if the collection is empty
    if not l:
        return None

    for e in l:
        w = f(e)
        if (None == v) or cmp(w, v):
            v = w
            item = e
            argmax = i
        i += 1
    return (argmax, v, item)


def argmax_f(l, f=lambda x: x):
    return argcmp_f(l, f=f)


def argmin_f(l, f=lambda x: x):
    return argcmp_f(l, cmp=lambda x, y: x < y, f=f)


def argmax_cqVector(v):
    return argmax_f(v.toTuple(), lambda e: abs(e))


def argmin_cqVector(v):
    return argmin_f(v.toTuple(), lambda e: abs(e))


def orthonormalVector(V):
    basis = [cq.Vector((1, 0, 0)), cq.Vector((0, 1, 0)), cq.Vector((0, 0, 1))]
    n = cq.Vector(V)
    if cq.Vector((0, 0, 0)) == n:
        raise ValueError("Input vector is null")

    (idx, _, _) = argmin_cqVector(n)

    e = basis[idx]
    return n.cross(e).normalized()


def ABcylinder(A, B, radius):
    l = cq.Vector(B - A)
    n = l.normalized()
    xDir = orthonormalVector(n)
    pl = cq.Plane(A, xDir, n)
    return cq.Workplane(pl).circle(radius).extrude(l.Length)

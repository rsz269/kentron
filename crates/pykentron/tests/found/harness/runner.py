"""runner.py MODE SCRIPT PROJECT OUT.json: runs a CadQuery script (as __main__) on pykentron
(MODE kentron) or CadQuery (MODE cq) and writes what it made (shown objects and module-level
CadQuery objects: volume, area, box, validity, counts) and how long it took."""
import builtins, json, os, runpy, sys, time, traceback

mode, script, project, out = sys.argv[1:5]
KENTRON_PY = os.path.expanduser('~/prj/kentron/crates/pykentron/python')
if mode == 'kentron':
    sys.path[:0] = [os.path.join(KENTRON_PY, 'pykentron', 'cadquery_shim'), KENTRON_PY]
    os.environ['CAD_KERNEL'] = 'kentron'
else:
    os.environ['CAD_KERNEL'] = 'occt'
# (The project's own modules after everything else: they must not shadow standard ones.)
for d in [os.path.dirname(os.path.abspath(script)), project, os.path.join(project, 'src')]:
    if d not in sys.path:
        sys.path.append(d)

import cadquery as cq

shown = []
# (Unnamed objects numbered among themselves: exports, recorded here too, are not counted, as
# a kernel may record one export twice through the methods it calls.)
unnamed = [0]
def show_object(obj, name=None, options=None, **kw):
    if name is None:
        name = f'shown{unnamed[0]}'
        unnamed[0] += 1
    shown.append((name, obj))
builtins.show_object = show_object
# (What viewers and CQ-editor provide: `show` collects as `show_object`, the rest nothing.)
builtins.show = lambda *objs, **kw: [show_object(o) for o in objs]
builtins.show_all = lambda *a, **k: None
builtins.set_defaults = lambda *a, **k: None
# (And the VS Code viewer's module, which scripts import to show what they make.)
import types
viewer = types.ModuleType('ocp_vscode')
viewer.show = builtins.show
viewer.show_object = show_object
viewer.__getattr__ = lambda name: (lambda *a, **k: None)
sys.modules.setdefault('ocp_vscode', viewer)
builtins.cq = cq
# (What a script exports counts as shown, named by its file.)
def recording(f, at=0):
    def wrapped(*args, **kw):
        try:
            obj = args[at] if len(args) > at else kw.get('w', kw.get('shape'))
            path = args[at + 1] if len(args) > at + 1 else kw.get('fname', kw.get('path', kw.get('fileName', 'export')))
            show_object(obj, name='export:' + os.path.basename(str(path)))
        except Exception:
            pass
        return f(*args, **kw)
    return wrapped
try:
    cq.exporters.export = recording(cq.exporters.export)
except Exception:
    pass
for cls, names in ((cq.Shape, ('export', 'exportStl', 'exportStep', 'exportBrep')), (cq.Assembly, ('save', 'export')), (cq.Workplane, ('export',))):
    for n in names:
        if hasattr(cls, n):
            setattr(cls, n, recording(getattr(cls, n)))
builtins.debug = lambda *a, **k: None
builtins.log = lambda *a, **k: None

def describe(obj):
    if isinstance(obj, cq.Assembly):
        shape = obj.toCompound()
    elif isinstance(obj, cq.Workplane):
        vals = [v for v in obj.vals() if isinstance(v, cq.Shape)]
        if not vals:
            return None
        shape = vals[0] if len(vals) == 1 else cq.Compound.makeCompound(vals)
    elif isinstance(obj, cq.Shape):
        shape = obj
    else:
        return None
    # (The tight box: OCCT's default one is grown by the shapes' tolerances and its curves'
    # control polygons.)
    if mode == 'kentron':
        bb = shape.BoundingBox()
        box = [bb.xmin, bb.ymin, bb.zmin, bb.xmax, bb.ymax, bb.zmax]
    else:
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        b = Bnd_Box()
        BRepBndLib.AddOptimal_s(shape.wrapped, b, False, False)
        box = list(b.Get())
    d = {'type': shape.ShapeType(), 'solids': len(shape.Solids()), 'faces': len(shape.Faces()), 'bbox': box}
    # (The volume of the solids only: OCCT gives wires and faces in a compound volumes of
    # their own.)
    volume = lambda: sum(s.Volume() for s in shape.Solids())
    for k, f in (('volume', volume), ('area', shape.Area), ('valid', shape.isValid)):
        try:
            d[k] = f()
        except Exception as e:
            d[k] = f'error: {e}'
    return d

result = {'mode': mode, 'script': script}
t, c = time.time(), time.process_time()
try:
    g = runpy.run_path(script, run_name='__main__')
    result['status'] = 'ok'
except SystemExit:
    g = {}
    result['status'] = 'ok'
except BaseException as e:
    g = {}
    result['status'] = 'error'
    result['error'] = f'{type(e).__name__}: {e}'[:400]
    result['trace'] = traceback.format_exc()[-1500:]
result['wall'] = time.time() - t
result['cpu'] = time.process_time() - c
objects = {}
for name, obj in shown + [(k, v) for k, v in g.items() if not k.startswith('_')]:
    try:
        d = describe(obj)
    except Exception as e:
        d = {'error': str(e)[:200]}
    if d is not None and name not in objects:
        objects[name] = d
result['objects'] = objects
with open(out, 'w') as f:
    json.dump(result, f, default=str)

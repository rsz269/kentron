"""hunt.py PROJECT_DIR [GLOB...]: up to 30 of the project's CadQuery scripts, each run on
pykentron and on CadQuery (fresh processes, 120 s limit), compared; results appended to
results/<project>.jsonl and a line per script printed."""
import json, os, re, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
HUNT = HERE.parent
PY = str(HUNT / 'venv' / 'bin' / 'python')
LIMIT = 120

def candidates(root: Path, globs):
    files = []
    pats = globs or ['**/*.py']
    for g in pats:
        files += [p for p in root.glob(g) if p.is_file()]
    skip = re.compile(r'(^|/)(test_|tests?/|setup\.py|conftest|docs?/conf|__init__\.py|\.github|venv|site-packages|build/)')
    out = []
    for p in sorted(set(files)):
        rel = str(p.relative_to(root))
        if skip.search(rel):
            continue
        try:
            text = p.read_text(errors='replace')
        except Exception:
            continue
        if re.search(r'^\s*(import cadquery|from cadquery|import cq_\w+|from cq_\w+)', text, re.M) and ('show_object' in text or re.search(r'^\w+\s*=\s*\(?\s*cq\.|^result\b', text, re.M) or 'cq.Workplane' in text):
            out.append(p)
    # (Examples first.)
    out.sort(key=lambda p: (0 if re.search(r'example|demo|sample', str(p), re.I) else 1, str(p)))
    return out[:30]

def run(mode, script, project, tag):
    out = HUNT / 'results' / 'raw' / f'{tag}.{mode}.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    work = HUNT / 'work' / f'{tag}.{mode}'
    work.mkdir(parents=True, exist_ok=True)
    # (Where scripts commonly export to: CadQuery writes nothing to a missing directory,
    # quietly.)
    for d in ('stl', 'step', 'out', 'output', 'export', 'exports', 'gltf', 'svg', 'dxf', 'docs'):
        (work / d).mkdir(exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', PYKENTRON_PROFILE='fast', MPLBACKEND='Agg')
    t = time.time()
    try:
        p = subprocess.run([PY, '-I', str(HERE / 'runner.py'), mode, str(script), str(project), str(out)], cwd=work, env=env, capture_output=True, text=True, timeout=LIMIT)
        if not out.exists():
            return {'status': 'crash', 'error': (p.stderr or p.stdout)[-400:], 'wall': time.time() - t, 'objects': {}}
    except subprocess.TimeoutExpired:
        return {'status': 'timeout', 'wall': LIMIT, 'objects': {}}
    return json.loads(out.read_text())

def compare(k, c):
    """Differences of kentron's objects from CadQuery's."""
    diffs = []
    if c['status'] != 'ok':
        return ['cadquery failed']
    if k['status'] != 'ok':
        return [f"kentron {k['status']}: {k.get('error', '')[:200]}"]
    # (Exports to temporary files are named at random (`tempfile`'s tmpXXXXXXXX): matched by
    # their place among them.)
    temp = re.compile(r'^(export:.*?)tmp[a-z0-9_]{8}(\.\w+)$')
    def untemp(objects):
        out, seen = {}, {}
        for name, o in objects.items():
            m = temp.match(name)
            if m:
                stem = f'{m.group(1)}tmp*{m.group(2)}'
                seen[stem] = seen.get(stem, 0) + 1
                name = f'{stem}#{seen[stem]}'
            out[name] = o
        return out
    k, c = dict(k, objects=untemp(k['objects'])), dict(c, objects=untemp(c['objects']))
    for name, co in c['objects'].items():
        ko = k['objects'].get(name)
        if ko is None:
            diffs.append(f'{name}: missing')
            continue
        for key in ('volume', 'area'):
            # (Volumes of solids only, areas of faces only: OCCT's of wires mean nothing. Nor
            # of a shape OCCT made invalid where kentron's is valid: no reference.)
            if (key == 'volume' and not co.get('solids')) or (key == 'area' and not co.get('faces')):
                continue
            if co.get('valid') is False and ko.get('valid') is True:
                continue
            a, b = ko.get(key), co.get(key)
            if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                if abs(a - b) > 1e-3 * max(abs(b), 1e-6) + 1e-6:
                    diffs.append(f'{name}: {key} {a:.6g} vs {b:.6g}')
            elif a != b:
                diffs.append(f'{name}: {key} {a} vs {b}')
        if ko.get('valid') is not True and co.get('valid') is True:
            diffs.append(f'{name}: invalid')
        if 'bbox' in ko and 'bbox' in co:
            size = max(1e-6, max(co['bbox'][3] - co['bbox'][0], co['bbox'][4] - co['bbox'][1], co['bbox'][5] - co['bbox'][2]))
            # (CadQuery's boxes are OCCT's, enlarged by tolerances and control points: ours
            # may lie inside one, short of it by a little.)
            within = all(k >= c - 1e-2 * size for k, c in zip(ko['bbox'][:3], co['bbox'][:3])) and all(k <= c + 1e-2 * size for k, c in zip(ko['bbox'][3:], co['bbox'][3:]))
            if max(abs(x - y) for x, y in zip(ko['bbox'], co['bbox'])) > (3e-2 if within else 1e-2) * size:
                diffs.append(f'{name}: bbox differs')
        if ko.get('solids') != co.get('solids') and co.get('solids'):
            diffs.append(f"{name}: solids {ko.get('solids')} vs {co.get('solids')}")
    return diffs

def main():
    project = Path(sys.argv[1]).resolve()
    name = project.name
    rows = []
    for script in candidates(project, sys.argv[2:]):
        tag = f'{name}__' + str(script.relative_to(project)).replace('/', '__')
        k = run('kentron', script, project, tag)
        c = run('cq', script, project, tag)
        diffs = compare(k, c)
        slower = c['status'] == 'ok' and k['status'] == 'ok' and k['wall'] > 2 * c['wall'] and k['wall'] - c['wall'] > 0.5
        row = {'project': name, 'script': str(script.relative_to(project)), 'kentron': {x: k.get(x) for x in ('status', 'wall', 'cpu', 'error')}, 'cq': {x: c.get(x) for x in ('status', 'wall', 'cpu', 'error')}, 'objects': len(c.get('objects', {})), 'diffs': diffs, 'slower': slower}
        rows.append(row)
        flag = 'CQFAIL' if c['status'] != 'ok' else ('DIFF' if diffs else ('SLOW' if slower else 'ok'))
        print(f"{flag:6s} k {k.get('wall', 0):6.2f}s cq {c.get('wall', 0):6.2f}s objs {row['objects']:2d} {row['script']}  {'; '.join(diffs)[:160]}", flush=True)
    with open(HUNT / 'results' / f'{name}.jsonl', 'w') as f:
        for r in rows:
            f.write(json.dumps(r) + '\n')

main()

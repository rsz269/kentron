"""analyze.py: the flagged cases of the first 20 eligible projects (in the order run), each
with what went wrong and the suspected area, as results/flagged.json and a summary."""
import json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

HUNT = Path(__file__).resolve().parent.parent
ORDER = ['cadquery-repo', 'cq-gridfinity', 'oasis', 'alphabets', 'cad-as-code-book', 'cq_warehouse', 'cadquery-contrib', 'cq_gears', 'cq-electronics', 'cq-kit',
         'cqterrain', 'cadqueryhelper', 'cadquery-plugins', 'cqMore', 'cadquery-projects', 'CadQuery-Examples', 'cadquery-models', 'spindle-assy-example', 'partcad',
         'cq_keyboard_gen', 'parametric-keycaps', 'cq-uav', 'streamlit-wing-design', 'cq-cam', 'cqparts_bucket',
         'openvmp-models', 'dollhouse', 'QSC', 'vibe-cading', 'cq-queryabolt', 'STEMFIE-X-CadQuery-Edition', 'pcb-jigify', 'rover', 'octahedroflake', 'cqindustry']

# Suspected area from an error message, else from what the script does.
ERRORS = [
    (r'not provided by pykentron|has no attribute|OCP\.|OCC', 'API coverage (OCP / missing methods)'),
    (r'fillet|chamfer|rolling-ball', 'fillets and chamfers'),
    (r'sweep|loft|ThruSections', 'sweeps and lofts'),
    (r'spline|approximat|filling|N-sided|interp', 'splines and surface fitting'),
    (r'text|font', 'text'),
    (r'offset|shell|thick', 'offsets and shells'),
    (r'boolean|could not be built|union|cut|intersect', 'booleans'),
    (r'export|import|STL|STEP|No such file', 'import/export'),
    (r'timeout', 'performance'),
]
OPS = [(r'\.text\(', 'text'), (r'\.sweep\(|sweep\(', 'sweeps and lofts'), (r'\.loft\(', 'sweeps and lofts'), (r'\.fillet\(|\.chamfer\(', 'fillets and chamfers'),
       (r'\.shell\(|offset', 'offsets and shells'), (r'spline|interpPlate|makeSplineApprox', 'splines and surface fitting'), (r'helix|Helix', 'sweeps and lofts'),
       (r'\.union\(|\.cut\(|\.intersect\(|combine', 'booleans')]

def area(row, script_text):
    text = ' '.join(row['diffs']) + ' ' + (row['kentron'].get('error') or '')
    if row['kentron']['status'] == 'timeout':
        return 'performance'
    if row['kentron']['status'] != 'ok':
        for pat, a in ERRORS:
            if re.search(pat, text, re.I):
                return a
        return 'other errors'
    if 'solids' in text and not re.search(r'volume|area', text):
        return 'compound/solid structure'
    found = [a for pat, a in OPS if re.search(pat, script_text)]
    if row['slower'] and not row['diffs']:
        return 'performance: ' + (found[0] if found else 'general')
    return found[0] if found else 'other geometry'

def main():
    eligible = []
    flagged = []
    for name in ORDER:
        f = HUNT / 'results' / f'{name}.jsonl'
        if not f.exists():
            continue
        rows = [json.loads(l) for l in open(f)]
        good = [r for r in rows if r['cq']['status'] == 'ok' and r['objects'] > 0]
        if not good:
            continue
        eligible.append(name)
        for r in good:
            if not (r['diffs'] or r['slower']):
                continue
            path = HUNT / 'repos' / name / r['script']
            text = path.read_text(errors='replace') if path.exists() else ''
            # (Unseeded randomness: no two runs make the same thing.)
            if (re.search(r'\brandom\b|randint|uniform\(', text) and not re.search(r'seed\(', text)) or 'random' in r['script'].lower():
                continue
            r['area'] = area(r, text)
            r['kind'] = 'fails' if r['kentron']['status'] != 'ok' else ('diverges' if r['diffs'] else 'slower')
            flagged.append(r)
        if len(eligible) == 20:
            break
    json.dump({'eligible': eligible, 'flagged': flagged}, open(HUNT / 'results' / 'flagged.json', 'w'), indent=1)
    print('eligible projects:', len(eligible), eligible)
    print('flagged cases:', len(flagged), Counter(r['kind'] for r in flagged))
    by = defaultdict(list)
    for r in flagged:
        by[r['area']].append(r)
    for a, rs in sorted(by.items(), key=lambda x: -len(x[1])):
        print(f'\n== {a}: {len(rs)}')
        for r in rs[:40]:
            k, c = r['kentron'], r['cq']
            what = (k.get('error') or '; '.join(r['diffs']))[:150] if r['kind'] != 'slower' else f"{k['wall']:.2f}s vs {c['wall']:.2f}s"
            print(f"  {r['kind']:8s} {r['project']}/{r['script']}: {what}")

main()

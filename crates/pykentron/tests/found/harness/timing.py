"""timing.py [FILTER...]: the passing scripts copied here (`passing.tsv`) run again on
pykentron, each in a fresh process: wall and CPU time against CadQuery's recorded ones,
and the totals. A script whose project or path contains a FILTER only, where any are given.
`--repeat N` takes the least of N runs (as noisy as the machine is busy)."""
import argparse, json, os, subprocess, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOUND = HERE.parent
ap = argparse.ArgumentParser()
ap.add_argument('filters', nargs='*')
ap.add_argument('--repeat', type=int, default=1)
ap.add_argument('--python', default=sys.executable, help='the interpreter to run the scripts with')
args = ap.parse_args()

rows = [line.rstrip('\n').split('\t') for line in open(FOUND / 'passing.tsv')][1:]
total = {'k': 0.0, 'kc': 0.0, 'c': 0.0, 'cc': 0.0}
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', PYKENTRON_PROFILE='fast', MPLBACKEND='Agg')
print(f"{'kentron s':>9} {'cpu s':>7} {'cadquery s':>10} {'cpu s':>7} {'ratio':>6}  script")
for project, script, source, licence, k_wall, k_cpu, c_wall, c_cpu, objects in rows:
    key = f'{project}/{script}'
    if args.filters and not any(f in key for f in args.filters):
        continue
    path = FOUND / project / script
    if not path.exists():
        continue
    best = None
    for _ in range(args.repeat):
        with tempfile.TemporaryDirectory() as work:
            out = Path(work) / 'out.json'
            # (Where scripts commonly export to.)
            for d in ('stl', 'step', 'out', 'output', 'export', 'exports', 'gltf', 'svg', 'dxf'):
                (Path(work) / d).mkdir()
            subprocess.run([args.python, '-I', str(HERE / 'runner.py'), 'kentron', str(path), str(FOUND / project), str(out)], cwd=work, env=env, capture_output=True, timeout=600)
            if not out.exists():
                break
            r = json.loads(out.read_text())
            if r['status'] != 'ok':
                break
            if best is None or r['wall'] < best['wall']:
                best = r
    if best is None:
        print(f"{'failed':>9} {'':>7} {float(c_wall):10.2f} {float(c_cpu):7.2f} {'':>6}  {key}")
        continue
    ratio = best['wall'] / max(float(c_wall), 1e-3)
    print(f"{best['wall']:9.2f} {best['cpu']:7.2f} {float(c_wall):10.2f} {float(c_cpu):7.2f} {ratio:6.2f}  {key}", flush=True)
    total['k'] += best['wall']
    total['kc'] += best['cpu']
    total['c'] += float(c_wall)
    total['cc'] += float(c_cpu)
print(f"{total['k']:9.2f} {total['kc']:7.2f} {total['c']:10.2f} {total['cc']:7.2f} {total['k'] / max(total['c'], 1e-3):6.2f}  total")

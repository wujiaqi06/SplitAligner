#!/usr/bin/env python3
"""Bounded, shell-free native-port smoke; also useful on macOS/Linux.

Use --work-dir with a new directory. Host/platform are recorded: running this
on another OS is explicitly not a Windows acceptance result.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess

if not __debug__:
    raise RuntimeError('Optimized Python disables this test suite assertions; run without -O/PYTHONOPTIMIZE.')

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument('--work-dir', required=True, type=Path)
parser.add_argument('--require-windows', action='store_true')
args = parser.parse_args()
if args.require_windows and os.name != 'nt':
    parser.error('--require-windows requires a native Windows host; WSL/macOS/Linux are not Windows acceptance')
work = args.work_dir.resolve()
work.mkdir(parents=True, exist_ok=False)
env = dict(os.environ, LC_ALL='C', PERL_HASH_SEED='0', PERL_PERTURB_KEYS='0')
results = []

def command(name, argv, cwd, success=True):
    p = subprocess.run(argv, cwd=cwd, env=env, text=True, encoding='utf-8', capture_output=True)
    (work / (name + '.stdout')).write_text(p.stdout, encoding='utf-8')
    (work / (name + '.stderr')).write_text(p.stderr, encoding='utf-8')
    assert (p.returncode == 0) == success, name + ': ' + p.stderr
    results.append({'case': name, 'exit': p.returncode, 'expected_success': success})
    return p

species = '((A,B),(C,(D,E)));\n'
genes = 'full\t((A:1,B:2):3,(C:4,(D:5,E:6):7):8);\n' \
        'pruned\t(A:0,(C:4,(D:5,E:6):7):8);\n'
science = ['matrix_no_fuse.txt', 'matrix_with_fuse.txt', 'primitive_state.tsv', 'fusion_cells.tsv']
snapshots = {}
for line_end, view in ((e,v) for e in ('LF','CRLF') for v in ('whole','direct')):
    case = line_end+'_'+view
    run = work / (case + ' space')
    run.mkdir()
    (run/'species.nwk').write_bytes(species.replace('\n','\r\n' if line_end=='CRLF' else '\n').encode())
    (run/'genes.nwk').write_bytes(genes.replace('\n','\r\n' if line_end=='CRLF' else '\n').encode())
    argv=['perl', str(ROOT/'SplitAligner.pl'), '--mode', 'matrix', '--species', 'species.nwk', '--gene', 'genes.nwk', '--label', 'case', '--fusion-view',view]
    command(case+'_matrix',argv,run)
    snapshot={f:hashlib.sha256((run/('case.'+f)).read_bytes()).hexdigest() for f in science}
    command(case+'_audit',['perl',str(ROOT/'scripts/audit_compact_gene.pl'),'--manifest','case.run_manifest.json','--audit-gene','full','--gene-input','genes.nwk','--species-input','species.nwk'],run)
    command(case+'_reject_existing',argv,run,False)
    command(case+'_force',argv+['--force'],run)
    assert snapshot=={f:hashlib.sha256((run/('case.'+f)).read_bytes()).hexdigest() for f in science}
    snapshots[(line_end,view)]=snapshot
for view in ('whole','direct'):
    assert snapshots[('LF',view)]==snapshots[('CRLF',view)], 'LF/CRLF scientific results differ'
for field in ('matrix_no_fuse.txt','primitive_state.tsv'):
    assert snapshots[('LF','whole')][field]==snapshots[('LF','direct')][field], 'view changed primitive output'
unicode_run=work/'unicode contents'
unicode_run.mkdir()
(unicode_run/'species.nwk').write_text(species.replace('A','猫'),encoding='utf-8')
(unicode_run/'genes.nwk').write_text(genes.replace('A','猫').replace('full','基因'),encoding='utf-8')
(unicode_run/'selection.txt').write_text('基因\n',encoding='utf-8')
command('unicode_matrix',['perl',str(ROOT/'SplitAligner.pl'),'--mode','matrix','--species','species.nwk','--gene','genes.nwk','--label','case'],unicode_run)
command('unicode_audit',['perl',str(ROOT/'scripts/audit_compact_gene.pl'),'--manifest','case.run_manifest.json','--audit-gene-file','selection.txt','--gene-input','genes.nwk','--species-input','species.nwk'],unicode_run)
(work/'RESULTS.json').write_text(json.dumps({'status':'PASS','platform':platform.platform(),'native_windows':os.name=='nt','cases':results,'LF_CRLF_scientific_equality':True},indent=2)+'\n')
print(json.dumps({'status':'PASS','native_windows':os.name=='nt','cases':len(results)}))

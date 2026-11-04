"""An absent shared split has one group diagnosis; its member axes stay fused."""
import csv
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

cli = sys.argv[1]
with TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root / 'ref').write_text('((a,b),x,(c,d));\n')
    genes = {
        'shown': '((a:1,b:1):0,c:1,d:1);',
        'no_length': '((a,b),c,d);',
        'alternative': '((a,c),b,d);',
        'unresolved': '(a,b,c,d);',
    }
    (root / 'genes').write_text(''.join(f'{g}\t{t}\n' for g,t in genes.items()))
    for role,scope in [('free','full'), ('free','support'), ('fixed','full')]:
        out = root / (role+scope)
        subprocess.run([cli, '--reference', str(root/'ref'), '--genes', str(root/'genes'),
                        '--role',role,'--scope',scope,'--out',str(out),
                        '--unmatched-splits','off'], check=True)
        ctx = json.loads((out/'context.json').read_text())
        members = {a['id'] for a in ctx['axes'] if not a['terminal']}
        assert len(members) == 2
        rows = [json.loads(line) for line in (out/'composites.jsonl').read_text().splitlines()]
        assert [row['gene_id'] for row in rows] == list(genes)
        for row in rows:
            gene = row['gene_id']
            recovered = gene in ('shown','no_length')
            assert set(row['members']) == members and row['recovered'] == recovered
            if role == 'free':
                assert row['group_state'] == ('Mapped' if recovered else 'NA_fuse_topo')
            else:
                assert 'group_state' not in row
            assert row['numeric_available'] == (gene == 'shown' and scope == 'full')
            assert row['value'] == (0 if row['numeric_available'] else None)
        if role == 'free':
            count_tool = Path(__file__).resolve().parents[1] / 'examples/count_free_topology.py'
            for unit, weight in [('primitive', 2), ('coordinate', 1)]:
                data = subprocess.check_output([sys.executable, str(count_tool), str(out),
                                                '--unit', unit], text=True)
                reports = list(csv.DictReader(data.splitlines(), delimiter='\t'))
                assert [r['gene_id'] for r in reports] == list(genes)
                for r in reports:
                    present = r['gene_id'] in ('shown', 'no_length')
                    expected = [5 + weight, 4, 1, weight if present else 0, 0,
                                0 if present else weight]
                    assert [int(r[k]) for k in ('Total','Mapped','NA_struct','NA_fuse',
                                                'NA_topo','NA_fuse_topo')] == expected
        if scope == 'full':
            with (out/'raw_states.tsv').open() as f:
                states = list(csv.DictReader(f,delimiter='\t'))
            assert len(states) == len(genes)
            for row in states:
                assert all(row['A' + str(m)] == 'F' for m in members)
                assert list(row.values()).count('S') == 1

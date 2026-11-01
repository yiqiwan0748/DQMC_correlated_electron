"""Combine explicitly supplied execution receipts without counting oracle-only runs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tomllib


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('receipts', type=Path, nargs='+')
    args = parser.parse_args()
    manifest_bytes = args.manifest.read_bytes()
    manifest = tomllib.loads(manifest_bytes.decode())
    expected_hash = hashlib.sha256(manifest_bytes).hexdigest()
    cases = {c['id']: dict(c, arms={}, execution_history=[]) for c in manifest['cases']}
    receipt_info = []
    for path in args.receipts:
        if 'pilot' in path.name:
            raise ValueError('Timing pilots must remain separate from campaign results')
        raw = path.read_bytes()
        receipt = tomllib.loads(raw.decode())
        if receipt['manifest_sha256'] != expected_hash:
            raise ValueError(f'Manifest differs in {path}')
        receipt_info.append({'path': str(path), 'sha256': hashlib.sha256(raw).hexdigest(),
                             'job_id': receipt.get('slurm_job_id', 'UNKNOWN'),
                             'oracle_sha256': receipt.get('oracle_sha256')})
        for row in receipt['cases']:
            case = cases[row['id']]
            case['execution_history'].append({'receipt': path.name, 'status': row['status']})
            for arm in row.get('arms', []):
                method = arm['method']
                old = case['arms'].get(method)
                if old and old['status'] == 'PASS' and arm['status'] == 'UNRESOLVED' and arm.get('reason') in (
                        'Monte Carlo arm was explicitly disabled', 'Requires full dependency-backed runner'):
                    continue
                case['arms'][method] = dict(arm, receipt=path.name)
    rows = []
    for case in cases.values():
        arms = [case['arms'].get(m, {'method': m, 'status': 'NOT_RUN'}) for m in case['target_methods']]
        statuses = [a['status'] for a in arms]
        status = ('PASS' if all(s == 'PASS' for s in statuses) else
                  'FAIL' if 'FAIL' in statuses else 'UNRESOLVED')
        case['status'] = status
        case['arms'] = arms
        rows.append({'id': case['id'], 'family': case['family'], 'observable': case['observable'],
                     'status': status, 'arms': '; '.join(f"{a['method']}:{a['status']}" for a in arms)})
    counts = {s: sum(c['status'] == s for c in cases.values()) for s in ('PASS', 'FAIL', 'UNRESOLVED')}
    arm_counts = {s: sum(a['status'] == s for c in cases.values() for a in c['arms'])
                  for s in ('PASS', 'FAIL', 'UNRESOLVED', 'NOT_RUN')}
    result = {'scope': 'Initial 100 finite-system external comparisons; broader distinct-idea roadmap is separate',
              'manifest_sha256': expected_hash, 'case_counts': counts, 'arm_counts': arm_counts,
              'receipt_sources': receipt_info, 'cases': list(cases.values())}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    with (args.output_dir / 'case_status.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({'case_counts': counts, 'arm_counts': arm_counts}))


if __name__ == '__main__':
    main()

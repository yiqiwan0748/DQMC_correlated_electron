#!/usr/bin/env python3
"""Portable replay of selected, unchanged saved-data science; Python 3.11+ stdlib.

Run beside PACKET_MANIFEST.json: python3 replay_math.py --role both --output NEW_DIR
This adapter binds public members to a temporary historical relative layout.
It does not execute an ISAAC controller, producer, sampler or runtime-guard main.
"""
import argparse
import contextlib
import datetime
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import platform
import shutil
import tarfile
import tempfile
import time

SCRIPT_IDENTITIES = {
    'x14': ('replay_x14_saved_math.py', 'a7264e2fc12a3af3d66d9251ce5fc70b9854450fa333eb95297b144e4e1e69e1', 'x14_local_saved_math_replay.json', 'cae8ff79874a2a4a4543e9b3b6b3e0c9192f7864af4ef7cfb52b6465221da93b'),
    'x19': ('review_x19_saved_data.py', '203214e1fc29df88e8a4fc28ea5130c2b1bff3c84d3b2a0cd039b569d43788f1', 'x19_independent_review.json', '1633b1845bcd43452ae7be11ae82775df08dfb9a65d5d22a8620811a247595b2'),
}
PROJECTED_IDENTITIES = {
    'results/x14_saved_compare_recovery_v1_pc1_20261024/execution.json': '54edea5c09c0f5c412f27f7caecc682081ac18e10e3d18867c4740b0d6953fdc',
    'results/x14_full_independent_saved_audit_v1_pc1_20261024/audit.json': '285613cd21ac3c285eb8905cd2c83612da73726ed5f01606cabc915b6fd23300',
}
PRIVATE_INVENTORY_CHECK = 'preserved runtime inventory current_observed75_registry_PRIVATE.json'


def digest(b):
    return hashlib.sha256(b).hexdigest()


def write_json(p, value):
    p.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def safe_name(s):
    p = PurePosixPath(s)
    if p.is_absolute() or '..' in p.parts or '\\' in s or not p.parts:
        raise ValueError('Unsafe packet path: ' + s)
    return p


def materialize(packet, destination):
    manifest = json.loads((packet / 'PACKET_MANIFEST.json').read_text())
    for f in manifest['auxiliary_files']:
        b = (packet / Path(*safe_name(f['path']).parts)).read_bytes()
        assert len(b) == f['bytes'] and digest(b) == f['sha256'], f['path']
    identities = {}
    count = 0
    for archive in manifest['archives']:
        path = packet / Path(*safe_name(archive['path']).parts)
        b = path.read_bytes()
        assert len(b) == archive['bytes'] and digest(b) == archive['sha256'], archive['path']
        with tarfile.open(path, 'r:gz') as t:
            members = t.getmembers()
            assert len(members) == len(archive['members'])
            assert [m.name for m in members] == sorted(m.name for m in members)
            assert len({m.name for m in members}) == len(members)
            for member, record in zip(members, archive['members']):
                safe_name(member.name)
                assert member.isfile() and member.mode == 0o644 and member.mtime == 0
                assert member.name == record['public_member']
                data = t.extractfile(member).read()
                data.decode('utf-8')
                assert len(data) == record['public_bytes'] and digest(data) == record['public_sha256']
                source = safe_name(record['source_identity'])
                assert source.parts[0] in {'source', 'results', 'analysis_inputs', 'source_releases'}
                transform = record['transformation']
                if transform['type'] == 'BYTE_EXACT':
                    assert record['public_sha256'] == record['original_sha256']
                    assert record['public_bytes'] == record['original_bytes']
                else:
                    assert transform['type'] == 'JSON_PRIVATE_FIELD_PROJECTION'
                    assert record['source_identity'] in PROJECTED_IDENTITIES
                    assert record['original_sha256'] == PROJECTED_IDENTITIES[record['source_identity']]
                    assert transform['original_file_sha256'] == record['original_sha256']
                    assert transform['all_other_decoded_fields_exact'] is True
                    assert len(transform['omitted_fields']) == 1
                    assert transform['omitted_fields'][0]['json_pointer'] == '/whole_user_queue_before_PRIVATE'
                    assert 'whole_user_queue_before_PRIVATE' not in json.loads(data)
                target = destination / Path(*source.parts)
                if target.exists():
                    assert target.read_bytes() == data, record['source_identity']
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(data)
                identities[record['source_identity']] = record
                count += 1
    assert set(PROJECTED_IDENTITIES).issubset(identities)
    return manifest, identities, count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='Fresh output directory; must not exist')
    parser.add_argument('--role', choices=['x14', 'x19', 'both'], default='both')
    parser.add_argument('--packet', type=Path, default=Path(__file__).resolve().parent, help='Packet directory (default: beside this adapter)')
    parser.add_argument('--independent-dir', type=Path, help='Exact original review scripts/reports (default: PACKET/independent)')
    args = parser.parse_args()
    packet = args.packet.resolve()
    originals = (args.independent_dir or packet / 'independent').resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError('Output must be fresh: ' + str(output))
    roles = ['x14', 'x19'] if args.role == 'both' else [args.role]
    for role in roles:
        script, script_sha, report, report_sha = SCRIPT_IDENTITIES[role]
        assert digest((originals / script).read_bytes()) == script_sha, 'Original review script identity: ' + role
        assert digest((originals / report).read_bytes()) == report_sha, 'Original review report identity: ' + role
    output.mkdir(parents=True)
    started = time.monotonic()
    result = dict(schema='WEEKEND_PUBLIC_SAVED_MATH_ADAPTER_V1', status='RUNNING', role=args.role,
                  python=platform.python_version(), original_pinned_runtime_replay=False,
                  scientific_completions=0, physical_gates_changed=False, reviews=[])
    try:
        with tempfile.TemporaryDirectory(prefix='weekend-public-saved-math-') as temporary:
            base = Path(temporary)
            snap = base / 'reconciliation_20261025/weekend_state/snapshot'
            manifest, identities, count = materialize(packet, snap)
            result['verified_regular_UTF8_members'] = count
            result['packet_manifest_sha256'] = digest((packet / 'PACKET_MANIFEST.json').read_bytes())
            # These aliases only reconstruct the original reviewer scripts' relative
            # input layout. Every copied byte is already verified by the manifest.
            old_raw = base / 'results/x14_full_v3_cancelled_6467601/results'
            shutil.copytree(snap / 'results/x14_exact_sign_full_v3_20261023T2115Z', old_raw)
            x19_failed = snap / 'source/x19_full_saved_audit_v4_pc1_20261024/auditor/preserved_v3_failure'
            x19_failed.mkdir()
            shutil.copyfile(snap / 'results/x19_full_saved_audit_v3_20261023T2110Z/audit.json', x19_failed / 'audit.json')
            execution = json.loads((old_raw / 'execution.json').read_text())
            assert execution['status'] == 'UNRESOLVED' and execution['scientific_completions'] == 0
            for role in roles:
                script_name, script_sha, report_name, report_sha = SCRIPT_IDENTITIES[role]
                work = base / 'reconciliation_20261025/weekend_science_review'
                work.mkdir(parents=True, exist_ok=True)
                script_path = work / script_name
                script_bytes = (originals / script_name).read_bytes()
                script_path.write_bytes(script_bytes)
                assert digest(script_path.read_bytes()) == script_sha
                capture = io.StringIO()
                # Exact original code: __file__ is rebound to the temporary relative
                # layout; no AST helper, physical number, gate or algorithm is edited.
                with contextlib.redirect_stdout(capture):
                    exec(compile(script_bytes, str(script_path), 'exec'), {'__name__': 'portable_original_review_' + role, '__file__': str(script_path)})
                actual = json.loads((work / report_name).read_text())
                expected = json.loads((originals / report_name).read_text())
                omissions = []
                if role == 'x14':
                    key = 'reconciliation_20261025/weekend_state/snapshot/results/x14_full_independent_saved_audit_v1_pc1_20261024/audit.json'
                    identity = identities['results/x14_full_independent_saved_audit_v1_pc1_20261024/audit.json']
                    assert expected['input_sha256'][key] == identity['original_sha256']
                    assert actual['input_sha256'][key] == identity['public_sha256']
                    # Receipt identity changes only because its private queue field
                    # was projected out; all saved scientific results are compared.
                    expected['input_sha256'][key] = identity['public_sha256']
                else:
                    retained = [c for c in expected['checks'] if c['name'] != PRIVATE_INVENTORY_CHECK]
                    assert len(retained) == len(expected['checks']) - 1
                    assert not (snap / 'source/x19_full_saved_audit_v4_pc1_20261024/auditor/current_observed75_registry_PRIVATE.json').exists()
                    expected['checks'] = retained
                    expected['check_count'] = len(retained)
                    omissions.append(dict(name=PRIVATE_INVENTORY_CHECK, reason='Private captured process registry is intentionally absent. Its historical ISAAC qualification is preserved by hash/report; this portable replay does not requalify any runtime. No physical check is omitted.'))
                for field in ['python', 'elapsed_seconds', 'completed_utc']:
                    actual.pop(field, None)
                    expected.pop(field, None)
                assert actual == expected, 'Every retained decoded scientific review field must match original reviewed report: ' + role
                shutil.copyfile(work / report_name, output / report_name)
                (output / (role + '_stdout.txt')).write_text(capture.getvalue())
                result['reviews'].append(dict(role=role, status='PASS', original_script_sha256=script_sha,
                                              original_review_report_sha256=report_sha,
                                              all_retained_decoded_scientific_review_fields_equal=True,
                                              explicitly_omitted_private_runtime_metadata_checks=omissions,
                                              check_count=actual['check_count'], output_sha256=digest((output / report_name).read_bytes())))
        result['status'] = 'PASS'
    except Exception as error:
        result['status'] = 'FAIL'
        result['error'] = type(error).__name__ + ': ' + str(error)
        raise
    finally:
        result['elapsed_seconds'] = time.monotonic() - started
        result['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        write_json(output / 'adapter_receipt.json', result)
    print(json.dumps({'status': result['status'], 'role': args.role, 'reviews': len(result['reviews']), 'elapsed_seconds': result['elapsed_seconds']}, sort_keys=True))


if __name__ == '__main__':
    main()

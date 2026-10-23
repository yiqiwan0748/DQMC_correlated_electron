#!/usr/bin/env python3
"""Independent stable-residual audit of the recovered, completed precision retest."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import tomllib


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inner(a, b):
    return sum(x.conjugate() * y for x, y in zip(a, b))


def mv(a, v):
    return [sum(x * y for x, y in zip(row, v)) for row in a]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('raw', 'tranche', 'original-tranche', 'base-ledger', 'output'):
        p.add_argument('--' + name, required=True, type=Path)
    a = p.parse_args()
    assert not a.output.exists(), 'Use a new output'
    checks = []

    def gate(name, ok):
        checks.append({'name': name, 'passed': bool(ok)})
        if not ok:
            raise ValueError(name)

    cfg = tomllib.loads((a.tranche / 'manifest.toml').read_text())
    original = tomllib.loads((a.original_tranche / 'manifest.toml').read_text())
    gate('Scientific gates identical to original eight-sweep attempt', cfg['gates'] == original['gates'])
    g = cfg['gates']
    transfer = json.loads((a.raw / 'pc_transfer.json').read_text())
    for row in transfer['files']:
        path = a.raw / row['local']
        gate('Recovered original bytes: ' + row['local'], path.stat().st_size == row['bytes'] and sha(path) == row['sha256'])
    native = a.raw / 'native_source'
    for rel, key in [('Project.toml', 'project_sha256'), ('Manifest.toml', 'environment_manifest_sha256'),
                     ('analysis/compute_s2_ed.jl', 'ed_source_sha256'),
                     ('analysis/audit_u4_l4_thermal_ed.jl', 'thermal_source_sha256'),
                     ('analysis/compute_s2_itensor_dmrg.jl', 'dmrg_source_sha256')]:
        gate('Native source closure: ' + rel, sha(native / rel) == cfg[key])
    target = tomllib.loads((a.raw / 'results/target.toml').read_text())
    oracle = json.loads((a.raw / 'results/oracle.json').read_text())
    for d in (target, oracle):
        gate('Executed manifest pin', d['manifest_sha256'] == sha(a.tranche / 'manifest.toml'))
    gate('Native adapter pin', target['adapter_sha256'] == sha(a.tranche / 'target.jl'))
    gate('Independent oracle pin', oracle['oracle_source_sha256'] == sha(a.tranche / 'oracle_quspin.py'))
    t, o = target['X09'], oracle['X09']
    expected = {(u, d) for u in range(16) for d in range(16) if u.bit_count() == d.bit_count() == 2}
    oq, tq = [tuple(x) for x in o['coordinates']], [tuple(x) for x in t['coordinates']]
    gate('Complete independent and user sector coordinates', len(oq) == len(tq) == 36 and set(oq) == set(tq) == expected)
    perm = [oq.index(q) for q in tq]
    h = [[o['Hamiltonian'][i][j] for j in perm] for i in perm]
    proj = [[o['ground_projector'][i][j] for j in perm] for i in perm]
    for key in ('Hamiltonian', 'MPO_qn_Hamiltonian', 'MPO_dense_Hamiltonian'):
        gate('Actual matrix agrees with QuSpin: ' + key, max(abs(t[key][i][j] - h[i][j]) for i in range(36) for j in range(36)) <= g['matrix_abs'])
    gate('Independent unique ground and gap', o['ground_dimension'] == 1 and o['gap'] > g['ground_gap_min'])
    gate('Hermitian oracle', max(abs(h[i][j] - h[j][i]) for i in range(36) for j in range(36)) <= g['matrix_abs'])
    expected_params = {k: cfg['X09'][k] for k in ('maxdim_schedule', 'cutoff', 'eigsolve_tol', 'eigsolve_krylovdim', 'eigsolve_maxiter', 'noise', 'initial_link_dimension')}
    expected_params['nsweeps'] = len(cfg['X09']['maxdim_schedule'])
    gate('Declared optimizer parameters and sixteen sweeps', t['optimizer_parameters'] == expected_params and expected_params['nsweeps'] == 16)
    labels = cfg['X09']['starts'] + cfg['X09']['challenge_states']
    gate('Four starts and both manufactured challenges retained', len(t['rows']) == 6 and {r['label'] for r in t['rows']} == set(labels))
    metrics = []
    for row in t['rows']:
        v = [complex(x, y) for x, y in zip(row['coeff_real'], row['coeff_imag'])]
        gate('Complete complex state: ' + row['label'], len(v) == 36)
        norm = inner(v, v).real
        hv = mv(h, v)
        energy = inner(v, hv).real / norm
        stable_variance = sum(abs(x - energy * y) ** 2 for x, y in zip(hv, v)) / norm
        fidelity = inner(v, mv(proj, v)).real / norm
        excess = energy - o['ground_energy']
        is_native = row['label'] in cfg['X09']['starts']
        certified = (stable_variance <= g['variance_max'] and -g['energy_abs'] <= excess <= g['energy_abs'] and fidelity >= g['ground_fidelity_min'])
        gate('Normalized physical sector: ' + row['label'], abs(norm - 1) <= g['norm_abs'])
        gate('Stable residual agrees with reported variance: ' + row['label'], abs(stable_variance - row['variance']) <= g['energy_abs'])
        gate('Native passes or manufactured excited challenge is rejected: ' + row['label'], certified == is_native)
        metrics.append({'label': row['label'], 'kind': row['kind'], 'stable_variance': stable_variance,
                        'energy_excess': excess, 'ground_projector_fidelity': fidelity, 'certified': certified})
    replay = json.loads((a.raw / 'fresh_comparison_replay.json').read_text())
    reported = json.loads((a.raw / 'results/comparison.json').read_text())
    gate('Fresh frozen-comparator replay exactly reproduces original JSON', replay == reported)
    gate('Both frozen comparison rows pass', {r['id']: r['status'] for r in replay['rows']} == {'X09': 'PASS', 'X12': 'PASS'})
    base = json.loads(a.base_ledger.read_text())
    gate('Base X09 unresolved and X12 already passed', {r['id']: r['status'] for r in base['rows'] if r['id'] in ('X09', 'X12')} == {'X09': 'UNRESOLVED', 'X12': 'PASS'})
    bound = [a.tranche / x for x in ('manifest.toml', 'target.jl', 'oracle_quspin.py', 'compare.py')]
    bound += [a.original_tranche / 'manifest.toml', a.base_ledger, a.raw / 'pc_transfer.json', a.raw / 'fresh_remote_transfer.json', a.raw / 'fresh_comparison_replay.json']
    bound += [a.raw / x['local'] for x in transfer['files']]
    result = {'schema': 'recovered_ground_precision_root_audit_v1', 'status': 'PASS', 'job_id': '6460656',
              'checks': checks, 'check_count': len(checks), 'scientific_status': {'X09': 'PASS', 'X12': 'PASS'},
              'new_distinct_passes': 1, 'X12_not_double_counted': True, 'metrics': metrics,
              'source_and_raw_inputs': [{'path': str(x), 'bytes': x.stat().st_size, 'sha256': sha(x)} for x in bound],
              'auditor_sha256': sha(Path(__file__)), 'scientific_solver_executed': False,
              'scope': 'The frozen four-site fixed-spin certification fixture only; original eight-sweep failure and all gates preserved.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': result['status'], 'checks': len(checks), 'new_distinct_passes': 1}))


if __name__ == '__main__':
    main()

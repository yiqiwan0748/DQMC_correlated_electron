#!/usr/bin/env python3
"""Export independent QuSpin Hubbard-chain references for external benchmarks.

This program imports no user ED, DQMC, DMRG or Renyi implementation. All fermion
Hamiltonians and partial traces are constructed by QuSpin. SciPy diagonalizes the
resulting small matrices. The finite-dtau target is explicitly V-then-K.
"""
from __future__ import annotations
import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.special import logsumexp
import quspin
from quspin.basis import spinful_fermion_basis_1d
from quspin.operators import hamiltonian

SOURCES = [
    'https://quspin.github.io/QuSpin/examples/example6.html',
    'https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html',
    'https://arxiv.org/abs/1804.06782',
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked_real(value, label, tolerance=2e-10):
    value = complex(value)
    if not np.isfinite(value.real) or not np.isfinite(value.imag) or abs(value.imag) > tolerance:
        raise ValueError(f'{label}: nonfinite or unexpectedly complex value {value}')
    return float(value.real)


def matrix(static, basis, *, hermitian=False):
    # QuSpin prints diagnostics for every construction; retain only final receipts.
    with contextlib.redirect_stdout(io.StringIO()):
        operator = hamiltonian(static, [], basis=basis, dtype=np.float64,
                               check_herm=hermitian, check_symm=False, check_pcon=False)
    return np.asarray(operator.toarray(), dtype=np.float64)


def chain_bonds(L, boundary):
    if boundary not in ('obc', 'pbc', 'apbc'):
        raise ValueError(f'Unsupported boundary: {boundary}')
    # L=2 PBC would be a two-directed-bond convention; exclude that ambiguity.
    if L < 3 and boundary != 'obc':
        raise ValueError('Use OBC for L<3; doubled PBC bonds are deliberately unsupported')
    bonds = [(i, i + 1, 1.0) for i in range(L - 1)]
    if boundary != 'obc':
        bonds.append((L - 1, 0, -1.0 if boundary == 'apbc' else 1.0))
    return bonds


def construct_model(model, ensemble):
    L = int(model['L'])
    if L < 1 or L > 6:
        raise ValueError('This dense benchmark oracle is restricted to 1 <= L <= 6')
    if ensemble != 'canonical_ground' and L > 4:
        raise ValueError('Grand-canonical dense thermal targets are restricted to L<=4')
    kwargs = {}
    if ensemble == 'canonical_ground':
        kwargs['Nf'] = (int(model['Nup']), int(model['Ndown']))
    basis = spinful_fermion_basis_1d(L, **kwargs)
    t, U, mu = (float(model.get(k, d)) for k, d in
                 [('t', 1.0), ('U', 0.0), ('mu_centered', 0.0)])
    if not all(np.isfinite(v) for v in (t, U, mu)):
        raise ValueError('Nonfinite model parameter')
    directed = []
    for i, j, phase in chain_bonds(L, model['boundary']):
        directed.extend([[-t * phase, i, j], [-t * phase, j, i]])
    Khop = matrix([['+-|', directed], ['|+-', directed]], basis, hermitian=True)
    number_up = [np.diag(matrix([['n|', [[1.0, i]]]], basis)) for i in range(L)]
    number_down = [np.diag(matrix([['|n', [[1.0, i]]]], basis)) for i in range(L)]
    n_up = np.array(number_up)
    n_down = np.array(number_down)
    n_total = (n_up + n_down).sum(axis=0)
    d_local = n_up * n_down
    d_total = d_local.sum(axis=0)
    # Build the interaction independently using QuSpin's n|n string, then check
    # its diagonal against the separately constructed one-species occupations.
    D = matrix([['n|n', [[1.0, i, i] for i in range(L)]]], basis, hermitian=True)
    if np.max(np.abs(D - np.diag(d_total))) > 1e-12:
        raise ValueError('QuSpin local occupation/interaction construction disagrees')
    Vcenter = U * (d_total - n_total / 2 + L / 4)
    K = Khop - np.diag(mu * n_total)
    Hgenerator = K + np.diag(Vcenter)
    if np.max(np.abs(Hgenerator - Hgenerator.T)) > 1e-12:
        raise ValueError('Reference Hamiltonian is not Hermitian')
    return {'basis': basis, 'L': L, 'U': U, 'mu': mu, 'Khop': Khop, 'K': K,
            'Vcenter': Vcenter, 'Hgenerator': Hgenerator,
            'n_up': n_up, 'n_down': n_down, 'n_total': n_total,
            'd_local': d_local, 'd_total': d_total}


def prepare_state(model, ensemble, beta=None, dtau=None):
    data = construct_model(model, ensemble)
    H = data['Hgenerator']
    if ensemble == 'canonical_ground':
        energies, vectors = eigh(H, check_finite=True)
        tolerance = 1e-10 * max(1.0, np.max(np.abs(energies)))
        degeneracy = int(np.sum(np.abs(energies - energies[0]) <= tolerance))
        if degeneracy != 1:
            raise ValueError(f'Ground state is {degeneracy}-fold degenerate; pure-state observables have no unique oracle')
        state = vectors[:, 0]
        residual = np.linalg.norm(H @ state - energies[0] * state)
        if residual > 1e-10:
            raise ValueError(f'Ground eigenpair residual failed: {residual}')
        metadata = {'ground_energy_centered_generator': float(energies[0]),
                    'ground_degeneracy': degeneracy,
                    'ground_residual': float(residual), 'state_kind': 'pure_state'}
    elif ensemble in ('grand_canonical_thermal', 'grand_canonical_trotter'):
        beta = float(beta)
        if not np.isfinite(beta) or beta <= 0:
            raise ValueError('beta must be finite and positive')
        if ensemble == 'grand_canonical_thermal':
            energies, vectors = eigh(H)
            logweights = -beta * energies
            probabilities = np.exp(logweights - logsumexp(logweights))
            state = (vectors * probabilities) @ vectors.T
            metadata = {'log_partition_function': float(logsumexp(logweights)),
                        'state_kind': 'grand_canonical_density_matrix',
                        'thermal_operator': 'exp(-beta Hc) / Tr(exp(-beta Hc))'}
        else:
            dtau = float(dtau)
            slices = int(round(beta / dtau)) if dtau > 0 else 0
            if not np.isfinite(dtau) or dtau <= 0 or slices < 1 or abs(beta - slices * dtau) > 1e-10:
                raise ValueError('beta/dtau must be a positive integer')
            kvals, kvecs = eigh(data['K'])
            Khalf = (kvecs * np.exp(-dtau * kvals / 2)) @ kvecs.T
            Khalf_inv = (kvecs * np.exp(dtau * kvals / 2)) @ kvecs.T
            Vstep = np.diag(np.exp(-dtau * data['Vcenter']))
            symmetric_step = Khalf @ Vstep @ Khalf
            lambdas, vectors = eigh(symmetric_step)
            if np.min(lambdas) <= 0:
                raise ValueError('Finite-step propagator lost positive eigenvalues')
            logweights = slices * np.log(lambdas)
            probabilities = np.exp(logweights - logsumexp(logweights))
            rho_symmetric = (vectors * probabilities) @ vectors.T
            # B=Vstep*Khalf^2=Khalf^{-1}*symmetric_step*Khalf.
            # Restore the actual asymmetric V K density at the measurement cut.
            state = Khalf_inv @ rho_symmetric @ Khalf
            metadata = {'log_partition_function': float(logsumexp(logweights)),
                        'state_kind': 'normalized_asymmetric_transfer_matrix',
                        'thermal_operator': '[exp(-dtau Vcenter) exp(-dtau (Khop-mu_c N))]^M / trace',
                        'trotter_order': 'V_K', 'time_slices': slices,
                        'density_Hermiticity_error': float(np.max(np.abs(state - state.T)))}
        if abs(np.trace(state) - 1) > 1e-10:
            raise ValueError('Density normalization failed')
        # A finite-dtau asymmetric density is similar to a positive matrix; it
        # need not be Hermitian in the ordinary occupation-basis metric.
    else:
        raise ValueError(f'Unsupported ensemble: {ensemble}')
    data.update(state=state, metadata=metadata)
    return data


def expect(data, operator):
    state = data['state']
    if state.ndim == 1:
        if np.ndim(operator) == 1:
            return checked_real(np.vdot(state, operator * state), 'Expectation')
        return checked_real(np.vdot(state, operator @ state), 'Expectation')
    if np.ndim(operator) == 1:
        return checked_real(np.dot(np.diag(state), operator), 'Expectation')
    return checked_real(np.trace(state @ operator), 'Expectation')




def quspin_reduced_density(data, sites):
    """Handle QuSpin's documented pure/mixed batch-shape distinction explicitly."""
    state = data['state']
    reduced = np.asarray(data['basis'].partial_trace(
        state, sub_sys_A=(sites, sites), return_rdm='A',
        enforce_pure=state.ndim == 1))
    original_shape = list(reduced.shape)
    dim_a = 4**len(sites)
    # Mixed-state partial_trace may retain a leading singleton sample axis.
    # np.trace on (1,dA,dA) traces the wrong axes and returns an array, so remove
    # exactly that documented singleton. Never squeeze a physical matrix axis.
    if reduced.shape == (1, dim_a, dim_a):
        reduced = reduced[0]
    if reduced.shape != (dim_a, dim_a):
        raise ValueError(f'Unexpected QuSpin partial-trace shape {original_shape}')
    data['metadata'].setdefault('partial_trace_original_shapes', {})[
        ','.join(map(str, sites))] = original_shape
    return reduced


def measure(data, case):
    observable = case['observable']
    L, U = data['L'], data['U']
    if observable == 'energy_uncentered':
        return expect(data, data['Khop'] + np.diag(U * data['d_total']))
    if observable == 'energy_centered':
        return expect(data, data['Khop'] + np.diag(data['Vcenter']))
    if observable == 'grand_energy_centered':
        return expect(data, data['Hgenerator'])
    if observable == 'kinetic_total':
        return expect(data, data['Khop'])
    if observable in ('double_occupancy_total', 'double_occupancy_per_site'):
        return expect(data, data['d_total']) / (L if observable.endswith('per_site') else 1)
    if observable in ('density_total', 'density_per_site'):
        return expect(data, data['n_total']) / (L if observable.endswith('per_site') else 1)
    if observable in ('local_n', 'local_double'):
        i = int(case['sites'][0])
        if not 0 <= i < L:
            raise ValueError('Invalid local site')
        return expect(data, data['n_up'][i] + data['n_down'][i] if observable == 'local_n' else data['d_local'][i])
    if observable in ('spin_z_correlation', 'onebody_up'):
        i, j = map(int, case['sites'])
        if not all(0 <= site < L for site in (i, j)):
            raise ValueError('Invalid correlation sites')
        if observable == 'spin_z_correlation':
            sz = (data['n_up'] - data['n_down']) / 2
            return expect(data, sz[i] * sz[j])
        if i == j:
            return expect(data, data['n_up'][i])
        Cij = matrix([['+-|', [[1.0, i, j]]]], data['basis'])
        return expect(data, Cij)
    if observable == 'S2':
        sites = list(map(int, case['subsystem_sites']))
        if len(set(sites)) != len(sites) or not all(0 <= i < L for i in sites):
            raise ValueError('Invalid subsystem sites')
        if not sites:
            return 0.0
        if len(sites) == L:
            purity = 1.0 if data['state'].ndim == 1 else np.trace(data['state'] @ data['state'])
        else:
            reduced = quspin_reduced_density(data, sites)
            if abs(np.trace(reduced) - 1) > 1e-9:
                raise ValueError('QuSpin reduced density matrix is not normalized')
            purity = np.trace(reduced @ reduced)
        purity = checked_real(purity, 'Subsystem purity')
        if not 0 < purity <= 1 + 1e-9:
            raise ValueError(f'Purity outside physical interval: {purity}')
        return float(-np.log(purity))
    raise ValueError(f'Unsupported observable: {observable}')


def source_receipt(manifest_path):
    distributions = {}
    distribution_records = {}
    for name in ('quspin', 'quspin-extensions', 'parallel-sparse-tools', 'numpy', 'scipy'):
        try:
            distribution = importlib.metadata.distribution(name)
            distributions[name] = distribution.version
            installed_record = distribution.read_text('RECORD')
            distribution_records[name] = (hashlib.sha256(installed_record.encode()).hexdigest()
                                          if installed_record is not None else None)
        except importlib.metadata.PackageNotFoundError:
            distributions[name] = 'not-installed-as-this-distribution-name'
    paths = {'exporter': str(Path(__file__).resolve()), 'quspin_entry': str(Path(quspin.__file__).resolve())}
    return {'oracle_implementation': 'QuSpin Hamiltonian and fermionic partial trace; NumPy/SciPy dense spectral arithmetic',
            'primary_sources': SOURCES, 'manifest_sha256': sha(manifest_path),
            'python_version': sys.version, 'platform': platform.platform(),
            'distributions': distributions, 'installed_distribution_RECORD_sha256': distribution_records,
            'source_paths': paths,
            'source_sha256': {key: sha(path) for key, path in paths.items()},
            'threads': {key: os.environ.get(key) for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
            'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
            'user_implementation_imported': False}


def selfcheck():
    # Analytic open Hubbard dimer checks guard the external package construction,
    # energy convention and spin normalization before exporting any oracle row.
    data = prepare_state({'L': 2, 't': 1, 'U': 4, 'mu_centered': 0,
                          'boundary': 'obc', 'Nup': 1, 'Ndown': 1}, 'canonical_ground')
    expected_energy = (4 - np.sqrt(4**2 + 16)) / 2
    expected_double = (1 - 4 / np.sqrt(4**2 + 16)) / 2
    observed_energy = measure(data, {'observable': 'energy_uncentered'})
    observed_double = measure(data, {'observable': 'double_occupancy_total'})
    expected_onebody = 2 / np.sqrt(4**2 + 16)
    expected_czz = -(1 - expected_double) / 4
    expected_s2 = -np.log((expected_double**2 + (1-expected_double)**2) / 2)
    comparisons = {
        'energy': (observed_energy, expected_energy),
        'double': (observed_double, expected_double),
        'onebody': (measure(data, {'observable': 'onebody_up', 'sites': [0, 1]}), expected_onebody),
        'Czz': (measure(data, {'observable': 'spin_z_correlation', 'sites': [0, 1]}), expected_czz),
        'S2': (measure(data, {'observable': 'S2', 'subsystem_sites': [0]}), expected_s2),
    }
    if any(abs(actual-expected) > 1e-11 for actual, expected in comparisons.values()):
        raise ValueError(f'Independent analytic Hubbard-dimer plumbing checks failed: {comparisons}')
    return {f'analytic_dimer_{key}_abs_error': abs(actual-expected)
            for key, (actual, expected) in comparisons.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    cases = manifest['cases']
    ids = [case['id'] for case in cases]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate case IDs')
    if args.output_dir.exists():
        raise FileExistsError('Use a new output directory; oracle output is never overwritten')
    args.output_dir.mkdir(parents=True)
    started = time.monotonic()
    receipt = source_receipt(args.manifest)
    try:
        receipt['selfcheck'] = selfcheck()
    except Exception as exc:
        receipt.update(status='FAILED_ORACLE_PLUMBING_GUARD', failure_reason=f'{type(exc).__name__}: {exc}')
        (args.output_dir / 'source_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        raise
    cache, results, values = {}, [], {}
    for case in cases:
        key = json.dumps({key: case.get(key) for key in ('model', 'ensemble', 'beta', 'dtau')}, sort_keys=True)
        begin = time.monotonic()
        try:
            if key not in cache:
                cache[key] = prepare_state(case['model'], case['ensemble'], case.get('beta'), case.get('dtau'))
            data = cache[key]
            value = measure(data, case)
            result = {'id': case['id'], 'status': 'ORACLE_READY', 'reference_value': value,
                      'case': case, 'basis_dimension': int(data['basis'].Ns),
                      'state_metadata': data['metadata'], 'elapsed_seconds': time.monotonic() - begin}
            values[case['id']] = value
        except Exception as exc:
            result = {'id': case['id'], 'status': 'UNRESOLVED', 'case': case,
                      'failure_reason': f'{type(exc).__name__}: {exc}',
                      'traceback': __import__('traceback').format_exc(),
                      'elapsed_seconds': time.monotonic() - begin}
        results.append(result)
        print(json.dumps({'id': result['id'], 'status': result['status'],
                          'reference_value': result.get('reference_value'),
                          'failure_reason': result.get('failure_reason')}), flush=True)
    output = args.output_dir / 'oracle_results.jsonl'
    output.write_text(''.join(json.dumps(result, sort_keys=True) + '\n' for result in results))
    toml_path = args.output_dir / 'oracle_values.toml'
    toml_path.write_text('[values]\n' + ''.join(f'"{case_id}" = {value:.17g}\n' for case_id, value in values.items()))
    receipt.update({'elapsed_seconds': time.monotonic() - started, 'case_count': len(cases),
                    'ready_count': len(values), 'unresolved_count': len(cases) - len(values),
                    'status': 'PASS_ORACLE_EXPORT' if len(values) == len(cases) else 'UNRESOLVED_ORACLE_EXPORT',
                    'output_sha256': {path.name: sha(path) for path in (output, toml_path)},
                    'note': 'ORACLE_READY means an independent reference is available; it is not a passing result for the user package.'})
    (args.output_dir / 'source_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))
    return 0 if len(values) == len(cases) else 2


if __name__ == '__main__':
    sys.exit(main())

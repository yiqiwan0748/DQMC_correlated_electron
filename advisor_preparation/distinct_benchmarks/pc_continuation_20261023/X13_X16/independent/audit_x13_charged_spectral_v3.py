#!/usr/bin/env python3
"""Independent X13 charged spectral receipt audit (NumPy and stdlib only).

This program never starts a native scientific solver, samples, submits, or writes
inside its inputs. Optional frozen comparator replay is supplemental, and is
unavailable until the complete approved source/archive/review closure passes.
An integrity PASS is distinct from the finite-cluster scientific decision.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tarfile
import tempfile
import tomllib

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import numpy as np

SCHEMA = 'X13_CHARGED_SPECTRAL_INDEPENDENT_AUDIT_V1'
RAW_SCHEMA = 'X13_CHARGED_SPECTRAL_RAW_V1'
RELEASE = 'x13_charged_spectral_v1'
CATALOGUE_SHA = '59a780c5994db50a56da260fc51ae639f22b674b3bc24f878cde1a6e7bc31818'
NATIVE_PINS = {
    'Project.toml': 'ea3b7eacfcf9a36e4652654a5e7b573ef27529edfbad2e272da2bafca1a7f58f',
    'Manifest.toml': '5cefaf1acdb4fa23a456c2c58f20fb1fcebd54b3cd27bc448cb65b0fcfecb90e',
    'analysis/compute_s2_ed.jl': 'f8f94cfc13e11bb8b5e64e8ac8f4ade5ed3a13b3976c83b5ab5db89d86ba5f28',
}
GATES = dict(matrix_abs=1e-10, eigen_residual_abs=1e-10,
             orthogonality_abs=1e-10, cluster_diameter_max=1e-8,
             cluster_separation_min=1e-6, cross_energy_abs=1e-9,
             cross_projector_abs=1e-9, cross_source_density_abs=1e-9,
             cross_weight_abs=1e-9, cross_threshold_abs=1e-9,
             car_abs=1e-10, occupation_abs=1e-10, zeroth_moment_abs=1e-10,
             first_moment_abs=1e-10, invisible_max=1e-12,
             visible_min=1e-9, negative_weight_min=-1e-12)
CASE_GRID = [('U2_pbc', 2., 'pbc', 1.), ('U2_apbc', 2., 'apbc', 1.),
             ('U4_pbc', 4., 'pbc', 2.), ('U4_apbc', 4., 'apbc', 2.)]
SOURCE_ROLES = ('manifest.toml', 'PROTOCOL.md', 'raw_io.py', 'target.jl',
                'reference.py', 'compare.py', 'run_candidate.py',
                'candidate.sbatch', 'README.md', 'test_candidate.py', 'PROVENANCE.json')
CUTOFF = 1792755300
EPS64 = 64 * np.finfo(float).eps


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite_tree(value):
    if isinstance(value, float):
        need(math.isfinite(value), 'Nonfinite scalar')
    elif isinstance(value, dict):
        for v in value.values():
            finite_tree(v)
    elif isinstance(value, list):
        for v in value:
            finite_tree(v)


def json_object(pairs):
    out = {}
    for k, v in pairs:
        need(k not in out, 'Duplicate JSON key: ' + k)
        out[k] = v
    return out


def load(path):
    path = Path(path)
    text = path.read_text(encoding='utf-8-sig')
    obj = tomllib.loads(text) if path.suffix == '.toml' else json.loads(text, object_pairs_hook=json_object)
    finite_tree(obj)
    return obj


def safe_relative(name):
    need(isinstance(name, str) and name and '\\' not in name and ':' not in name, 'Relative POSIX source/receipt role')
    p = PurePosixPath(name)
    need(not p.is_absolute() and all(x not in ('', '.', '..') for x in name.split('/')), 'No path traversal/aliases')
    devices = {'CON', 'PRN', 'AUX', 'NUL', *(f'COM{i}' for i in range(1, 10)), *(f'LPT{i}' for i in range(1, 10))}
    need(all(not x.endswith((' ', '.')) and x.split('.')[0].upper() not in devices for x in p.parts), 'No Windows path aliases')
    return p


class Audit:
    def __init__(self):
        self.checks = []
        self.inputs = {}
        self.numerical = []
        self.comparisons = []

    def check(self, name, ok, detail=None):
        self.checks.append(dict(name=name, passed=bool(ok), detail=detail))
        return bool(ok)

    def qualified(self, name, value, limit):
        ok = math.isfinite(float(value)) and value <= limit
        self.numerical.append(dict(name=name, passed=bool(ok), value=float(value), limit=limit))
        return bool(ok)

    def compare(self, name, value, limit):
        ok = math.isfinite(float(value)) and value <= limit
        self.comparisons.append(dict(name=name, passed=bool(ok), value=float(value), limit=limit))
        return bool(ok)

    def file(self, path):
        p = Path(path).resolve()
        need(p.is_file() and not Path(path).is_symlink(), 'Required regular receipt: ' + str(path))
        h = sha(p)
        if str(p) in self.inputs:
            need(self.inputs[str(p)] == h, 'Input changed while auditing')
        self.inputs[str(p)] = h
        return p

    def read(self, path):
        return load(self.file(path))

    def integrity(self):
        return all(c['passed'] for c in self.checks)


def maxabs(a):
    a = np.asarray(a)
    return float(np.max(np.abs(a))) if a.size else 0.


def close(a, b, atol):
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)) and maxabs(a - b) <= atol


def numeric_list(value):
    need(isinstance(value, list) and all(type(x) in (int, float) and math.isfinite(x) for x in value), 'Finite numeric array, no string/boolean coercion')
    return value


def coordinate_pairs(value, L=4):
    need(isinstance(value, list) and all(isinstance(c, list) and len(c) == 2 and all(type(x) is int and 0 <= x < (1 << L) for x in c) for c in value), 'Exact integer occupation-coordinate pairs')
    return [tuple(c) for c in value]


def matrix(obj, rows=None, columns=None):
    r, c = obj['rows'], obj['columns']
    need(type(r) is int and type(c) is int and 0 <= r <= 256 and 0 <= c <= 256, 'Dense matrix dimensions')
    need((rows is None or rows == r) and (columns is None or columns == c), 'Exact dense matrix shape')
    real, imag = obj['real'], obj['imag']
    need(isinstance(real, list) and isinstance(imag, list) and len(real) == len(imag) == r * c, 'Complete separate row-major real/imag data')
    numeric_list(real)
    numeric_list(imag)
    z = np.asarray(real, dtype=float) + 1j * np.asarray(imag, dtype=float)
    need(np.all(np.isfinite(z)), 'Finite complex matrix')
    return z.reshape(r, c)


def sparse(obj, dimension=256):
    need(obj['rows'] == obj['columns'] == dimension, 'Full-Fock sparse shape')
    fields = [obj[k] for k in ('row0', 'column0', 'real', 'imag')]
    need(all(isinstance(v, list) for v in fields) and len({len(v) for v in fields}) == 1, 'Complete sparse triplets')
    z = np.zeros((dimension, dimension), complex)
    seen = set()
    for i, j, x, y in zip(*fields):
        need(type(i) is int and type(j) is int and 0 <= i < dimension and 0 <= j < dimension and (i, j) not in seen, 'Unique in-range sparse entries')
        need(type(x) in (int, float) and type(y) in (int, float) and math.isfinite(x) and math.isfinite(y) and (x != 0 or y != 0), 'Finite numeric nonzero sparse entries')
        seen.add((i, j))
        z[i, j] = complex(x, y)
    return z


def coordinate_grid(L=4):
    return [(u, d) for u in range(1 << L) for d in range(1 << L)]


def annihilator(L, mode):
    """Independent occupation algebra; down modes follow all up modes."""
    n = 1 << L
    C = np.zeros((n*n, n*n), complex)
    spin, j = divmod(mode, L)
    for u in range(n):
        for d in range(n):
            mask = u if spin == 0 else d
            if mask & (1 << j):
                lower = (mask & ((1 << j) - 1)).bit_count()
                parity = lower + (u.bit_count() if spin else 0)
                new = mask ^ (1 << j)
                i = new*n + d if spin == 0 else u*n + new
                C[i, u*n+d] = -1. if parity % 2 else 1.
    return C


def canonical_hamiltonian(L, t, U, mu, boundary):
    """Bitwise Hubbard K=kinetic+U D-mu N, without candidate/native imports."""
    need(boundary in ('pbc', 'apbc'), 'Boundary condition')
    n = 1 << L
    H = np.zeros((n*n, n*n), complex)
    for u in range(n):
        for d in range(n):
            col = u*n+d
            H[col, col] = U*(u & d).bit_count() - mu*(u.bit_count()+d.bit_count())
            for spin, mask in enumerate((u, d)):
                for j in range(L):
                    k = (j+1) % L
                    seam = -1 if boundary == 'apbc' and j == L-1 else 1
                    for to, source in ((j, k), (k, j)):
                        if mask >> source & 1 and not mask >> to & 1:
                            lo, hi = sorted((to, source))
                            between = ((1 << hi)-1) ^ ((1 << (lo+1))-1)
                            sign = -1 if (mask & between).bit_count() % 2 else 1
                            moved = mask ^ (1 << to) ^ (1 << source)
                            row = moved*n+d if spin == 0 else u*n+moved
                            H[row, col] += -t*seam*sign
    return H


def car_error(operators):
    n = operators[0].shape[0]
    I = np.eye(n)
    err = 0.
    for a, A in enumerate(operators):
        for b, B in enumerate(operators):
            err = max(err, maxabs(A@B+B@A), maxabs(A@B.conj().T+B.conj().T@A-(I if a == b else 0)))
    return err


def creation_bridge(raw_C, raw_coordinates, L=4):
    """Algebraic sign derivation from the vacuum, with no fitted gauge.

    Canonical ket |n0...n7> = (c0†)^n0 ... (c7†)^n7 |vac>.
    Operators act on the ket right-to-left. Thus reversed occupied-mode
    application fixes every charged-sector phase, including species parity.
    """
    coords = [tuple(c) for c in raw_coordinates]
    expected = coordinate_grid(L)
    need(len(coords) == len(set(coords)) == len(expected) and set(coords) == set(expected), 'Complete raw coordinate bijection')
    D = len(coords)
    T = np.zeros((D, D), complex)
    vac = coords.index((0, 0))
    for column, (u, d) in enumerate(expected):
        v = np.zeros(D, complex)
        v[vac] = 1.
        occupied = [j for j in range(L) if u >> j & 1] + [L+j for j in range(L) if d >> j & 1]
        for mode in reversed(occupied):
            v = raw_C[mode].conj().T @ v
        nonzero = np.flatnonzero(np.abs(v) > 1e-12)
        need(len(nonzero) == 1 and coords[nonzero[0]] == (u, d) and abs(abs(v[nonzero[0]])-1) <= 1e-12, 'Unique normalized analytic creation-string ket')
        T[:, column] = v
    need(maxabs(T.conj().T@T-np.eye(D)) <= 1e-12, 'Signed permutation is unitary')
    return T


def keyed(rows, key, expected):
    values = [row[key] for row in rows]
    need(len(values) == len(set(values)) == len(expected) and set(values) == set(expected), 'Complete unique ' + key)
    return dict(zip(values, rows))


def block_grid(L=4):
    return [(N, up, N-up) for N in (3, 4, 5) for up in range(L+1) if 0 <= N-up <= L]


def block_name(N, up, down):
    return f'N{N}_up{up}_dn{down}'


def clusters(records, diameter=1e-8):
    records = sorted(records, key=lambda r: (r['energy'], r['block_key'], r['eigen_index0']))
    groups = []
    for r in records:
        if not groups or r['energy']-groups[-1][0]['energy'] > diameter:
            groups.append([])
        groups[-1].append(r)
    out = []
    for group in groups:
        e = [r['energy'] for r in group]
        out.append(dict(energy=math.fsum(e)/len(e), energy_min=min(e), energy_max=max(e),
                        members=[dict(block_key=r['block_key'], eigen_index0=r['eigen_index0']) for r in group]))
    return out


def cluster_qualified(groups, separation=1e-6):
    return all(b['energy_min']-a['energy_max'] >= separation for a, b in zip(groups, groups[1:]))


def visibility(weight):
    return 'INVISIBLE' if weight <= 1e-12 else ('VISIBLE' if weight >= 1e-9 else 'UNRESOLVED')


def weight_qualified(saved, reconstructed):
    return (math.isfinite(saved) and math.isfinite(reconstructed)
            and saved >= -1e-12 and reconstructed >= -1e-12
            and visibility(reconstructed) != 'UNRESOLVED'
            and visibility(saved) == visibility(reconstructed))


def eigensystem_errors(H, E, V):
    H, E, V = np.asarray(H), np.asarray(E), np.asarray(V)
    need(H.ndim == 2 and H.shape[0] == H.shape[1] and V.shape == H.shape and E.shape == (len(H),), 'Complete eigensystem dimensions')
    need(np.all(np.isfinite(H)) and np.all(np.isfinite(E)) and np.all(np.isfinite(V)) and np.all(np.diff(E) >= 0), 'Finite ascending eigensystem')
    return dict(Hermiticity=maxabs(H-H.conj().T),
                orthogonality=maxabs(V.conj().T@V-np.eye(len(E))),
                residual=maxabs(H@V-V*E))


def complete_pole_measure(groups, vectors, S, dg, E0, final_N, mu):
    """Full projector trace, including zeros; no visible-state filtering."""
    need(S.ndim == 2 and S.shape[1] == dg and dg > 0, 'Full uniform ground-ensemble source')
    poles = []
    for cluster in groups:
        W = np.column_stack([vectors[m['block_key'], m['eigen_index0']] for m in cluster['members']])
        weight = float(np.sum(np.abs(W.conj().T@S)**2)/dg)
        cost = cluster['energy']-E0
        poles.append(dict(energy=cluster['energy'], energy_physical=cluster['energy']+final_N*mu,
                          cost=cost, omega=cost if final_N == 5 else -cost, weight=weight,
                          multiplicity=len(cluster['members']), visibility=visibility(weight)))
    return poles


def measure_errors(saved, computed):
    need(len(saved) == len(computed), 'Every charged cluster including invisible poles retained')
    out = []
    for a, b in zip(saved, computed):
        need(a['multiplicity'] == b['multiplicity'], 'Full charged pole multiplicity')
        out.append(max(abs(a[k]-b[k]) for k in ('energy', 'energy_physical', 'cost', 'omega', 'weight')))
    return out


def channel_grid(boundary, L=4):
    theta = 0. if boundary == 'pbc' else math.pi
    return [dict(channel_id=f'{kind}{j}_{spin}_{action}', kind=kind, index=j, spin=spin, action=action,
                 k=(2*math.pi*j+theta)/L if kind == 'k' else 0.)
            for kind in ('site', 'k') for j in range(L) for spin in ('up', 'down') for action in ('remove', 'add')]


def mode_operator(operators, channel, L=4):
    offset = 0 if channel['spin'] == 'up' else L
    if channel['kind'] == 'site':
        c = operators[offset+channel['index']]
    else:
        c = sum((complex(math.cos(-channel['k']*j), math.sin(-channel['k']*j))*operators[offset+j] for j in range(L)), np.zeros_like(operators[0]))/math.sqrt(L)
    return c.conj().T if channel['action'] == 'add' else c


def audit_quspin_bridge(audit, raw, canonical_C):
    b = raw['metadata']['quspin_bridge']
    raw_C = [sparse(c) for c in b['raw_site_annihilators']]
    need(len(raw_C) == 8, 'All eight raw QuSpin charged operators')
    states = b['raw_basis_states']
    need(len(states) == 256 and all(type(s) is int for s in states) and set(states) == set(range(256)), 'Unreduced full QuSpin integer basis')
    # Standard spinful integer encoding: left up sites, right down sites;
    # site0 is the most significant bit within each species. This checks the
    # number-derived coordinates independently, rather than trusting the map.
    def reverse4(n):
        return sum(((n >> j) & 1) << (3-j) for j in range(4))
    coords = coordinate_pairs(b['raw_coordinates'])
    expected_coords = [(reverse4(s >> 4), reverse4(s & 15)) for s in states]
    audit.check('QuSpin raw integer/number coordinate bridge', coords == expected_coords)
    for mode, C in enumerate(raw_C):
        number = C.conj().T@C
        expected = np.array([float((u if mode < 4 else d) >> (mode % 4) & 1) for u, d in coords])
        audit.check('QuSpin own raw number/operator coordinates mode' + str(mode), close(number, np.diag(expected), 1e-10))
    audit.qualified('QuSpin full raw CAR including crossspin', car_error(raw_C), GATES['car_abs'])
    T = creation_bridge(raw_C, coords)
    mapping = np.argmax(np.abs(T), axis=0).tolist()
    phases = T[mapping, np.arange(256)]
    audit.check('QuSpin analytic canonical-to-raw indices', b['canonical_to_raw'] == mapping)
    audit.check('QuSpin analytic canonical-to-raw phases', close(phases, np.asarray(numeric_list(b['phases_real']))+1j*np.asarray(numeric_list(b['phases_imag'])), 1e-12))
    audit.check('QuSpin declared nonfitted bridge', b['derivation'] == 'vacuum plus canonical creation strings, descending application order; not fitted to target')
    for mode, (raw_op, expected) in enumerate(zip(raw_C, canonical_C)):
        audit.qualified('QuSpin charged bridge canonical mode' + str(mode), maxabs(T.conj().T@raw_op@T-expected), GATES['matrix_abs'])
    return T


def audit_case(audit, rawcase, operators, label, qbridge=None):
    case_id, U, boundary, mu = next(r for r in CASE_GRID if r[0] == rawcase['case_id'])
    prefix = label + '/' + case_id + '/'
    need(rawcase['U'] == U and rawcase['mu'] == mu and rawcase['boundary'] == boundary, 'Exact fixed case model')
    expected_H = canonical_hamiltonian(4, 1., U, mu, boundary)
    expected_blocks = {block_name(*k): k for k in block_grid()}
    blocks = keyed(rawcase['blocks'], 'block_key', expected_blocks)
    vectors, spectra, ownH = {}, {3: [], 4: [], 5: []}, np.zeros((256, 256), complex)
    coords = coordinate_grid()
    for key, (N, up, down) in expected_blocks.items():
        b = blocks[key]
        need(b['N'] == N and b['Nup'] == up and b['Ndown'] == down, 'Exact charged/spin block identity')
        ids = [i for i, (u, d) in enumerate(coords) if u.bit_count() == up and d.bit_count() == down]
        bc = coordinate_pairs(b['coordinates'])
        need(len(bc) == len(set(bc)) == len(ids) and set(bc) == {coords[i] for i in ids}, 'Complete block coordinates')
        # Arbitrary row order is allowed only through its explicit coordinates.
        ids = [u*16+d for u, d in bc]
        dim = len(ids)
        H, V = matrix(b['Hamiltonian'], dim, dim), matrix(b['eigenvectors'], dim, dim)
        E = np.asarray(numeric_list(b['eigenvalues']), float)
        need(E.shape == (dim,) and np.all(np.isfinite(E)) and np.all(np.diff(E) >= 0), 'Complete ascending eigenspectrum')
        errors = eigensystem_errors(H, E, V)
        audit.qualified(prefix+key+'/Hermiticity', errors['Hermiticity'], GATES['matrix_abs'])
        audit.qualified(prefix+key+'/independent Hubbard H', maxabs(H-expected_H[np.ix_(ids, ids)]), GATES['matrix_abs'])
        audit.qualified(prefix+key+'/orthogonality', errors['orthogonality'], GATES['orthogonality_abs'])
        audit.qualified(prefix+key+'/all eigen residuals', errors['residual'], GATES['eigen_residual_abs'])
        ownH[np.ix_(ids, ids)] = H
        for j, e in enumerate(E):
            v = np.zeros(256, complex)
            v[ids] = V[:, j]
            vectors[key, j] = v
            spectra[N].append(dict(energy=float(e), block_key=key, eigen_index0=j))
    groups = {N: clusters(spectra[N]) for N in (3, 4, 5)}
    saved_spectra = keyed(rawcase['sector_spectra'], 'N', (3, 4, 5))
    for N in (3, 4, 5):
        actual = saved_spectra[N]['clusters']
        need(len(actual) == len(groups[N]), 'Complete full-sector energy clusters')
        for j, (left, right) in enumerate(zip(actual, groups[N])):
            audit.check(prefix+f'N{N} cluster{j} members', left['members'] == right['members'])
            audit.check(prefix+f'N{N} cluster{j} energy data', all(abs(left[f]-right[f]) <= 1e-10 for f in ('energy', 'energy_min', 'energy_max')))
        audit.numerical.append(dict(name=prefix+f'N{N} unambiguous full-diameter clusters', passed=cluster_qualified(groups[N]), detail='Distinct adjacent cluster gap must be at least1e-6; all clusters diameter at most1e-8.'))
    gc = groups[4][0]
    G = np.column_stack([vectors[m['block_key'], m['eigen_index0']] for m in gc['members']])
    dg, E0 = G.shape[1], gc['energy']
    saved_ground = rawcase['ground']
    need(saved_ground['dimension'] == dg and dg > 0, 'Complete global canonical ground degeneracy')
    SG = matrix(saved_ground['vectors'], 256, dg)
    Pg = G@G.conj().T
    audit.check(prefix+'saved own ground projector', close(SG@SG.conj().T, Pg, 1e-10))
    audit.check(prefix+'saved own ground orthonormality', close(SG.conj().T@SG, np.eye(dg), 1e-10))
    audit.check(prefix+'saved own ground energies', close(numeric_list(saved_ground['energies']), [blocks[m['block_key']]['eigenvalues'][m['eigen_index0']] for m in gc['members']], 1e-10))
    audit.check(prefix+'saved ground K/physical energy', abs(saved_ground['energy']-E0) <= 1e-10 and abs(saved_ground['energy_physical']-E0-4*mu) <= 1e-10)
    if qbridge is not None:
        qH = sparse(rawcase['raw_quspin_hamiltonian'])
        keep = [i for i, (u, d) in enumerate(coords) if u.bit_count()+d.bit_count() in (3, 4, 5)]
        audit.qualified(prefix+'raw QuSpin full Hamiltonian bridge', maxabs((qbridge.conj().T@qH@qbridge)[np.ix_(keep, keep)]-ownH[np.ix_(keep, keep)]), GATES['matrix_abs'])
        # Raw full H also covers all other sectors; no H-only bridge can certify C.
        audit.qualified(prefix+'independent all-sector raw QuSpin H', maxabs(qbridge.conj().T@qH@qbridge-expected_H), GATES['matrix_abs'])
    channels = keyed(rawcase['channels'], 'channel_id', [c['channel_id'] for c in channel_grid(boundary)])
    computed = {}
    for desc in channel_grid(boundary):
        cid = desc['channel_id']
        row = channels[cid]
        audit.check(prefix+cid+'/exact channel', all(row[k] == v if k != 'k' else abs(row[k]-v) <= EPS64 for k, v in desc.items()))
        N = 5 if desc['action'] == 'add' else 3
        need(row['final_N'] == N, 'Charged source final sector')
        O = mode_operator(operators, desc)
        # The saved ground columns need not equal independently assembled columns.
        # Audit source columns against that backend's own raw ground columns.
        savedS = matrix(row['source_vectors'], 256, dg)
        audit.check(prefix+cid+'/own raw source vector columns', close(savedS, O@SG, 1e-10))
        S = O@G
        density = S@S.conj().T/dg
        savedpoles = row['poles']
        poles = complete_pole_measure(groups[N], vectors, S, dg, E0, N, mu)
        errors = measure_errors(savedpoles, poles)
        for j, (saved, pole, error) in enumerate(zip(savedpoles, poles, errors)):
            audit.check(prefix+cid+f'/pole{j} own raw formula', error <= 1e-10)
            weight = pole['weight']
            saved_label = visibility(saved['weight'])
            audit.numerical.append(dict(name=prefix+cid+f'/pole{j} visibility qualification',
                                        passed=weight_qualified(saved['weight'], weight),
                                        weight=weight, saved_weight=saved['weight'], visibility=pole['visibility'], saved_visibility=saved_label))
        c = mode_operator(operators, dict(desc, action='remove'))
        occupation = float(np.sum(np.abs(c@G)**2)/dg)
        zero = float(np.sum(np.abs(S)**2)/dg)
        first = float(np.vdot(S, ownH@S-E0*S).real/dg)
        spectral_zero = math.fsum(p['weight'] for p in poles)
        spectral_first = math.fsum(p['weight']*p['cost'] for p in poles)
        audit.check(prefix+cid+'/own occupation formula', abs(row['occupation']-occupation) <= GATES['occupation_abs'])
        audit.check(prefix+cid+'/own zeroth formula', abs(row['zeroth_moment']-zero) <= GATES['zeroth_moment_abs'])
        audit.check(prefix+cid+'/own first formulas', abs(row['cost_first_moment']-first) <= GATES['first_moment_abs'] and abs(row['signed_first_moment']-(first if N == 5 else -first)) <= GATES['first_moment_abs'])
        audit.check(prefix+cid+'/thermodynamic threshold from every final state', abs(row['charged_threshold']-poles[0]['cost']) <= GATES['matrix_abs'])
        audit.qualified(prefix+cid+'/zeroth sumrule', abs(spectral_zero-zero), GATES['zeroth_moment_abs'])
        audit.qualified(prefix+cid+'/first sumrule', abs(spectral_first-first), GATES['first_moment_abs'])
        audit.qualified(prefix+cid+'/occupation/addrem normalization', abs(zero-(occupation if N == 3 else 1.-occupation)), GATES['occupation_abs'])
        visible = [p for p in poles if p['visibility'] == 'VISIBLE']
        computed[cid] = dict(descriptor=desc, final_N=N, source_density=density,
                             occupation=occupation, zeroth_moment=zero, cost_first_moment=first,
                             charged_threshold=poles[0]['cost'], lowest_visible_cost=visible[0]['cost'] if visible else None,
                             lowest_charged_visibility=poles[0]['visibility'], poles=poles)
    for desc in channel_grid(boundary):
        if desc['action'] == 'remove':
            a = computed[desc['channel_id']]
            b = computed[desc['channel_id'].removesuffix('remove')+'add']
            audit.qualified(prefix+desc['channel_id']+'/paired zeroth sum one', abs(a['zeroth_moment']+b['zeroth_moment']-1.), GATES['zeroth_moment_abs'])
    return dict(case_id=case_id, blocks=blocks, groups=groups, ground_projector=Pg,
                full_covered_H=ownH, ground_dimension=dg, ground_energy=E0, channels=computed)


def audit_backend(audit, raw, label, source, cfg, freeze_sha):
    meta = raw['metadata']
    audit.check(label+'/schema and backend', raw['schema'] == meta['schema'] == RAW_SCHEMA and meta['backend'] == label)
    audit.check(label+'/canonical release and literal ID', meta['release'] == RELEASE and meta['catalogue_id'] == 'X13')
    audit.check(label+'/fixed ensemble and Fock order', meta['initial_ensemble'] == cfg['initial_ensemble'] and meta['fock_order'] == cfg['fock_order'])
    audit.check(label+'/native and complete runtime source pins', meta['native_sha256'] == NATIVE_PINS['analysis/compute_s2_ed.jl'] and meta['manifest_sha256'] == sha(source/'manifest.toml') and meta['source_freeze_sha256'] == freeze_sha and meta['source_sha256'] == {n: sha(source/n) for n in SOURCE_ROLES})
    versions = {'julia': '1.12.1'} if label == 'native_ed_new_adapter' else {'quspin': '1.0.1', 'numpy': '2.4.6', 'scipy': '1.17.1'}
    audit.check(label+'/actual pinned runtime versions', meta['versions'] == versions)
    finished = raw['completion']['finished_unix']
    audit.check(label+'/complete finite raw receipt before cutoff', raw['completion']['status'] == 'COMPLETE' and math.isfinite(finished) and 0 < finished < CUTOFF and meta.get('scientific_sampling') is False)
    need(len(raw['site_annihilators']) == 8, 'Every site and spin annihilator')
    actual_C = [sparse(c) for c in raw['site_annihilators']]
    canonical_C = [annihilator(4, mode) for mode in range(8)]
    for mode, (C, expected) in enumerate(zip(actual_C, canonical_C)):
        audit.qualified(label+'/canonical bit-parity site C'+str(mode), maxabs(C-expected), GATES['matrix_abs'])
    audit.qualified(label+'/full canonical same/crossspin CAR', car_error(actual_C), GATES['car_abs'])
    T = audit_quspin_bridge(audit, raw, canonical_C) if label == 'independent_quspin' else None
    cases = keyed(raw['cases'], 'case_id', [r[0] for r in CASE_GRID])
    derived = {}
    for cid, *_ in CASE_GRID:
        derived[cid] = audit_case(audit, cases[cid], actual_C, label, T)
    return derived


def compare_backends(audit, target, reference):
    rows = []
    for cid, *_ in CASE_GRID:
        a, b = target[cid], reference[cid]
        audit.compare(cid+'/cross complete covered Hamiltonian', maxabs(a['full_covered_H']-b['full_covered_H']), GATES['matrix_abs'])
        audit.compare(cid+'/cross ground dimension', float(a['ground_dimension'] != b['ground_dimension']), 0.)
        audit.compare(cid+'/cross global ground projector', maxabs(a['ground_projector']-b['ground_projector']), GATES['cross_projector_abs'])
        for N in (3, 4, 5):
            ea = sorted(e for block in a['blocks'].values() if block['N'] == N for e in block['eigenvalues'])
            eb = sorted(e for block in b['blocks'].values() if block['N'] == N for e in block['eigenvalues'])
            audit.compare(cid+f'/cross complete N{N} eigenenergies', maxabs(np.asarray(ea)-eb), GATES['cross_energy_abs'])
        for channel_id, x in a['channels'].items():
            y = b['channels'][channel_id]
            prefix = cid+'/'+channel_id+'/'
            audit.compare(prefix+'cross source projector covariance', maxabs(x['source_density']-y['source_density']), GATES['cross_source_density_abs'])
            same_count = len(x['poles']) == len(y['poles'])
            audit.compare(prefix+'cross complete pole count', float(not same_count), 0.)
            for j, (p, q) in enumerate(zip(x['poles'], y['poles'])):
                audit.compare(prefix+f'pole{j} cross multiplicity', float(p['multiplicity'] != q['multiplicity']), 0.)
                audit.compare(prefix+f'pole{j} cross K/cost/signed energy', max(abs(p[k]-q[k]) for k in ('energy', 'cost', 'omega')), GATES['cross_energy_abs'])
                audit.compare(prefix+f'pole{j} cross weight', abs(p['weight']-q['weight']), GATES['cross_weight_abs'])
                audit.compare(prefix+f'pole{j} cross visibility', float(p['visibility'] != q['visibility']), 0.)
            audit.compare(prefix+'cross thermodynamic threshold', abs(x['charged_threshold']-y['charged_threshold']), GATES['cross_threshold_abs'])
            va, vb = x['lowest_visible_cost'], y['lowest_visible_cost']
            audit.compare(prefix+'cross first visible threshold', float((va is None) != (vb is None)) if va is None or vb is None else abs(va-vb), 0. if va is None or vb is None else GATES['cross_threshold_abs'])
            rows.append(dict(case_id=cid, **x['descriptor'],
                             target={k: v for k, v in x.items() if k not in ('source_density', 'descriptor')},
                             reference={k: v for k, v in y.items() if k not in ('source_density', 'descriptor')},
                             hidden_lowest=x['lowest_charged_visibility'] == 'INVISIBLE'))
    return rows


def source_closure(audit, args):
    source = args.release_dir.resolve()
    freeze = audit.read(args.freeze)
    need(re.fullmatch('[0-9a-f]{64}', args.freeze_sha256) and re.fullmatch('[0-9a-f]{64}', args.source_review_sha256), 'Explicit externally approved freeze/review anchors')
    need(sha(args.freeze) == args.freeze_sha256, 'Approved external freeze SHA')
    need(freeze.get('catalogue_id') == 'X13' and freeze.get('release') == RELEASE, 'Exact frozen release/catalogue identity')
    if 'release_name' in freeze:
        need(freeze['release_name'] == freeze['release'], 'Freeze release aliases agree')
    files = freeze['files']
    need(isinstance(files, dict) and set(SOURCE_ROLES) <= set(files) and 'SHA256SUMS' in files and 'source_freeze.json' not in files, 'Noncircular complete frozen inventory')
    actual = {p.relative_to(source).as_posix() for p in source.rglob('*') if p.is_file() and p.name != 'source_freeze.json'}
    need(actual == set(files), 'Exact release inventory, no hidden/missing files')
    need(not any(p.is_symlink() for p in source.rglob('*')), 'No release symlinks')
    aliases = set()
    for name, expected in files.items():
        rel = safe_relative(name)
        need(name.casefold() not in aliases and re.fullmatch('[0-9a-f]{64}', expected), 'Unique valid frozen hash role')
        aliases.add(name.casefold())
        need(sha(audit.file(source/rel)) == expected, 'Frozen source changed: '+name)
    need(audit.file(source/'source_freeze.json').read_bytes() == Path(args.freeze).read_bytes(), 'Executed administrative source freeze exact')
    checksums = (source/'SHA256SUMS').read_bytes()
    need(b'\r' not in checksums, 'LF-only checksum inventory')
    listed = {}
    for line in checksums.decode('utf-8').splitlines():
        m = re.fullmatch(r'([0-9a-f]{64})  (.+)', line)
        need(m is not None and m[2] not in listed, 'Canonical SHA256SUMS')
        listed[m[2]] = m[1]
    need(listed == {k: v for k, v in files.items() if k != 'SHA256SUMS'}, 'Checksum manifest covers all source bytes exactly')
    review = audit.read(args.source_review)
    need(sha(args.source_review) == args.source_review_sha256 == freeze['independent_source_review_sha256'], 'Approved independent source review anchor')
    need(PurePosixPath(freeze['independent_source_review']).name == args.source_review.name, 'Exact approved source review identity')
    if 'peer_review_sha256' in freeze:
        need(freeze['peer_review_sha256'] == freeze['independent_source_review_sha256'], 'Freeze review aliases agree')
    if review.get('schema') == 'X13_INDEPENDENT_LOCAL_CANDIDATE_PEER_REVIEW_V1':
        need(review.get('status') == 'NO_MATERIAL_BLOCKER_IN_LOCAL_SOURCE_PENDING_ROOT_APPROVAL_AND_REQUIRED_PREFLIGHT'
             and review.get('canonical_release') == RELEASE and review.get('catalogue_id') == 'X13'
             and review.get('scientific_execution') is False
             and review['source_file_sha256'] == listed
             and all(r['passed'] is True for r in review['checks']),
             'Root-anchored source-only peer review binds complete snapshot; no execution approval inferred')
    else:
        need(review.get('status') == 'APPROVED_FOR_SOURCE_FREEZE_ONLY' and review['files'] == listed,
             'Explicit source-freeze-only review binds complete pre-checksum snapshot')
    archive = audit.file(args.source_archive)
    need(sha(archive) == freeze['archive_sha256'] and archive.stat().st_size == freeze['archive_bytes'], 'Exact approved source archive bytes/length')
    need(PurePosixPath(freeze['archive']).name == archive.name, 'Exact approved source archive identity')
    verify_source_archive(archive, files)
    for name, expected in NATIVE_PINS.items():
        need(sha(audit.file(args.native_source_dir/safe_relative(name))) == expected, 'Pinned native source/environment '+name)
    catalogue = audit.read(args.catalogue)
    need(sha(args.catalogue) == CATALOGUE_SHA, 'Original100-question catalogue exact SHA')
    # Catalogue shape is deliberately not inferred from status ledgers.
    items = catalogue if isinstance(catalogue, list) else catalogue.get('benchmarks', catalogue.get('ideas', []))
    row = [r for r in items if r.get('id') == 'X13']
    need(len(row) == 1 and row[0]['acceptance'] == 'Match both thresholds and visibility weights; do not require an invisible state to create a spectral peak.', 'Literal X13 question and acceptance')
    cfg = audit.read(source/'manifest.toml')
    need(cfg['schema'] == 'X13_CHARGED_SPECTRAL_PROTOCOL_V1' and cfg['catalogue_id'] == 'X13' and cfg['release'] == RELEASE, 'Fixed protocol identity')
    need(cfg['gates'] == GATES, 'All fixed numerical/visibility gates unchanged')
    need(cfg['L'] == 4 and cfg['t'] == 1. and cfg['initial_N'] == 4 and cfg['total_N'] == [3, 4, 5], 'Fixed model and complete charged sectors')
    need([(c['case_id'], c['U'], c['boundary'], c['mu']) for c in cfg['cases']] == CASE_GRID, 'Exact four-case model grid')
    need(cfg['native_sha256'] == NATIVE_PINS['analysis/compute_s2_ed.jl'] and cfg['native_relative_path'] == 'analysis/compute_s2_ed.jl', 'Unchanged native API source')
    need(cfg['science_cutoff_unix'] == CUTOFF and cfg['wall_cap_seconds'] == 1200 and cfg['cpus'] == 1 and cfg['memory_mb'] == 3500, 'Fixed allocation and whole-science deadline')
    need(cfg['initial_ensemble'] == 'Uniform canonical complete global N4 ground projector Pg/dg over every spin block' and cfg['fock_order'] == 'index0=up_mask*16+down_mask; site j is bit j; modes up0..up3,down0..down3', 'Exact initial ensemble and canonical mode order')
    provenance = audit.read(source/'PROVENANCE.json')
    need(provenance['schema'] == 'X13_ADAPTER_PROVENANCE_V1' and provenance['catalogue_id'] == 'X13' and provenance['catalogue_sha256'] == CATALOGUE_SHA and provenance['native_source_sha256'] == NATIVE_PINS, 'Preserved catalogue/native source provenance')
    for p in source.iterdir():
        if p.suffix in ('.py', '.jl', '.sbatch', '.toml', '.sh'):
            need(b'\r' not in p.read_bytes(), 'Raw LF-only executable/config '+p.name)
    audit.check('approved complete source/archive/review/native/catalogue closure before replay', True)
    return cfg, freeze


def archive_role(name, files):
    """Root's frozen archive is flat; allow one exact unambiguous layout.

    The source-freeze/archive hash was checked before this parser is called.
    No paths are extracted, and extra roles/aliases/mixed layouts are refused.
    """
    p = safe_relative(name)
    if name in files:
        return name, 'flat'
    if len(p.parts) > 1 and p.parts[0] == RELEASE:
        role = '/'.join(p.parts[1:])
        need(role in files, 'Unexpected canonical-prefixed archive role')
        return role, 'canonical_prefix'
    raise ValueError('Archive member outside exact frozen source roles')


def verify_source_archive(archive, files):
    with tarfile.open(archive, 'r:gz') as tar:
        seen, layouts = set(), set()
        for m in tar.getmembers():
            need(m.isfile() and not m.issym() and not m.islnk(), 'Only regular source archive files')
            name, layout = archive_role(m.name, files)
            need(name not in seen, 'Unique source archive roles')
            seen.add(name)
            layouts.add(layout)
            need(len(layouts) == 1, 'Mixed flat/prefixed archive layout refused')
            need(hashlib.sha256(tar.extractfile(m).read()).hexdigest() == files[name], 'Archive source member hash')
        need(seen == set(files), 'Complete source archive closure')


def duration_seconds(value):
    days = 0
    if '-' in value:
        d, value = value.split('-', 1)
        days = int(d)
    fields = [float(x) for x in value.split(':')]
    need(1 <= len(fields) <= 3 and all(x >= 0 for x in fields), 'Finite scheduler duration')
    return days*86400 + sum(v*60**i for i, v in enumerate(reversed(fields)))


def scheduler_epoch(value):
    """This dated ISAAC campaign's sacct clock is EDT (UTC−4).

    Current verified scheduler receipts are naive local timestamps. Only the
    campaign's Oct22/23 dates are accepted; an explicit timezone is honored.
    The assumption is reported, and never inferred from scientific timestamps.
    """
    value = dt.datetime.fromisoformat(value.replace('Z', '+00:00'))
    if value.tzinfo is None:
        need(value.date() in (dt.date(2026, 10, 22), dt.date(2026, 10, 23)), 'Dated campaign Slurm local clock')
        value = value.replace(tzinfo=dt.timezone(dt.timedelta(hours=-4)))
    return value.timestamp()


def bind_execution_clock(ex, scheduler_rows):
    need(ex['actual_batch_start_unix'] <= ex['worker_start_unix'] <= ex['processes'][0]['start_unix'], 'Ordered batch/worker/first science start')
    times = [ex['actual_batch_start_unix'], ex['worker_start_unix'], ex['finished_unix']]
    times.extend(t for p in ex['processes'] for t in (p['start_unix'], p['end_unix']))
    need(all(type(x) in (int, float) and math.isfinite(x) for x in times), 'Finite exact execution epoch fields')
    for row in scheduler_rows:
        start, end = scheduler_epoch(row['Start']), scheduler_epoch(row['End'])
        need(start < end and abs((end-start)-duration_seconds(row['Elapsed'])) <= 1, 'Scheduler start/end/elapsed agreement')
        need(all(start-1 <= t <= end+1 for t in times), 'Every batch/worker/process/finish epoch inside actual scheduler interval (1s rounding)')
    return True


def transfer_closure(audit, args, freeze):
    jobdir = args.job_dir.resolve()
    record = audit.read(jobdir/'transfer.json')
    need(str(record['job_id']) == args.job_id and record.get('terminal_snapshot') is True and record.get('two_pass_inventory_and_bytes_verified') is True and record.get('remote_writes') is False and record.get('binaries_copied') is False and not record.get('missing_paths'), 'Complete stable read-only terminal transfer')
    if 'transfer_status' in record:
        need(record['transfer_status'] == 'COMPLETE', 'Complete original transfer publication')
    seen, aliases = set(), set()
    for row in record['files']:
        rel = safe_relative(row['local'])
        need(row['status'] == 'RETRIEVED' and str(rel).casefold() not in aliases, 'Unique retrieved text receipt')
        aliases.add(str(rel).casefold())
        p = audit.file(jobdir/rel)
        need(p.is_relative_to(jobdir), 'Transferred file stays in job directory')
        need(p.stat().st_size == row['bytes'] and sha(p) == row['sha256'], 'Exact original raw transfer '+str(rel))
        seen.add(str(rel))
    required = {'results/target.toml', 'results/reference.json', 'results/comparison.json', 'results/execution.json', 'results/preflight.json'}
    required.update('results/'+name+suffix for name in ('target', 'reference', 'comparator') for suffix in ('.stdout.txt', '.stderr.txt', '.resources.txt', '.process.json'))
    required.update({'results/scontrol.txt', 'results/process_exit_codes.txt'})
    required.update('source/'+n for n in freeze['files'])
    required.add('source/source_freeze.json')
    required.update('native_source/'+n for n in NATIVE_PINS)
    need(required <= seen, 'Every source/native/raw role present in original transfer')
    for name, h in freeze['files'].items():
        need(sha(jobdir/'source'/name) == h, 'Actual executed retrieved source '+name)
    need(sha(jobdir/'source/source_freeze.json') == args.freeze_sha256, 'Actual executed external source freeze')
    for name, h in NATIVE_PINS.items():
        need(sha(jobdir/'native_source'/name) == h, 'Actual retrieved native environment '+name)
    scheduler_raw = audit.file(jobdir/'scheduler.txt').read_bytes()
    # New receipt scheduler bytes use LF. Prior administrative CRLF receipts
    # remain faithful through exact transfer hashes; no original is rewritten.
    need(scheduler_raw.decode('utf-8-sig').replace('\r\n', '\n') == record['sacct'].replace('\r\n', '\n'), 'Scheduler receipt bound to original terminal snapshot')
    scontrol = record['scontrol']
    control_raw = audit.file(jobdir/'scontrol.txt').read_bytes().decode('utf-8-sig').replace('\r\n', '\n')
    expected_control = 'exit_code='+str(scontrol['exit_code'])+'\n'+scontrol['stdout']+scontrol['stderr']
    need(control_raw == expected_control.replace('\r\n', '\n'), 'scontrol exact supplementary receipt')
    rows = list(csv.DictReader(io.StringIO(record['sacct']), delimiter='|'))
    main = [r for r in rows if r.get('JobIDRaw') == args.job_id]
    batch = [r for r in rows if r.get('JobIDRaw') == args.job_id+'.batch']
    need(len(main) == len(batch) == 1, 'Unique actual job main/batch scheduler rows')
    m, b = main[0], batch[0]
    need(m['State'] in ('COMPLETED', 'FAILED') and b['State'] in ('COMPLETED', 'FAILED') and m['ExitCode'] in ('0:0', '1:0', '2:0') and b['ExitCode'] == m['ExitCode'], 'Terminal exact comparison returncode, no timeout/preemption')
    need(m['State'] == b['State'] == ('COMPLETED' if m['ExitCode'] == '0:0' else 'FAILED'), 'Scheduler state matches exact comparator outcome')
    need(m['ReqCPUS'] == m['AllocCPUS'] == '1' and m['ReqMem'] == '3500M' and duration_seconds(m['Timelimit']) == 1200 and 0 < duration_seconds(m['Elapsed']) <= 1201, 'Exact bounded one-CPU 3500MB20minute allocation')
    for beforeafter in ('scheduler_before', 'scheduler_after'):
        snap = record[beforeafter]
        need(snap.get('target_queue_lines') == [] and all(snap['main'].get(k) == m.get(k) for k in ('JobIDRaw', 'State', 'ExitCode', 'Start', 'End')), 'Stable terminal before/after scheduler identity')
    audit.check('all raw/source/native transfer bytes and terminal scheduler closure', True)
    return record, (m, b)


def verify_actual_start_admission(receipt, cfg):
    """Bind the six-second TERM reserve to the already frozen manifest.

    candidate.sbatch separately uses a five-second KILL grace. These are
    different values; neither changes the preregistered 1200-second cap.
    """
    need(cfg['termination_guard_seconds'] == 6 and cfg['wall_cap_seconds'] == 1200
         and cfg['science_cutoff_unix'] == CUTOFF, 'Exact frozen TERM guard/cap/cutoff configuration')
    need(receipt['science_cutoff_unix'] == CUTOFF
         and receipt['cap_end_unix'] == receipt['actual_batch_start_unix']+1200 <= CUTOFF
         and receipt['termination_guard_seconds'] == cfg['termination_guard_seconds'],
         'Actual-start fullcap deadline admission with frozen six-second TERM reserve')
    return cfg['termination_guard_seconds']


def preliminary_execution(audit, args, cfg, freeze, scheduler):
    """Runner-specific field contract is explicit; no log prose counts as proof."""
    result = args.job_dir/'results'
    ex = audit.read(result/'execution.json')
    pre = audit.read(result/'preflight.json')
    # Completed source/physics/compare processes are required even if the last
    # comparator deliberately returns1(FAIL) or2(UNRESOLVED).
    need(str(ex['job_id']) == args.job_id and str(pre['job_id']) == args.job_id, 'Executed job identity')
    need(ex.get('schema') == pre.get('schema') == 'X13_EXECUTION_RECEIPT_V1', 'Explicit complete X13 execution schema')
    need(ex.get('scientific_sampling') is False and ex.get('catalogue_id') == 'X13', 'Finite exact X13 execution, no MC')
    for receipt in (pre, ex):
        need(receipt['source_freeze_sha256'] == args.freeze_sha256 and receipt['native_source_sha256'] == NATIVE_PINS and receipt['source_sha256'] == freeze['files'], 'Actual complete executed source/native inventory')
        need(receipt['release'] == RELEASE and receipt['catalogue_id'] == 'X13' and receipt['cpus'] == 1 and receipt['memory_mb'] == 3500, 'Actual source/resource identity')
        need(receipt['threads'] == {k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'JULIA_NUM_THREADS')}, 'Actual single-thread environment receipt')
        verify_actual_start_admission(receipt, cfg)
    need(pre['status'] == 'STARTED' and pre['processes'] == [] and pre['process_exit_codes'] == {}, 'Pre-science source preflight receipt')
    need(ex['status'] == 'PROCESSES_COMPLETE' and ex['source_unchanged'] is True and ex['native_unchanged'] is True, 'Complete unchanged-source execution receipt')
    for key in ('actual_batch_start_unix', 'worker_start_unix', 'cap_end_unix', 'science_cutoff_unix', 'process_group_id'):
        need(pre[key] == ex[key], 'Preflight/execution immutable timing/process group identity')
    need(pre['active_queue_snapshot'] == ex['active_queue_snapshot'], 'Immutable actual-start whole-user queue receipt')
    queue = [line.split('|') for line in pre['active_queue_snapshot'].splitlines() if line.strip()]
    need(1 <= len(queue) <= 2 and all(len(r) == 3 and r[1] in ('RUNNING', 'COMPLETING', 'CONFIGURING', 'SUSPENDED', 'RESIZING') and r[2].isdigit() and int(r[2]) > 0 for r in queue), 'Actual-start active allocation roles')
    need(len({r[0] for r in queue}) == len(queue) and args.job_id in {r[0] for r in queue} and sum(int(r[2]) for r in queue) <= 2, 'At most two whole-user active jobs/CPUs including this job')
    processes = ex['processes']
    need([p['name'] for p in processes] == ['target', 'reference', 'comparator'], 'Exact sequential native/reference/comparison roles')
    need(all(p['exit_code'] == 0 for p in processes[:2]) and processes[2]['exit_code'] in (0, 1, 2), 'Both complete backends; comparator independent outcome')
    need(all(p['process_group_inherited'] is True and ex['actual_batch_start_unix'] <= p['start_unix'] < p['end_unix'] <= ex['cap_end_unix']-cfg['termination_guard_seconds'] for p in processes), 'All actual processes before cutoff with inherited bounded process group')
    need(all(x['end_unix'] <= y['start_unix'] for x, y in zip(processes, processes[1:])), 'Serial process timing')
    need(ex['process_exit_codes'] == {r['name']: r['exit_code'] for r in processes}, 'Exact process exitcode vector')
    need(processes[-1]['end_unix'] <= ex['finished_unix'] <= ex['cap_end_unix']-cfg['termination_guard_seconds'], 'Final receipt before outer termination reserve')
    need(scheduler[0]['ExitCode'] == str(processes[2]['exit_code'])+':0', 'Scheduler exact comparator exit')
    bind_execution_clock(ex, scheduler)
    required_outputs = {'preflight.json', 'scontrol.txt', 'target.toml', 'reference.json', 'comparison.json'}
    required_outputs.update(name+suffix for name in ('target', 'reference', 'comparator') for suffix in ('.stdout.txt', '.stderr.txt', '.resources.txt', '.process.json'))
    need(set(ex['output_sha256']) == required_outputs, 'Exact complete17 execution output hash roles')
    for name, h in ex['output_sha256'].items():
        need(sha(audit.file(result/name)) == h, 'Executed raw output hash '+name)
    for process in processes:
        name = process['name']
        need(audit.read(result/(name+'.process.json')) == process, 'Immutable per-process receipt')
        text = audit.file(result/(name+'.resources.txt')).read_text(encoding='utf-8')
        exitcodes = re.findall(r'^\s*Exit status:\s*(\d+)\s*$', text, re.M)
        need(exitcodes == [str(process['exit_code'])], 'GNUtime process exit status')
        walls = re.findall(r'^\s*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)\s*$', text, re.M)
        need(len(walls) == 1 and 0 <= duration_seconds(walls[0]) <= process['end_unix']-process['start_unix']+1, 'GNUtime bounded actual process wall')
        rss = re.findall(r'^\s*Maximum resident set size \(kbytes\):\s*(\d+)\s*$', text, re.M)
        need(len(rss) == 1 and 0 < int(rss[0]) <= 3500*1024, 'GNUtime process memory bound')
    exit_text = audit.file(result/'process_exit_codes.txt').read_bytes()
    need(exit_text == ''.join(name+'_exit='+str(ex['process_exit_codes'][name])+'\n' for name in ('target', 'reference', 'comparator')).encode('ascii'), 'Exact post-execution code text')
    audit.check('actual runtime/execution/fullcap/serial/cutoff closure', True)
    return ex, pre


def status(audit):
    if not audit.integrity() or not all(r['passed'] for r in audit.numerical):
        return 'UNRESOLVED'
    return 'PASS' if all(r['passed'] for r in audit.comparisons) else 'FAIL'


def replay_comparator(audit, args, saved_status):
    # Called only after source_closure, transfer_closure and primary independent
    # reconstruction. The replay uses an empty scratch sibling of the output.
    with tempfile.TemporaryDirectory(prefix='x13_frozen_comparator_', dir=args.output.resolve().parent) as d:
        output = Path(d)/'comparison.json'
        command = [sys.executable, '-B', str(args.release_dir/'compare.py'),
                   '--target', str(args.job_dir/'results/target.toml'),
                   '--reference', str(args.job_dir/'results/reference.json'),
                   '--manifest', str(args.release_dir/'manifest.toml'),
                   '--freeze-sha256', args.freeze_sha256, '--output', str(output)]
        proc = subprocess.run(command, cwd=args.release_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, timeout=180, env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})
        need(output.is_file(), 'Frozen comparator produced fresh output')
        replay = load(output)
        saved = load(args.job_dir/'results/comparison.json')
        need(replay['scientific_status'] == saved_status == saved['scientific_status'] and proc.returncode == {'PASS': 0, 'FAIL': 1, 'UNRESOLVED': 2}[saved_status], 'Frozen replay labels and exitcodes match independent status')
        # Absolute paths differ only if receipt sources moved. Compare all
        # science/check arrays recursively; input paths remain separately bound.
        left = {k: v for k, v in saved.items() if k != 'input_sha256'}
        right = {k: v for k, v in replay.items() if k != 'input_sha256'}
        recursive_roundoff(left, right)
        return dict(performed=True, supplemental_only=True, gates_and_labels_exact=True,
                    metrics_within_64_epsilon=True, saved_sha256=sha(args.job_dir/'results/comparison.json'),
                    replay_sha256=sha(output), returncode=proc.returncode)


def recursive_roundoff(a, b):
    if isinstance(a, dict) or isinstance(b, dict):
        need(isinstance(a, dict) and isinstance(b, dict) and set(a) == set(b), 'Replay exact object keys')
        for key in a:
            recursive_roundoff(a[key], b[key])
    elif isinstance(a, list) or isinstance(b, list):
        need(isinstance(a, list) and isinstance(b, list) and len(a) == len(b), 'Replay exact list shape')
        for x, y in zip(a, b):
            recursive_roundoff(x, y)
    elif isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        need(a == b, 'Replay exact nonnumeric value')
    else:
        need(math.isfinite(a) and math.isfinite(b) and abs(a-b) <= EPS64*max(1., abs(a), abs(b)), 'Replay scalar64epsilon')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--job-dir', type=Path, required=True)
    p.add_argument('--job-id', required=True)
    p.add_argument('--release-dir', type=Path, required=True)
    p.add_argument('--native-source-dir', type=Path)
    p.add_argument('--freeze', type=Path, required=True)
    p.add_argument('--freeze-sha256', required=True, help='Explicit root-approved source-freeze hash; no unfrozen admission')
    p.add_argument('--source-review', type=Path, required=True)
    p.add_argument('--source-review-sha256', required=True)
    p.add_argument('--source-archive', type=Path, required=True)
    p.add_argument('--catalogue', type=Path, default=Path(__file__).resolve().parent/'benchmark_ideas_100.json')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--replay', action='store_true', help='Supplemental frozen comparator, after independent audit; requires candidate runtime dependencies locally')
    args = p.parse_args(argv)
    args.job_dir, args.release_dir = args.job_dir.resolve(), args.release_dir.resolve()
    args.native_source_dir = (args.native_source_dir or args.job_dir/'native_source').resolve()
    need(re.fullmatch(r'[0-9]+', args.job_id), 'Numeric scheduler job ID')
    need(not args.output.exists() and args.output.parent.is_dir(), 'Fresh explicit output in existing writable directory')
    # Never write an audit into its own frozen source/transfer roles.
    need(not args.output.resolve().is_relative_to(args.release_dir) and args.output.resolve() not in {args.freeze.resolve(), args.source_review.resolve()}, 'Output cannot overwrite source/provenance')
    audit = Audit()
    audit.file(__file__)
    channels, reconstruction, frozen_replay = [], {}, dict(performed=False, supplemental_only=True)
    execution = None
    candidate_status = None
    try:
        cfg, freeze = source_closure(audit, args)
        record, scheduler = transfer_closure(audit, args, freeze)
        execution, _ = preliminary_execution(audit, args, cfg, freeze, scheduler)
        target = audit.read(args.job_dir/'results/target.toml')
        reference = audit.read(args.job_dir/'results/reference.json')
        saved = audit.read(args.job_dir/'results/comparison.json')
        need(saved['schema'] == 'X13_CHARGED_SPECTRAL_COMPARISON_V1' and saved['catalogue_id'] == 'X13', 'Saved candidate comparison identity')
        candidate_status = saved['scientific_status']
        for raw, process in zip((target, reference), execution['processes'][:2]):
            need(process['start_unix'] <= raw['completion']['finished_unix'] <= process['end_unix'], 'Actual raw completion inside its observed process')
        need(execution['processes'][0]['start_unix'] <= target['metadata']['started_unix'] <= target['completion']['finished_unix'], 'Native adapter start/finish in observed process')
        remote_result = PurePosixPath(record['result_root'])
        remote_source = PurePosixPath(record['source_root'])
        expected_comparison_inputs = {str(remote_result/'target.toml'): sha(args.job_dir/'results/target.toml'),
                                      str(remote_result/'reference.json'): sha(args.job_dir/'results/reference.json'),
                                      str(remote_source/'manifest.toml'): sha(args.release_dir/'manifest.toml')}
        need(saved['input_sha256'] == expected_comparison_inputs, 'Saved comparator exact original source/raw SHA map')
        t = audit_backend(audit, target, 'native_ed_new_adapter', args.release_dir, cfg, args.freeze_sha256)
        r = audit_backend(audit, reference, 'independent_quspin', args.release_dir, cfg, args.freeze_sha256)
        for mode, (tc, rc) in enumerate(zip(target['site_annihilators'], reference['site_annihilators'])):
            audit.compare('cross site annihilator'+str(mode), maxabs(sparse(tc)-sparse(rc)), GATES['matrix_abs'])
        channels = compare_backends(audit, t, r)
        science = status(audit)
        audit.check('saved candidate outcome agrees with independent full raw audit', saved['scientific_status'] == science)
        need(len(channels) == 128 and all(len(d['blocks']) == 13 for d in t.values()) and all(len(d['blocks']) == 13 for d in r.values()), 'Exact complete128channels/104spinblocks across backends')
        reconstruction = dict(complete_cases=4, spinblocks_per_backend=52, channels=128,
                              initial_ensemble='complete canonical N4 uniform ground projector',
                              all_invisible_poles_retained=True, hidden_state_witness_required=False,
                              hidden_lowest_channels=sum(x['hidden_lowest'] for x in channels),
                              ground_dimensions={cid: dict(target=t[cid]['ground_dimension'], reference=r[cid]['ground_dimension']) for cid, *_ in CASE_GRID},
                              channel_receipts=channels)
        if args.replay:
            need(audit.integrity(), 'No replay after failed primary integrity')
            frozen_replay = replay_comparator(audit, args, science)
    except Exception as e:
        audit.check('complete independently auditable receipt', False, type(e).__name__+': '+str(e))
    for path, expected in list(audit.inputs.items()):
        audit.check('bound input unchanged '+path, Path(path).is_file() and sha(path) == expected)
    final_science = status(audit)
    result = dict(schema=SCHEMA, catalogue_id='X13', job_id=args.job_id,
                  created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                  auditor_sha256=sha(__file__), integrity_status='PASS' if audit.integrity() else 'UNRESOLVED',
                  audit_status='PASS' if audit.integrity() else 'UNRESOLVED', scientific_status=final_science,
                  candidate_status=candidate_status,
                  scheduler_clock=dict(naive_sacct_timezone='America/New_York EDT UTC-04:00 on2026-10-22/23', rounding_tolerance_seconds=1),
                  scientific_sampling=False, scientific_completions=0, ledger_writes=False,
                  approved_freeze_sha256=args.freeze_sha256,
                  approved_source_review_sha256=args.source_review_sha256,
                  input_sha256=audit.inputs, check_count=len(audit.checks),
                  failure_count=sum(not r['passed'] for r in audit.checks),
                  qualification_failure_count=sum(not r['passed'] for r in audit.numerical),
                  cross_failure_count=sum(not r['passed'] for r in audit.comparisons),
                  failures=[r for r in audit.checks if not r['passed']],
                  checks=audit.checks, qualification_checks=audit.numerical,
                  cross_checks=audit.comparisons, independent_reconstruction=reconstruction,
                  frozen_replay=frozen_replay,
                  scope='Finite L4 exact charged-sector/native Hamiltonian adapter only; no native dynamical solver, bulk extrapolation, or ledger promotion.')
    raw = (json.dumps(result, allow_nan=False, ensure_ascii=True, indent=2)+'\n').encode('utf-8')
    with args.output.open('xb') as f:
        f.write(raw)
    print(json.dumps({k: result[k] for k in ('audit_status', 'scientific_status', 'failure_count', 'check_count')}))
    return 0 if result['audit_status'] == 'PASS' else 2


if __name__ == '__main__':
    raise SystemExit(main())


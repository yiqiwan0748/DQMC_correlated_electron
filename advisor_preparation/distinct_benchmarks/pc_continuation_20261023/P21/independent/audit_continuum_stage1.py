#!/usr/bin/env python3
"""Independent, local-only X01 continuum stage1 receipt audit.

No sampler, native Julia entry point, or comparator function is imported. The
frozen comparator is replayed as a separate, analysis-only subprocess. Raw-bin
jackknives, rational fit weights, Student critical value, intervals, mixing
metrics and a bit-Fock exact reference are reconstructed independently here.
An UNRESOLVED scientific result can have a successful integrity audit.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
from fractions import Fraction as F
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
import tempfile
import tomllib

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np

SCHEMA = 'X01_CONTINUUM_STAGE1_INDEPENDENT_AUDIT_V1'
APPROVED_RELEASE = 'continuum_production_v1c'
APPROVED_FREEZE_SHA = '659b4456ba007aed1618e2f937fd79f87133ff89fbf62e07bfb5aacc69b35c7f'
APPROVED_ARCHIVE_SHA = '597f263d66efbecf283421a72cf3dbc7806b96b996904bedba2e09fe3e868a0b'
APPROVED_MANIFEST_SHA = 'bfd9f808b9b25bf4e93393d3a05b2e958cf903ff4a43660744585100f1ac09ea'
APPROVED_TARGET_SHA = 'ff967af02086422e30cb69986f76626bef6881dcd0435ae156d1cd17f6d4a591'
H = [0.25, 0.125, 0.0625, 0.03125]
EXPECTED_GATES = dict(energy_margin=.02, S2_margin=.01,
    endpoint_maximum_halfwidth_fraction=.5, fit_difference_margin_fraction=.5,
    cumulative_deltaG_max=1e-8, oracle_matrix_abs=1e-10,
    oracle_scalar_abs=1e-10, oracle_positive_eigenvalue_min=-1e-12,
    positive_phase_abs=1e-10, blocking_error_relative_change_max=.5,
    split_Rhat_max=1.05, effective_bins_min_per_chain=16.,
    effective_bins_min_total=128., acceptance_min_exclusive=0.,
    acceptance_max_exclusive=1.)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    p = Path(path)
    return tomllib.loads(p.read_text(encoding='utf-8')) if p.suffix == '.toml' else json.loads(p.read_text(encoding='utf-8-sig'))


def finite(value):
    if isinstance(value, dict):
        return all(finite(v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return all(finite(v) for v in value)
    return not isinstance(value, float) or math.isfinite(value)


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def differences(a, b, path='', out=None):
    """Path-sensitive recursive comparison; gates/labels remain exact."""
    out = [] if out is None else out
    if len(out) >= 50:
        return out
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            out.append(path + ': key grid differs')
        for key in set(a) & set(b):
            differences(a[key], b[key], path + '/' + str(key), out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(path + ': length differs')
        for i, (x, y) in enumerate(zip(a, b)):
            differences(x, y, path + '/' + str(i), out)
    elif isinstance(a, bool) or isinstance(b, bool):
        if a is not b:
            out.append(path + ': Boolean differs')
    elif isinstance(a, (float, int)) and isinstance(b, (float, int)):
        if not math.isclose(a, b, rel_tol=2e-10, abs_tol=2e-12):
            out.append(path + ': numeric mismatch ' + repr(a) + ' vs ' + repr(b))
    elif a != b:
        out.append(path + ': value differs')
    return out


class Audit:
    def __init__(self):
        self.checks = []
        self.inputs = {}

    def check(self, name, ok, detail=None):
        self.checks.append(dict(name=name, passed=bool(ok), detail=plain(detail)))
        return bool(ok)

    def file(self, path):
        p = Path(path)
        if not self.check('required file ' + str(p), p.is_file()):
            raise ValueError('Missing required receipt: ' + str(p))
        self.inputs[str(p.resolve())] = sha(p)
        return p

    def load(self, path):
        return read(self.file(path))

    def equal(self, name, a, b):
        d = differences(plain(a), plain(b))
        return self.check(name, not d, d or None)


def paired_jackknife(numerators, denominators, entropy=False):
    """Explicit delete-one-bin resampling with paired simultaneous node deletion."""
    n, d = np.asarray(numerators, float), np.asarray(denominators, float)
    if n.shape != d.shape or n.ndim != 2 or not np.isfinite(n).all() or not np.isfinite(d).all():
        raise ValueError('Invalid paired raw-bin grid')
    B = n.shape[1]
    if B < 2:
        raise ValueError('At least two bins required')
    def statistic(mask):
        ns = [math.fsum(row[mask]) for row in n]
        ds = [math.fsum(row[mask]) for row in d]
        if min(ds) <= 0:
            raise ValueError('Nonpositive jackknife phase denominator')
        ratios = [x/y for x, y in zip(ns, ds)]
        if entropy:
            if min(ratios) <= 0:
                raise ValueError('Nonpositive entropy ratio')
            return -math.fsum(math.log(x) for x in ratios)
        if len(ratios) != 1:
            raise ValueError('Scalar statistic requires one stream')
        return ratios[0]
    allbins = np.ones(B, dtype=bool)
    full = statistic(allbins)
    leave = []
    for j in range(B):
        mask = allbins.copy()
        mask[j] = False
        leave.append(statistic(mask))
    avg = math.fsum(leave)/B
    return dict(plug_in=full, bias_corrected=B*full-(B-1)*avg,
        standard_error=math.sqrt((B-1)/B*math.fsum((x-avg)**2 for x in leave)),
        bins=B, leave_one_out_min=min(leave), leave_one_out_max=max(leave))


def merge_bins(a):
    return (np.asarray(a, float)[:, ::2]+np.asarray(a, float)[:, 1::2])/2


def fraction_inverse(matrix):
    n = len(matrix)
    a = [[F(v) for v in row]+[F(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if a[i][k])
        a[k], a[pivot] = a[pivot], a[k]
        q = a[k][k]
        a[k] = [v/q for v in a[k]]
        for i in range(n):
            if i != k:
                q = a[i][k]
                a[i] = [x-q*y for x, y in zip(a[i], a[k])]
    return [row[n:] for row in a]


def fit_design():
    """Exact rational normal equations, independently of comparator SVD/pinv."""
    # h=1/4,1/8,... therefore h²=1/16,1/64,...
    x = [F(1, 16*4**i) for i in range(4)]
    def weights(indices, degree):
        inv = fraction_inverse([[sum(x[i]**(j+k) for i in indices) for k in range(degree+1)] for j in range(degree+1)])
        result = [[F(0)]*4 for _ in range(degree+1)]
        for j in range(degree+1):
            for i in indices:
                result[j][i] = sum(inv[j][k]*x[i]**k for k in range(degree+1))
        return result
    full, fine, coarse, quartic = weights(range(4), 1), weights(range(1, 4), 1), weights(range(3), 1), weights(range(4), 2)
    residuals = [[F(i == j)-full[0][j]-x[i]*full[1][j] for j in range(4)] for i in range(4)]
    result = dict(primary=full[0], fine=fine[0], quartic=quartic[0],
        primary_minus_fine=[a-b for a, b in zip(full[0], fine[0])],
        primary_minus_quartic=[a-b for a, b in zip(full[0], quartic[0])],
        residuals=residuals,
        heldout_coarse=[F(j == 0)-fine[0][j]-x[0]*fine[1][j] for j in range(4)],
        heldout_fine=[F(j == 3)-coarse[0][j]-x[3]*coarse[1][j] for j in range(4)])
    return {k: np.asarray(v, float) for k, v in result.items()}


def student7_cdf(t):
    """Integrate cos^6 with a recurrence, rather than comparator trig expansion."""
    z = math.atan(t/math.sqrt(7))
    integral = z
    for n in (2, 4, 6):
        integral = math.sin(z)*math.cos(z)**(n-1)/n+(n-1)/n*integral
    return .5+16/(5*math.pi)*integral


def student7_critical():
    target = 1-.01/(3*26)/2
    lo, hi = 0., 100.
    for _ in range(100):
        mid = (lo+hi)/2
        if student7_cdf(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def interval(name, kind, values, margin, critical):
    a = np.asarray(values, float)
    if a.shape != (8,) or not np.isfinite(a).all():
        raise ValueError('Exactly eight finite bundle contrasts required')
    mean = math.fsum(a)/8
    se = math.sqrt(math.fsum((v-mean)**2 for v in a)/(7*8))
    half = critical*se
    lower, upper = mean-half, mean+half
    return dict(name=name, kind=kind, margin=margin,
        equivalent=lower >= -margin and upper <= margin,
        wholly_outside=upper < -margin or lower > margin,
        mean=mean, standard_error=se, halfwidth=half, lower=lower, upper=upper, bundles=8, df=7)


def rank_normalize(a):
    shape = a.shape
    flat = a.ravel()
    order = np.argsort(flat, kind='stable')
    z = np.empty(flat.size)
    start = 0
    while start < len(order):
        end = start+1
        while end < len(order) and flat[order[end]] == flat[order[start]]:
            end += 1
        rank = (start+1+end)/2
        z[order[start:end]] = statistics.NormalDist().inv_cdf((rank-.375)/(len(flat)+.25))
        start = end
    return z.reshape(shape)


def rhat(a):
    W = np.var(a, axis=1, ddof=1).mean()
    B = a.shape[1]*np.var(a.mean(axis=1), ddof=1)
    return (1. if B == 0 else None) if W == 0 else max(1., math.sqrt(((a.shape[1]-1)*W+B)/(a.shape[1]*W)))


def ess(a):
    """Independent NumPy autocovariances and initial-positive/monotone pairs."""
    a = np.asarray(a, float)
    scale = np.max(abs(a))
    if scale:
        a = a/scale
    m, n = a.shape
    means = a.mean(axis=1)
    W = np.var(a, axis=1, ddof=1).mean()
    if W <= 0:
        return None
    centered = a-means[:, None]
    ac = np.array([np.sum(centered[:, :n-lag]*centered[:, lag:])/(m*n) for lag in range(n)])
    if m == 1:
        rho = ac/ac[0]
    else:
        B = n*np.var(means, ddof=1)
        vp = (n-1)/n*W+B/n
        rho = np.r_[1., 1-(W-ac[1:])/vp]
    pairs = []
    for j in range(0, n-1, 2):
        pair = rho[j]+rho[j+1]
        if not math.isfinite(pair) or pair <= 0:
            break
        pairs.append(min(pair, pairs[-1]) if pairs else pair)
    return min(float(m*n), m*n/max(1., -1+2*math.fsum(pairs)))


def own_mixing(chains):
    a = np.asarray(chains, float)
    if a.shape != (8, 64) or not np.isfinite(a).all():
        raise ValueError('Mixing grid must be8x64 finite')
    split = a.reshape(8, 2, 32).reshape(16, 32)
    ranked = rank_normalize(split)
    scaled = split/max(np.max(abs(split)), 1e-300)
    ordered = np.sort(scaled.ravel())
    low, high = ordered[255:257]
    folded_values = np.where(scaled <= low, low-scaled, scaled-high)+(high-low)/2
    folded = rank_normalize(folded_values)
    constant = np.ptp(a, axis=1) == 0
    rhats = [rhat(ranked), rhat(folded)] if np.ptp(a) else [None, None]
    total = [ess(split), ess(ranked)]
    per = []
    for row in a:
        pair = [ess(row[None, :]), ess(rank_normalize(row[None, :]))]
        per.append(None if None in pair else min(pair))
    R = None if None in rhats else max(rhats)
    E = None if None in total else min(total)
    good = not constant.any() and R is not None and R <= 1.05 and E is not None and E >= 128 and all(v is not None and v >= 16 for v in per)
    return dict(eligible=bool(good), rhat=R, rank_split_rhat=rhats[0], folded_rank_split_rhat=rhats[1], ess_total=E, ess_per_chain=per)


def exact_reference():
    """Independent full256-state bit-Fock matrices and fermionic partial trace.

    Mode order is up0..3,down0..3. The cut reorders modes to upA,downA,
    upB,downB with sign(-1)^(NdownA*NupB); this is not a bosonic trace.
    No native ED/QuSpin/comparator code is imported.
    """
    K = np.zeros((256, 256))
    double, number = np.zeros(256), np.zeros(256)
    permutation, signs = [], []
    for state in range(256):
        u, d = state & 15, state >> 4
        double[state], number[state] = (u & d).bit_count(), state.bit_count()
        a, b = (u & 3)+4*(d & 3), (u >> 2)+4*(d >> 2)
        permutation.append(16*a+b)
        signs.append((-1)**((d & 3).bit_count()*(u >> 2).bit_count()))
        for spin in range(2):
            for i in range(4):
                for j in ((i+1) % 4, (i-1) % 4):
                    p, q = i+4*spin, j+4*spin
                    if state >> q & 1 and not state >> p & 1:
                        z = state ^ (1 << q)
                        parity = (state & ((1 << q)-1)).bit_count()+(z & ((1 << p)-1)).bit_count()
                        K[z ^ (1 << p), state] -= (-1)**parity
    physical = K+np.diag(4*double)
    V = 4*(double-.5*number+1)
    perm = np.argsort(permutation)
    signs = np.asarray(signs)[perm]
    def describe(rho):
        reordered = rho[np.ix_(perm, perm)]*signs[:, None]*signs[None, :]
        cut = np.einsum('abcb->ac', reordered.reshape(16, 16, 16, 16))
        purity = float(np.trace(cut@cut))
        return dict(physical_energy=float(np.einsum('ij,ji', physical, rho)), Nmean=float(np.dot(number, np.diag(rho))),
            purity=purity, S2=-math.log(purity), trace=float(np.trace(rho)), cut_trace=float(np.trace(cut)))
    w, Q = np.linalg.eigh(K+np.diag(V))
    p = np.exp(-(w-w[0])); p /= sum(p)
    continuum = describe((Q*p)@Q.T)
    wK, QK = np.linalg.eigh(K)
    rows = []
    for h in H:
        d = np.exp(-h*V/2)
        kinetic = (QK*np.exp(-h*wK))@QK.T
        S = d[:, None]*kinetic*d[None, :]
        w, Q = np.linalg.eigh(S)
        if min(w) <= 0:
            raise ValueError('Independent transfer spectrum is nonpositive')
        p = np.exp(round(1/h)*(np.log(w)-math.log(max(w)))); p /= sum(p)
        sym = (Q*p)@Q.T
        asym = d[:, None]*sym/d[None, :]
        rows.append(dict(dtau=h, slices=round(1/h), asymmetric=describe(asym), symmetric=describe(sym)))
    return dict(hilbert_dimension=256, continuum=continuum, finite_dtau=rows)


def derived_seed(root, k, replica):
    mask = (1 << 64)-1
    x = root ^ (((k+1)*0x9e3779b97f4a7c15) & mask) ^ ((replica*0xbf58476d1ce4e5b9) & mask)
    z = (x+0x9e3779b97f4a7c15) & mask
    z = ((z ^ (z >> 30))*0xbf58476d1ce4e5b9) & mask
    z = ((z ^ (z >> 27))*0x94d049bb133111eb) & mask
    return (z ^ (z >> 31)) % ((1 << 63)-1)


def expected_seed_rows(namespace):
    rows = []
    for i, h in enumerate(H, 1):
        for b in range(1, 9):
            def label(stream):
                s = f'{namespace}|dtau_index={i}|bundle={b}|stream={stream}'
                return int.from_bytes(hashlib.sha256(s.encode()).digest()[:8], 'big') & ((1 << 63)-1)
            e, root = label('energy'), label('entropy')
            rows.append(dict(dtau_index=i, dtau=h, bundle=b, energy_seed=e, entropy_root_seed=root,
                node_replica_seeds=[[derived_seed(root, k, r) for r in (1, 2)] for k in range(4)]))
    return rows


def duration_seconds(value):
    value = value.strip()
    days = 0
    if '-' in value:
        ds, value = value.split('-', 1); days = int(ds)
    parts = [float(x) for x in value.split(':')]
    return days*86400+sum(v*60**i for i, v in enumerate(reversed(parts)))


def resource_receipt(audit, path, timeout):
    text = audit.file(path).read_text(encoding='utf-8')
    def item(label):
        match = re.search(r'^\s*'+re.escape(label)+r':\s*(.+)$', text, re.M)
        if not match:
            raise ValueError('Missing time resource field ' + label)
        return match.group(1)
    wall = duration_seconds(item('Elapsed (wall clock) time (h:mm:ss or m:ss)'))
    rss = int(item('Maximum resident set size (kbytes)'))
    exit_code = int(item('Exit status'))
    audit.check(path.name+' resource caps', 0 < wall <= timeout+1 and 0 < rss <= 3500*1024, dict(wall_seconds=wall, max_rss_kb=rss, exit_code=exit_code))
    return dict(wall_seconds=wall, max_rss_kb=rss, exit_code=exit_code)


def audit_receipts(a, audit):
    source, native, result = a.release_dir.resolve(), a.native_source_dir.resolve(), a.result_dir.resolve()
    cfg = audit.load(source/'manifest.toml')
    audit.check('approved original scientific manifest digest', sha(source/'manifest.toml') == APPROVED_MANIFEST_SHA)
    checks = dict(schema='X01_CONTINUUM_PRODUCTION_V1', phase='X01_CONTINUUM_STAGE1', stage_index=1,
        catalogue_id='X01', reserved_maximum_stages=3, automatic_escalation=False,
        scientific_completion_allowed=True, production_campaign_implemented=True,
        L=4, t=1., U=4., beta=1., mu_centered=0., mu_physical=2., boundary='pbc', nA=2,
        production_dtau=H, reference_dtau=H, bundles_per_dtau=8, burn_sweeps=2048,
        measured_sweeps=16384, bins=64, n_stab=2, sweeps_per_measurement=1,
        lambda_increments=4, lambda_indices=[0, 1, 2, 3], lambda_values=[0., .25, .5, .75],
        delta_lambda=.25, cut_sites_one_based=[1, 2],
        seed_namespace='X01_CONTINUUM_PRODUCTION_V1_STAGE1_20261023_16384')
    audit.equal('original fixed stage1 model and sampling parameters', {k: cfg.get(k) for k in checks}, checks)
    audit.equal('original scientific and numerical thresholds', cfg['gates'], EXPECTED_GATES)
    st = cfg['statistics']
    audit.check('original26-CI3-stage7df allocation', st['independent_units'] == 8 and st['degrees_of_freedom'] == 7 and st['family_alpha'] == .01 and st['reserved_stages'] == 3 and st['intervals_per_stage'] == 26 and math.isclose(st['per_interval_two_sided_alpha'], .01/(3*26), rel_tol=0, abs_tol=1e-18))
    critical = student7_critical()
    audit.check('independent Student7 critical reconstruction', abs(st['student_t_critical']-critical) < 2e-9, dict(critical=critical, declared=st['student_t_critical'], cdf=student7_cdf(st['student_t_critical'])))
    # Bind every frozen source file, without accepting traversal paths.
    sha_lines = audit.file(source/'SHA256SUMS').read_text().splitlines()
    frozen = {}
    for line in sha_lines:
        match = re.fullmatch(r'([0-9a-f]{64})\s+\*?(.+)', line.strip())
        if not match:
            raise ValueError('Malformed SHA256SUMS')
        want, name = match.groups()
        path = (source/name).resolve()
        if not path.is_relative_to(source) or name in frozen:
            raise ValueError('Unsafe or duplicate frozen path')
        frozen[name] = want
        audit.check('release frozen hash '+name, sha(audit.file(path)) == want)
    for name in ('manifest.toml', 'seeds.toml', 'excluded_pilot_seeds.toml', 'target_production.jl', 'oracle_quspin.py', 'compare_production.py', 'mixing_diagnostics.py', 'production.sbatch', 'generate_seeds.py'):
        audit.check('required source covered by freeze '+name, name in frozen)
    audit.check('corrected symmetric Blom source denominator', '(len(values) + 1.0 / 4.0)' in (source/'mixing_diagnostics.py').read_text())
    driver = (source/'target_production.jl').read_bytes()
    audit.check('approved mechanical retry ASCII division and exact driver', driver.isascii() and driver.count(b'bs=div(CFG["measured_sweeps"],B)') == 1 and sha(source/'target_production.jl') == APPROVED_TARGET_SHA)
    if a.freeze:
        freeze = audit.load(a.freeze)
        audit.check('approved v1c pre-execution freeze', sha(a.freeze) == APPROVED_FREEZE_SHA and freeze.get('schema') == 'X01_CONTINUUM_V1C_SOURCE_FREEZE' and freeze.get('archive_sha256') == APPROVED_ARCHIVE_SHA and len(freeze.get('files', {})) == 23)
        audit.check('mechanical retry excludes all cancelled v1b partial results', freeze.get('mechanical_retry_of_job') == '6460728' and freeze.get('same_predeclared_seed_grid_retry') is True and freeze.get('all_prior_partial_outputs_excluded') is True and freeze.get('additional_independent_chains') == 0)
        if a.freeze_sha256:
            audit.check('externally supplied freeze digest', sha(a.freeze) == a.freeze_sha256.lower())
        files = freeze.get('files', {})
        if isinstance(files, list):
            files = {v['path']: v['sha256'] for v in files}
        matched = 0
        for name, want in files.items():
            if isinstance(want, dict):
                want = want.get('sha256')
            options = [source/name, source.parent/name, Path(a.freeze).parent/name]
            p = next((p for p in options if p.is_file()), None)
            audit.check('outer freeze hash '+name, p is not None and sha(p) == want)
            matched += p is not None
        audit.check('outer freeze contains release binding', matched > 0)
    else:
        audit.check('external pre-execution freeze receipt supplied', False, 'Supply --freeze and preferably --freeze-sha256; SHA256SUMS alone is internally consistent, not an independently anchored freeze.')
    for name, want in cfg['source_hashes'].items():
        audit.check('native source pin '+name, sha(audit.file(native/name)) == want)
    manifest = audit.load(native/'Manifest.toml')
    dep = manifest['deps']['JDQMCFramework'][0]
    audit.check('native JDQMCFramework version/tree pin', dep['version'] == cfg['jdqmcframework_version'] and dep['git-tree-sha1'] == cfg['jdqmcframework_tree_sha1'])
    seeds, pilot = audit.load(source/'seeds.toml'), audit.load(source/'excluded_pilot_seeds.toml')
    expected = expected_seed_rows(cfg['seed_namespace'])
    audit.equal('independent integer-only288seed regeneration', seeds['bundles'], expected)
    actual = lambda rows: [v for r in rows for v in [r['energy_seed']]+[s for pair in r['node_replica_seeds'] for s in pair]]
    labels = lambda rows: [r[k] for r in rows for k in ('energy_seed', 'entropy_root_seed')]
    actual_new, old, streams, oldstreams = actual(expected), actual(pilot['bundles']), labels(expected), labels(pilot['bundles'])
    audit.check('288actual positive unique;64labels disjoint from all pilot initializations', len(actual_new) == len(set(actual_new)) == 288 and min(actual_new) > 0 and len(set(streams)) == 64 and not (set(actual_new+streams) & set(old+oldstreams)))
    raw, oracle, saved = audit.load(result/'target.toml'), audit.load(result/'oracle.json'), audit.load(result/'comparison.json')
    audit.check('all retained raw/reference fields finite', finite(raw) and finite(oracle))
    meta, om = raw['metadata'], oracle['metadata']
    audit.check('job identity across raw/oracle/comparator', meta['job_id'] == om['job_id'] == saved['job_id'] == a.job_id)
    audit.check('driver/source/manifest/seed raw closure', meta['manifest_sha256'] == sha(source/'manifest.toml') and meta['adapter_sha256'] == sha(source/'target_production.jl') and meta['seeds_sha256'] == sha(source/'seeds.toml') and meta['source_hashes'] == cfg['source_hashes'])
    audit.check('target pinned Julia/JDQMC and single threads', meta['julia_version'] == cfg['julia_version'] and meta['jdqmcframework_version'] == cfg['jdqmcframework_version'] and meta['julia_threads'] == meta['blas_threads'] == 1)
    keys = ('production_dtau', 'bundles_per_dtau', 'burn_sweeps', 'measured_sweeps', 'bins', 'n_stab', 'sweeps_per_measurement', 'lambda_indices', 'lambda_values', 'delta_lambda', 'cut_sites_one_based')
    audit.equal('raw complete sampling parameter binding', meta['sampling'], {k: cfg[k] for k in keys})
    audit.check('target complete with zero automatic completions', raw['completion']['status'] == 'PRODUCTION_STAGE1_COMPLETE' and raw['completion']['errors'] == [] and raw['completion']['scientific_completions'] == meta['scientific_completions'] == saved['scientific_completions'] == 0 and meta['phase'] == cfg['phase'] and meta['schema'] == cfg['schema'])
    audit.check('oracle independent source/runtime/model closure', om['manifest_sha256'] == sha(source/'manifest.toml') and om['adapter_sha256'] == sha(source/'oracle_quspin.py') and om['phase'] == 'REFERENCE_ONLY' and om['scientific_completions'] == 0 and om['versions'] == {k: cfg[k+'_version'] for k in ('quspin', 'numpy', 'scipy')} and oracle['hilbert_dimension'] == 256)
    for prefix, diagnostic in [('continuum', oracle['diagnostics'])]+[(str(r['dtau']), r['diagnostics']) for r in oracle['finite_dtau']]:
        for name, value in diagnostic.items():
            audit.check('oracle numerical diagnostic '+prefix+'/'+name, math.isfinite(value) and (value >= -1e-12 if name.endswith('_min') else abs(value) <= 1e-10), value)
    independent = exact_reference()
    audit.equal('independent full-Fock continuum reference', oracle['continuum'], independent['continuum'])
    stripped = [{k: r[k] for k in ('dtau', 'slices', 'asymmetric', 'symmetric')} for r in oracle['finite_dtau']]
    audit.equal('independent four-step transfer/reference and fermionic cut', stripped, independent['finite_dtau'])
    exits = dict(line.split('=', 1) for line in audit.file(result/'process_exit_codes.txt').read_text().splitlines() if line.strip())
    audit.check('oracle/target successful process exits', exits.get('oracle_exit') == exits.get('target_exit') == '0')
    audit.check('comparator exit corresponds to declared candidate result', exits.get('comparator_exit') == ('0' if saved['scientific_status'] == 'PASS' else '2'))
    resources = {name: resource_receipt(audit, result/(name+'.resources.txt'), cfg['resources'][name+'_timeout_seconds']) for name in ('oracle', 'target', 'comparator')}
    for name, driver in [('oracle', 'oracle_quspin.py'), ('target', 'target_production.jl'), ('comparator', 'compare_production.py')]:
        command = (result/(name+'.resources.txt')).read_text()
        audit.check(name+' timed command binds frozen driver and timeout', driver in command and str(cfg['resources'][name+'_timeout_seconds']) in command and 'timeout' in command)
    for name in ('oracle', 'target'):
        audit.check(name+' resource exit zero', resources[name]['exit_code'] == 0)
    audit.check('comparator resource exit matches process', resources['comparator']['exit_code'] == int(exits['comparator_exit']))
    audit.check('positive bounded target driver elapsed', 0 < raw['completion']['elapsed_seconds'] <= 3300)
    audit.check('positive bounded oracle adapter elapsed', 0 < om['elapsed_seconds'] <= 90)
    scheduler = a.scheduler or result.parent/'scheduler.txt'
    rows = list(csv.DictReader(io.StringIO(audit.file(scheduler).read_text(encoding='utf-8-sig')), delimiter='|'))
    job = [r for r in rows if r.get('JobIDRaw') == a.job_id]
    batch = [r for r in rows if r.get('JobIDRaw') == a.job_id+'.batch']
    audit.check('unique scheduler main and batch identities', len(job) == len(batch) == 1)
    if job and batch:
        j, b = job[0], batch[0]
        expected_exit = '0:0' if saved['scientific_status'] == 'PASS' else '2:0'
        audit.check('scheduler terminal exit reflects candidate status', j['ExitCode'] == b['ExitCode'] == expected_exit and j['State'] == b['State'] == ('COMPLETED' if saved['scientific_status'] == 'PASS' else 'FAILED'))
        wall = duration_seconds(j['Elapsed'])
        audit.check('scheduler oneCPU3500M60min admission', j['ReqCPUS'] == '1' and j['ReqMem'] == '3500M' and 0 < wall <= 3601)
        audit.check('resource wall fits scheduler wall', sum(r['wall_seconds'] for r in resources.values()) <= wall+3)
    transfer_path = a.transfer or result.parent/'transfer.json'
    transfer = audit.load(transfer_path)
    audit.check('transfer binds same job and read-only retrieval', str(transfer['job_id']) == a.job_id and transfer.get('remote_writes') is False and transfer.get('binaries_copied') is False)
    names = [r['local'] for r in transfer['files']]
    audit.check('transfer paths unique', len(names) == len(set(names)))
    for row in transfer['files']:
        p = (Path(transfer_path).parent/row['local']).resolve()
        if not p.is_relative_to(Path(transfer_path).parent.resolve()):
            raise ValueError('Unsafe transfer receipt path')
        audit.check('transfer byte/hash '+row['local'], row.get('status') == 'RETRIEVED' and p.is_file() and p.stat().st_size == row['bytes'] and sha(p) == row['sha256'])
    if 'sacct' in transfer:
        audit.check('scheduler exact transfer binding', audit.file(scheduler).read_text(encoding='utf-8-sig') == transfer['sacct'])
    return cfg, raw, oracle, saved, expected, resources


def reconstruct(a, audit, cfg, raw, oracle, saved, seeds):
    result = a.result_dir.resolve()
    er, sr = raw['energy_rows'], raw['entropy_node_rows']
    energy = {(r['dtau_index'], r['bundle']): r for r in er}
    nodes = {(r['dtau_index'], r['bundle'], r['lambda_index']): r for r in sr}
    expected = {(i, b) for i in range(1, 5) for b in range(1, 9)}
    audit.check('complete unique32energy128lambda raw grid', len(er) == len(energy) == 32 and set(energy) == expected and len(sr) == len(nodes) == 128 and set(nodes) == {(i, b, k) for i, b in expected for k in range(4)})
    Y = {ob: np.zeros((8, 4)) for ob in ('energy', 'S2')}
    bundles, groups, blocks, runtime_sum = [], {}, [], 0.
    for seed in seeds:
        i, b = seed['dtau_index'], seed['bundle']
        E, S = energy[i, b], [nodes[i, b, k] for k in range(4)]
        for k, r in [(None, E)]+list(enumerate(S)):
            name = f'dt{i}_bundle{b}_'+('energy' if k is None else 'lambda'+str(k))
            audit.equal(name+' side/aggregate identity', audit.load(result/(name+'.toml')), r)
            audit.check(name+' completed grid/seed parameters', r['status'] == 'COMPLETED' and r['dtau'] == H[i-1] and r['dtau_index'] == i and r['bundle'] == b and (r['seed'] == seed['energy_seed'] if k is None else r['root_seed'] == seed['entropy_root_seed'] and [r['seed1'], r['seed2']] == seed['node_replica_seeds'][k] and r['lambda_index'] == k and r['lambda_value'] == k/4 and r['delta_lambda'] == .25 and r['cut_sites_one_based'] == [1, 2]))
            arrays = ['numerator_bins', 'denominator_bins']+(['N_numerator_bins'] if k is None else ['native_avg_den_bins', 'native_avg_num_bins'])
            for key in arrays:
                arr = np.asarray(r[key], float)
                audit.check(name+'64finite bins '+key, arr.shape == (64,) and np.isfinite(arr).all())
            audit.check(name+' positive phase bins', max(abs(v-1) for v in r['denominator_bins']) <= 1e-10)
            delta = r['cumulative_deltaG' if k is None else 'native_cumulative_deltaG']
            audit.check(name+' stability and acceptance', 0 <= delta <= 1e-8 and 0 < r['acceptance'] < 1)
            audit.check(name+' finite positive bounded runtime', math.isfinite(r['elapsed_seconds']) and 0 < r['elapsed_seconds'] <= 3300)
            runtime_sum += r['elapsed_seconds']
            if k is None:
                history = [r['deltaG_after_burnin']]+r['cumulative_deltaG_bins']
                audit.check(name+'65point cumulative history includes burn', len(history) == 65 and min(history) >= 0 and all(x <= y for x, y in zip(history, history[1:])) and history[-1] == delta)
        en, ed = [E['numerator_bins']], [E['denominator_bins']]
        sn, sd = [r['numerator_bins'] for r in S], [r['denominator_bins'] for r in S]
        estimates, merged = {}, {}
        for ob, n, d, entropy in [('energy', en, ed, False), ('S2', sn, sd, True)]+[('node'+str(k), [sn[k]], [sd[k]], True) for k in range(4)]:
            estimates[ob] = paired_jackknife(n, d, entropy)
            merged[ob] = paired_jackknife(merge_bins(n), merge_bins(d), entropy)
            v, w = estimates[ob]['standard_error'], merged[ob]['standard_error']
            eligible = v > 0 and w > 0 and abs(w/v-1) <= .5
            blocks.append(dict(dtau_index=i, bundle=b, observable=ob, eligible=eligible, se64=v, se32=w))
        for ob in Y:
            Y[ob][b-1, i-1] = estimates[ob]['bias_corrected']
        bundles.append(dict(dtau_index=i, bundle=b, dtau=H[i-1], estimates=estimates, merged32bin=merged))
        groups[i, b] = (np.asarray(en[0]), np.asarray(ed[0]), np.asarray(sn), np.asarray(sd))
    audit.check('sum individual API runtimes bounded by complete driver runtime', 0 < runtime_sum <= raw['completion']['elapsed_seconds']+1, runtime_sum)
    audit.equal('independent paired jackknives and64to32blocking estimates', saved['bundles'], bundles)
    mixing = []
    for i in range(1, 5):
        rows = [groups[i, b] for b in range(1, 9)]
        sn, sd = np.asarray([r[2] for r in rows]), np.asarray([r[3] for r in rows])
        ratios = np.sum(sn, axis=(0, 2))/np.sum(sd, axis=(0, 2))
        means = np.mean(sd, axis=(0, 2))
        influence = -(sn-ratios[None, :, None]*sd)/(ratios[None, :, None]*means[None, :, None])
        streams = dict(energy=np.asarray([r[0]/r[1] for r in rows]), S2=np.sum(influence, axis=1))
        streams.update({'node'+str(k): influence[:, k, :] for k in range(4)})
        for ob, chains in streams.items():
            own = own_mixing(chains)
            original = next(r['diagnostics'] for r in saved['mixing'] if r['dtau_index'] == i and r['observable'] == ob)
            audit.equal(f'independent mixing metrics dt{i}/{ob}', {k: original[k] for k in own}, own)
            mixing.append(dict(dtau_index=i, observable=ob, diagnostics=own))
    D, critical = fit_design(), cfg['statistics']['student_t_critical']
    audit.equal('exact rational fixed OLS contrast vectors', saved['fit_design'], D)
    intervals, exactfits, precision, bias, heldouts = [], {}, [], [], []
    for ob, values in Y.items():
        margin, key = (.02, 'physical_energy') if ob == 'energy' else (.01, 'S2')
        y0 = oracle['continuum'][key]
        ref = np.asarray([r['asymmetric'][key] for r in oracle['finite_dtau']])
        exactfits[ob] = {k: float(D[k]@ref-y0) for k in ('primary', 'fine', 'quartic')}
        exactfits[ob].update({k: float(D[k]@ref) for k in ('heldout_coarse', 'heldout_fine')})
        for fit in ('primary', 'fine', 'quartic'):
            r = interval(ob+'_'+fit+'_continuum', 'endpoint', values@D[fit]-y0, margin, critical)
            intervals.append(r)
            precision.append(dict(name=r['name'], eligible=r['halfwidth'] <= margin/2, halfwidth=r['halfwidth'], maximum=margin/2))
        for contrast in ('primary_minus_fine', 'primary_minus_quartic'):
            intervals.append(interval(ob+'_'+contrast, 'fit_stability', values@D[contrast], margin/2, critical))
        for j in range(4):
            intervals.append(interval(f'{ob}_primary_residual_h{j+1}', 'fit_residual', values@D['residuals'][j], margin, critical))
        own = []
        for j in range(4):
            r = interval(f'{ob}_finite_h{j+1}', 'finite_h', values[:, j]-ref[j], margin, critical)
            intervals.append(r); own.append(r)
        for contrast in ('heldout_coarse', 'heldout_fine'):
            vals = values@D[contrast]
            heldouts.append(dict(observable=ob, contrast=contrast, mean=math.fsum(vals)/8,
                standard_error=float(np.std(vals, ddof=1)/math.sqrt(8)), bundles=8,
                interpretation='Descriptive only; no acceptance interval, scientific gate or additional inferential claim.'))
        lower = own[0]['lower']+ref[0]-own[3]['upper']-ref[3]
        upper = own[0]['upper']+ref[0]-own[3]['lower']-ref[3]
        bias.append(dict(observable=ob, lower=float(lower), upper=float(upper), exact_coarse_minus_finest=float(ref[0]-ref[3]),
            method='Minkowski difference of existing simultaneous intervals; no new statistical look', eligible=bool(upper < 0 and ref[0] < ref[3])))
    audit.equal('independent26simultaneous confidence intervals', saved['intervals'], intervals)
    audit.equal('independent deterministic fit admissibility', saved['exact_fit_diagnostics'], exactfits)
    audit.equal('independent required endpoint halfwidths', saved['endpoint_precision'], precision)
    audit.equal('independent derived coarse-minus-finest bias intervals', saved['derived_bias_detection'], bias)
    audit.equal('heldouts descriptive only and excluded from26gate intervals', saved['heldout_descriptive'], heldouts)
    counts = {k: sum(r['kind'] == k for r in intervals) for k in ('endpoint', 'fit_stability', 'fit_residual', 'finite_h')}
    audit.equal('26CI family exact allocation', counts, dict(endpoint=6, fit_stability=4, fit_residual=8, finite_h=8))
    fits_ok = all(r['equivalent'] for r in intervals if r['kind'] in ('fit_stability', 'fit_residual'))
    endpoints = [r for r in intervals if r['kind'] == 'endpoint']
    finiteh = [r for r in intervals if r['kind'] == 'finite_h']
    # Do not accept the saved comparator's eligibility flag as evidence. Every
    # source/grid/runtime/reference invariant above has its own independent
    # check; fixed exact-reference fit admissibility is reconstructed here.
    numerical_raw_ok = all(r['passed'] for r in audit.checks) and all(abs(exactfits[ob][name]) <= (.02 if ob == 'energy' else .01)/4 for ob in Y for name in ('primary', 'fine', 'quartic'))
    independent_mix_ok = all(r['diagnostics']['eligible'] for r in mixing) and all(r['eligible'] for r in blocks)
    precision_ok = all(r['eligible'] for r in precision)
    bias_ok = all(r['eligible'] for r in bias)
    components = dict(raw_source_numerical_mixing_eligible=bool(numerical_raw_ok and independent_mix_ok),
        precision_adequate=bool(precision_ok), fit_adequate=bool(fits_ok),
        finite_h_equivalent=all(r['equivalent'] for r in finiteh),
        continuum_equivalent=all(r['equivalent'] for r in endpoints),
        finite_step_bias_demonstrated=bool(bias_ok))
    audit.equal('independent scientific decision components', saved['decision_components'], components)
    if not numerical_raw_ok or not independent_mix_ok or not precision_ok:
        status = 'UNRESOLVED'
    elif any(r['wholly_outside'] for r in finiteh):
        status = 'FAIL'
    elif not fits_ok or not bias_ok:
        status = 'UNRESOLVED'
    elif any(r['wholly_outside'] for r in endpoints):
        status = 'FAIL'
    elif all(r['equivalent'] for r in finiteh+endpoints):
        status = 'PASS'
    else:
        status = 'UNRESOLVED'
    audit.check('independent scientific decision matches frozen candidate', status == saved['scientific_status'], dict(independent=status, saved=saved['scientific_status']))
    return dict(scientific_status=status, interval_counts=counts, fit_design=D, bundle_estimates=Y,
        intervals=intervals, endpoint_precision=precision, mixing=mixing, blocking=blocks,
        derived_bias_detection=bias, heldout_descriptive=heldouts, decision_components=components,
        eligibility_disclosure='Source/runtime/grid/reference integrity, fixed-reference fit admissibility and mixing eligibility are reconstructed independently; complete comparator output is additionally replayed. These checks do not prove mixing or tail coverage.')


def replay(a, audit, saved):
    # Keep temporary analysis output beside the explicitly requested audit
    # output. The managed Windows runtime's default temp may not be writable.
    with tempfile.TemporaryDirectory(prefix='x01_continuum_replay_', dir=a.output.parent.resolve()) as temp:
        output = Path(temp)/'comparison.json'
        cmd = [sys.executable, str(a.release_dir.resolve()/'compare_production.py'), '--source', str(a.release_dir.resolve()),
            '--native-root', str(a.native_source_dir.resolve()), '--result-dir', str(a.result_dir.resolve()),
            '--target-exit-code', '0', '--oracle-exit-code', '0', '--output', str(output)]
        env = os.environ.copy(); env['PYTHONDONTWRITEBYTECODE'] = '1'
        run = subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=env)
        audit.check('frozen analysis-only comparator replay expected exit', run.returncode == (0 if saved['scientific_status'] == 'PASS' else 2), dict(exit_code=run.returncode, stdout=run.stdout[-3000:], stderr=run.stderr[-3000:]))
        audit.check('replay produced a retained comparable result', output.is_file())
        replayed = read(output)
        # Local and remote absolute input keys differ after transfer; byte hashes
        # are compared as multisets and every copied file is separately bound.
        audit.equal('replayed comparator input hash multiset', sorted(saved['input_sha256'].values()), sorted(replayed['input_sha256'].values()))
        left = {k: v for k, v in saved.items() if k != 'input_sha256'}
        right = {k: v for k, v in replayed.items() if k != 'input_sha256'}
        audit.equal('complete frozen comparator output replay', left, right)
        return dict(exit_code=run.returncode, replay_sha256=sha(output), saved_sha256=sha(a.result_dir/'comparison.json'),
            numeric_replay_tolerance=dict(absolute=2e-12, relative=2e-10), gates_and_labels_exact=True)


def self_test():
    # Mathematical tests only: deterministic arrays, no scientific chain run.
    audit = Audit()
    n = [[float(j*j+3*j+7) for j in range(64)]]
    d = [[1.+(j % 3)/10 for j in range(64)]]
    q = paired_jackknife(n, d)
    assert math.isfinite(q['standard_error']) and q['standard_error'] > 0
    designs = fit_design()
    x = np.square(H)
    for name in ('primary', 'fine', 'quartic'):
        assert abs(sum(designs[name])-1) < 1e-14 and abs(designs[name]@x) < 1e-16
    assert abs(designs['quartic']@np.square(x)) < 1e-16
    assert abs(student7_critical()-7.58185668916078) < 2e-9
    refs = exact_reference()
    assert abs(refs['continuum']['physical_energy']+.6218669203630302) < 1e-12
    assert abs(refs['continuum']['S2']-1.9307833675726933) < 1e-12
    # Deliberately constant scientific chains must never appear eligible.
    assert not own_mixing(np.ones((8, 64)))['eligible']
    # Correct Blom transform is symmetric: denominator S+1/4. The earlier
    # source proposal's S-1/4 typo is never accepted by this independent audit.
    ranks = rank_normalize(np.arange(64., dtype=float).reshape(1, 64))
    assert np.max(abs(ranks+ranks[:, ::-1])) < 1e-14
    print('PASS: deterministic mathematical auditor self-test; zero scientific sampling')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--result-dir', type=Path)
    p.add_argument('--release-dir', '--source', type=Path)
    p.add_argument('--native-source-dir', '--native-root', type=Path)
    p.add_argument('--job-id')
    p.add_argument('--output', type=Path)
    p.add_argument('--scheduler', type=Path)
    p.add_argument('--transfer', type=Path)
    p.add_argument('--freeze', type=Path)
    p.add_argument('--freeze-sha256')
    p.add_argument('--self-test', action='store_true')
    a = p.parse_args()
    if a.self_test:
        self_test(); return
    for name in ('result_dir', 'release_dir', 'native_source_dir', 'job_id', 'output'):
        if getattr(a, name) is None:
            p.error('--'+name.replace('_', '-')+' is required')
    if a.output.exists():
        raise ValueError('Preserve historical audit outputs; choose a new --output')
    a.output.parent.mkdir(parents=True, exist_ok=True)
    audit = Audit()
    out = dict(schema=SCHEMA, catalogue_id='X01', scientific_completions=0,
        approved_release=APPROVED_RELEASE, approved_source_freeze_sha256=APPROVED_FREEZE_SHA,
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), auditor_sha256=sha(__file__),
        job_id=a.job_id, remote_calls=False, scientific_sampling=False,
        limitations=['Eight independent bundle vectors support the fixed26interval family;64bins do not create64independent bundles.',
            'Nominal Bonferroni99percent coverage assumes approximately normal independent bundle contrasts. Extreme Student7tails/nonlinearS2coverage are not calibrated.',
            'Independent rank/foldedRhat,ESS and blocking calculations are rejection diagnostics, not proofs of equilibrium or tail coverage.',
            'The frozen native kernels are byte-bound; this audit does not independently reimplement the sampler.',
            'Independent bit-Fock references use a fermionic cut and the original physical energy; no finite-step correction is subtracted from fitted sample endpoints.'])
    try:
        cfg, raw, oracle, saved, seeds, resources = audit_receipts(a, audit)
        out['candidate_scientific_status'] = saved['scientific_status']
        out['independent_reconstruction'] = reconstruct(a, audit, cfg, raw, oracle, saved, seeds)
        out['frozen_comparator_replay'] = replay(a, audit, saved)
        out['resources'] = resources
    except Exception as error:
        audit.check('audit completed without missing/malformed receipts', False, type(error).__name__+': '+str(error))
    failures = [r for r in audit.checks if not r['passed']]
    out.update(audit_status='PASS' if not failures else 'UNRESOLVED',
        scientific_status=out.get('independent_reconstruction', {}).get('scientific_status', 'UNRESOLVED') if not failures else 'UNRESOLVED',
        catalogue_status='PENDING_ROOT_REVIEW' if not failures else 'PENDING_INDEPENDENT_AUDIT',
        check_count=len(audit.checks), failure_count=len(failures), failures=failures,
        checks=audit.checks, input_sha256=audit.inputs)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open('x', encoding='utf-8') as f:
        json.dump(plain(out), f, indent=2, sort_keys=True, allow_nan=False); f.write('\n')
    print(f"audit={out['audit_status']}; science={out['scientific_status']}; checks={len(audit.checks)}; failures={len(failures)}; zero automatic completions")
    if failures:
        sys.exit(2)


if __name__ == '__main__':
    main()

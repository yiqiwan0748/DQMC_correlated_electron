#!/usr/bin/env python3
"""Independent full-grid X16 receipt auditor candidate: NumPy and stdlib only.

No native/reference solver, sampling, network, admission, ledger or canonical
writes. Primary math imports no candidate module. Optional comparator replay is
supplemental and unavailable before complete original provenance/raw closure.
"""
from __future__ import annotations
import argparse,csv,datetime as dt,hashlib,io,json,math,os,re,sys,tarfile,tempfile,tomllib
from pathlib import Path,PurePosixPath
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
SCHEMA='X16_INDEPENDENT_RAW_AUDIT_V1'
RAW_SCHEMA='X16_PAIR_BINDING_RAW_V1'
RELEASE='x16_pair_binding_v1'
CATALOGUE_SHA='59a780c5994db50a56da260fc51ae639f22b674b3bc24f878cde1a6e7bc31818'
INHERITED_HELPER_SOURCE_SHA='b3b69b4e4fbcc1b009ef0b24ec9eb45303264ac6f792607a9c8772f2212bd80c'
NATIVE_PINS={'Project.toml':'ea3b7eacfcf9a36e4652654a5e7b573ef27529edfbad2e272da2bafca1a7f58f','Manifest.toml':'5cefaf1acdb4fa23a456c2c58f20fb1fcebd54b3cd27bc448cb65b0fcfecb90e','analysis/compute_s2_ed.jl':'f8f94cfc13e11bb8b5e64e8ac8f4ade5ed3a13b3976c83b5ab5db89d86ba5f28'}
GATES=dict(matrix_abs=1e-10,eigen_residual_abs=1e-10,orthogonality_abs=1e-10,cluster_diameter_max=1e-8,cluster_separation_min=1e-6,cross_energy_abs=1e-9,cross_projector_abs=1e-9,cross_pair_density_abs=1e-9,cross_binding_abs=4e-9,car_abs=1e-10,pair_operator_abs=1e-10,pair_density_hermitian_abs=1e-10,pair_density_eigenvalue_min=-1e-10,pair_trace_abs=1e-10,particle_number_abs=1e-10,onsite_identity_abs=1e-10,pair_diagonal_identity_abs=1e-10,maximum_block_energy_budget=1e-8,maximum_binding_halfwidth=5e-8,require_all_binding_signs_resolved=True)
CASE_GRID=[('U2_pbc',2.,'pbc',1.),('U2_apbc',2.,'apbc',1.),('U4_pbc',4.,'pbc',2.),('U4_apbc',4.,'apbc',2.)]
SOURCE_ROLES=('manifest.toml','PROTOCOL.md','PROVENANCE.json','raw_io.py','target.jl','reference.py','compare.py','run_candidate.py','pilot.sbatch','production.sbatch','README.md','test_candidate.py','catalogue_X16.json')
CUTOFF=1792755300
PRIOR_JOBS=('6460990','6460991','6460992','6460994','6460996','6461010','6462017','6462170','6462171','6462172','6462173','6462174')
REMOTE_ROOT='/lustre/isaac24/scratch/ywang287/advisor-validation-20261022'
ACTIVE={'RUNNING','COMPLETING','CONFIGURING','SUSPENDED','RESIZING'}
TERMINAL={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED','SPECIAL_EXIT'}
EPS64=float(64*np.finfo(float).eps)
OUTPUT_ROLES={'preflight.json','active_queue.txt','previous_twelve_terminal.txt','scontrol.txt','target.toml','reference.json','comparison.json','execution.json','process_exit_codes.txt'}|{n+'.'+s for n in ('target','reference','comparator')for s in ('stdout.txt','stderr.txt','resources.txt','process.json')}
ORIGINAL_AUDITOR_SHA='5c3332567aa8ddc8600418d587d1c6f483deb02133874ec2c0595a60e1c999ba'
ORIGINAL_WORKFLOW_SHA='a064d86fe2d3c41235ef48fb74d29b1268c7dbc50c4397a74a0b3340022410d4'
ORIGINAL_PILOT='6462342'
ORIGINAL_FAILED_AUDIT_SHA='85878a438e7d095f6cc665ab4cebc98106c0f8133aba84f90fde8006a5f623ad'
STAGE3_CLOSURE_SHA='aee287fbb68366bd703a6d83debec7bd9f1fe887371058f8fba7f60ea27c6044'
CANCELLED_BEFORE_START={'6462172':'2026-10-23T03:12:42','6462173':'2026-10-23T03:12:45'}

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


def eigensystem_errors(H, E, V):
    H, E, V = np.asarray(H), np.asarray(E), np.asarray(V)
    need(H.ndim == 2 and H.shape[0] == H.shape[1] and V.shape == H.shape and E.shape == (len(H),), 'Complete eigensystem dimensions')
    need(np.all(np.isfinite(H)) and np.all(np.isfinite(E)) and np.all(np.isfinite(V)) and np.all(np.diff(E) >= 0), 'Finite ascending eigensystem')
    return dict(Hermiticity=maxabs(H-H.conj().T),
                orthogonality=maxabs(V.conj().T@V-np.eye(len(E))),
                residual=maxabs(H@V-V*E))


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


def sparse(obj, dimension=256, columns=None):
    """Complete rectangular source/operator triplets, without tolerance pruning."""
    columns = dimension if columns is None else columns
    need(type(obj['rows']) is int and type(obj['columns']) is int and obj['rows'] == dimension and obj['columns'] == columns and 0 <= dimension <= 256 and 0 <= columns <= 256, 'Exact rectangular sparse shape')
    fields = [obj[k] for k in ('row0', 'column0', 'real', 'imag')]
    need(all(isinstance(v, list) for v in fields) and len({len(v) for v in fields}) == 1, 'Complete sparse triplet fields')
    z = np.zeros((dimension, columns), complex); seen = set()
    for i,j,x,y in zip(*fields):
        need(type(i) is int and type(j) is int and 0 <= i < dimension and 0 <= j < columns and (i,j) not in seen, 'Unique in-range triplet coordinates')
        need(type(x) in (int,float) and type(y) in (int,float) and math.isfinite(x) and math.isfinite(y), 'Finite exact numeric triplet values')
        z[i,j] = complex(x,y); seen.add((i,j))
    return z


def sparse_product(A,B):
    need(A.ndim == B.ndim == 2 and A.shape[1] == B.shape[0], 'Operator composition shapes')
    out=np.zeros((A.shape[0],B.shape[1]),complex)
    for k in range(B.shape[0]):
        cols=np.flatnonzero(B[k])
        if len(cols):out[:,cols] += A[:,k,None]*B[k,cols]
    return out


def car_error(operators):
    n=operators[0].shape[0];ident=np.eye(n);error=0.
    for i,a in enumerate(operators):
        for j,b in enumerate(operators):
            error=max(error,maxabs(sparse_product(a,b)+sparse_product(b,a)),maxabs(sparse_product(a,b.conj().T)+sparse_product(b.conj().T,a)-(ident if i==j else 0)))
    return error


def block_grid(L=4):
    return [(N,u,N-u) for N in (2,3,4,5,6) for u in range(L+1) if 0 <= N-u <= L]


def pair_grid(modes=8):return [(a,b) for a in range(modes) for b in range(a+1,modes)]


def typed_members(value):
    need(isinstance(value,list) and all(isinstance(r,dict) and set(r)=={'block_key','eigen_index0'} and isinstance(r['block_key'],str) and type(r['eigen_index0']) is int and r['eigen_index0']>=0 for r in value),'Exact typed complete cluster member records')
    return value


def scalar_formula(audit,name,saved,computed,scale=None):
    need(type(saved) in (int,float) and type(computed) in (int,float) and math.isfinite(saved) and math.isfinite(computed), 'Finite scalar formula inputs')
    scale=max(abs(saved),abs(computed),1e-300) if scale is None else scale
    need(type(scale) in (int,float) and math.isfinite(scale) and scale>0, 'Positive declared formula conditioning scale')
    return audit.check(name,abs(saved-computed)<=EPS64*scale,dict(saved=saved,recomputed=computed,tolerance=EPS64*scale,conditioning_scale=scale))


def energy_budget_from_norms(n,gram,residual,H_fro,D_max):
    """Polar similarity/Weyl enclosure, with frozen explicit FP allowances."""
    need(type(n) is int and n>0 and all(type(x) in (int,float) and math.isfinite(x) and x>=0 for x in (gram,residual,H_fro,D_max)), 'Complete nonnegative finite eigensystem norms')
    delta=gram+EPS64*n;pad=EPS64*n*max(1.,H_fro,D_max)
    need(delta<1, 'Full square eigenvector matrix must be invertible')
    root=math.sqrt(1-delta)
    eta=(residual+pad)/root+2*D_max*delta/(root*(1+root))+pad
    return dict(dimension=n,gram_fro=gram,residual_fro=residual,H_fro=H_fro,D_max=D_max,delta=delta,roundoff=pad,energy_budget=eta)


def energy_budget(H,E,V):
    n=len(E);need(H.shape==V.shape==(n,n) and np.asarray(E).shape==(n,), 'Whole square eigensystem budget')
    gram=float(np.linalg.norm(V.conj().T@V-np.eye(n),'fro'))
    residual=float(np.linalg.norm(H@V-V*E,'fro'))
    return energy_budget_from_norms(n,gram,residual,float(np.linalg.norm(H,'fro')),maxabs(E))


def diagnostic_budget_range(own):
    """Propagate two implementations' declared primitive-norm roundoff only.

    Both implementations are checked against the same saved primitive norms.
    This interval qualifies saved comparator diagnostics, never scientific
    energy/sign/precision gates, which use our independently computed budget.
    """
    n=own['dimension'];dh=2*EPS64*max(1.,own['H_fro']);dd=2*EPS64*max(1.,own['D_max']);dg=2*EPS64*n
    upper_pad=EPS64*n*max(1.,own['H_fro']+dh,own['D_max']+dd);dr=2*upper_pad
    low=energy_budget_from_norms(n,max(0.,own['gram_fro']-dg),max(0.,own['residual_fro']-dr),max(0.,own['H_fro']-dh),max(0.,own['D_max']-dd))['energy_budget']
    high=energy_budget_from_norms(n,own['gram_fro']+dg,own['residual_fro']+dr,own['H_fro']+dh,own['D_max']+dd)['energy_budget']
    return low,high


def budget_receipt(audit,label,saved,own):
    keys=('gram_fro','residual_fro','H_fro','D_max','delta','roundoff','energy_budget')
    numeric_list([saved[k] for k in keys]);need(type(saved['dimension']) is int and saved['dimension']==own['dimension'], 'Budget dimension matches complete block')
    audit.check(label+' saved nonnegative norms/budget',all(saved[k]>=0 for k in keys))
    n=own['dimension']
    audit.check(label+' independent primitive gram norm',abs(saved['gram_fro']-own['gram_fro'])<=EPS64*n)
    audit.check(label+' independent primitive residual norm',abs(saved['residual_fro']-own['residual_fro'])<=own['roundoff'])
    for k in ('H_fro','D_max'):audit.check(label+' independent primitive '+k,abs(saved[k]-own[k])<=EPS64*max(1.,own[k]))
    formula=energy_budget_from_norms(n,saved['gram_fro'],saved['residual_fro'],saved['H_fro'],saved['D_max'])
    for k in ('delta','roundoff','energy_budget'):scalar_formula(audit,label+' own saved formula '+k,saved[k],formula[k])
    audit.qualified(label+' own complete spectral energy budget',own['energy_budget'],GATES['maximum_block_energy_budget'])
    audit.qualified(label+' saved complete spectral energy budget',saved['energy_budget'],GATES['maximum_block_energy_budget'])


def physical_interval(mu,N,E,lower,upper):
    pad=EPS64*max(1.,abs(E),abs(mu*N))
    return dict(energy=E+mu*N,lower=lower+mu*N-pad,upper=upper+mu*N+pad,physical_shift_roundoff=pad)


def binding_interval(name,odd,low,high):
    scale=max(1.,2*abs(odd['energy_physical'])+abs(low['energy_physical'])+abs(high['energy_physical']))
    pad=EPS64*scale
    value=2*odd['energy_physical']-low['energy_physical']-high['energy_physical']
    lower=2*odd['energy_lower']-low['energy_upper']-high['energy_upper']-pad
    upper=2*odd['energy_upper']-low['energy_lower']-high['energy_lower']+pad
    sign='BOUND' if lower>0 else 'NONBINDING' if upper<0 else 'UNRESOLVED'
    return dict(name=name,odd_N=odd['N'],even_N=[low['N'],high['N']],value=value,lower=lower,upper=upper,halfwidth=(upper-lower)/2,roundoff=pad,classification=sign,conditional_numerical_budget=True)


def binding_receipt(audit,label,saved,own,saved_ground):
    odd=own['odd_N'];low,high=own['even_N'];formula=binding_interval(own['name'],saved_ground[odd],saved_ground[low],saved_ground[high])
    need(type(saved['odd_N']) is int and saved['odd_N']==odd and saved['even_N']==[low,high] and all(type(x) is int for x in saved['even_N']) and saved['conditional_numerical_budget'] is True, 'Exact fixed binding descriptor and conditional scope')
    scale=max(1.,2*abs(saved_ground[odd]['energy_physical'])+abs(saved_ground[low]['energy_physical'])+abs(saved_ground[high]['energy_physical']))
    for k in ('value','lower','upper'):scalar_formula(audit,label+' exact saved propagation '+k,saved[k],formula[k],scale)
    scalar_formula(audit,label+' saved halfwidth own endpoints',saved['halfwidth'],(saved['upper']-saved['lower'])/2)
    scalar_formula(audit,label+' saved roundoff own ground inputs',saved['roundoff'],formula['roundoff'])
    audit.check(label+' saved ordered enclosure',saved['lower']<=saved['value']<=saved['upper'])
    saved_sign='BOUND' if saved['lower']>0 else 'NONBINDING' if saved['upper']<0 else 'UNRESOLVED'
    audit.check(label+' saved sign exact own interval',saved['classification']==saved_sign)
    audit.qualified(label+' independent and saved sign agree',0. if saved['classification']==own['classification'] else 1.,0.)
    audit.qualified(label+' own binding halfwidth',own['halfwidth'],GATES['maximum_binding_halfwidth'])


def audit_quspin_bridge(audit,raw,canonical_C):
    b=raw['metadata']['quspin_bridge'];raw_C=[sparse(v) for v in b['raw_site_annihilators']]
    need(len(raw_C)==8,'All eight QuSpin raw fermion operators')
    states=b['raw_basis_states'];need(len(states)==256 and all(type(x) is int for x in states) and set(states)==set(range(256)),'Complete unreduced QuSpin integer basis')
    reverse=lambda x:sum(((x>>j)&1)<<(3-j) for j in range(4))
    coords=coordinate_pairs(b['raw_coordinates']);audit.check('QuSpin integer bit coordinates independently derived',coords==[(reverse(s>>4),reverse(s&15)) for s in states])
    for i,c in enumerate(raw_C):
        expected=np.array([float((u if i<4 else d)>>(i%4)&1) for u,d in coords])
        audit.check('QuSpin raw number coordinate '+str(i),close(sparse_product(c.conj().T,c),np.diag(expected),GATES['matrix_abs']))
    audit.qualified('QuSpin raw full CAR',car_error(raw_C),GATES['car_abs'])
    T=creation_bridge(raw_C,coords);mapping=np.argmax(abs(T),axis=0).tolist();phase=T[mapping,np.arange(256)]
    need(all(type(i) is int for i in b['canonical_to_raw']),'Integer bridge map without coercion')
    audit.check('QuSpin analytic creation-string coordinate bridge',b['canonical_to_raw']==mapping)
    audit.check('QuSpin analytic creation-string phase bridge',close(phase,np.array(numeric_list(b['phases_real']))+1j*np.array(numeric_list(b['phases_imag'])),1e-12))
    audit.check('QuSpin nonfitted creation-string derivation',b['derivation']=='vacuum plus canonical creation strings, descending application order; not fitted to target')
    for i,(c,want) in enumerate(zip(raw_C,canonical_C)):
        audit.qualified('QuSpin EXACT canonical C bridge '+str(i),0. if np.array_equal(T.conj().T@c@T,want) else 1.,0.)
    return T


def audit_case(audit,c,case,operators,pairs,label,bridge=None):
    prefix=label+'/'+case['case_id'];coords=coordinate_grid();nominal=canonical_hamiltonian(4,1.,case['U'],case['mu'],case['boundary'])
    need(all(c[k]==case[k] for k in ('case_id','U','mu','boundary')) and all(type(c[k]) in (int,float) for k in ('U','mu')),'Fixed typed original case metadata')
    grid=block_grid();need([r['block_key'] for r in c['blocks']]==[block_name(*x) for x in grid] and len(c['blocks'])==19,'Complete ordered19-spinblock grid')
    vectors={};values={};budgets={};saved_budgets={};budget_ranges={};records={N:[] for N in range(2,7)}
    if bridge is not None:audit.qualified(prefix+' EXACT original raw QuSpin H bridge',0. if np.array_equal(bridge.conj().T@sparse(c['raw_quspin_hamiltonian'])@bridge,nominal) else 1.,0.)
    for row,(N,u,d) in zip(c['blocks'],grid):
        key=block_name(N,u,d);name=prefix+'/'+key;coord=[p for p in coords if p[0].bit_count()==u and p[1].bit_count()==d];ids=[a*16+b for a,b in coord];n=len(ids)
        need(all(type(row[k]) is int for k in ('N','Nup','Ndown')) and (row['N'],row['Nup'],row['Ndown'])==(N,u,d) and coordinate_pairs(row['coordinates'])==coord,'Full exact block coordinate basis')
        H=matrix(row['Hamiltonian'],n,n);V=matrix(row['eigenvectors'],n,n);E=np.array(numeric_list(row['eigenvalues']));need(E.shape==(n,),'Every ordered block eigenvalue')
        ownH=nominal[np.ix_(ids,ids)];audit.qualified(name+' EXACT nominal integer H',0. if np.array_equal(H,ownH) else 1.,0.)
        errors=eigensystem_errors(ownH,E,V)
        for k,gate in [('Hermiticity','matrix_abs'),('residual','eigen_residual_abs'),('orthogonality','orthogonality_abs')]:audit.qualified(name+' '+k,errors[k],GATES[gate])
        budget=energy_budget(ownH,E,V);budget_receipt(audit,name,row['numerical_energy_budget'],budget);budgets[key]=budget['energy_budget'];saved_budgets[key]=row['numerical_energy_budget']['energy_budget'];budget_ranges[key]=diagnostic_budget_range(budget)
        for j,e in enumerate(E):
            full=np.zeros(256,complex);full[ids]=V[:,j];vectors[key,j]=full;values[key,j]=float(e);records[N].append(dict(energy=float(e),block_key=key,eigen_index0=j))
    need([r['N'] for r in c['sector_spectra']]==list(range(2,7)) and all(type(r['N']) is int for r in c['sector_spectra']),'Complete ordered five sectors')
    spectra={N:clusters(records[N]) for N in range(2,7)}
    for row in c['sector_spectra']:
        N=row['N'];groups=spectra[N];need(len(row['clusters'])==len(groups) and sum(len(g['members']) for g in groups)==math.comb(8,N),'Full unfiltered fixed-N spectrum cardinality')
        audit.check(prefix+f'/N{N} complete exact cluster memberships',[typed_members(r['members']) for r in row['clusters']]==[r['members'] for r in groups]);bands=[]
        for i,(saved,own) in enumerate(zip(row['clusters'],groups)):
            audit.check(prefix+f'/N{N}/cluster{i} own saved energy extrema',close(numeric_list([saved[k] for k in ('energy','energy_min','energy_max')]),[own[k] for k in ('energy','energy_min','energy_max')],GATES['matrix_abs']))
            lo=min(values[m['block_key'],m['eigen_index0']]-budgets[m['block_key']] for m in own['members']);hi=max(values[m['block_key'],m['eigen_index0']]+budgets[m['block_key']] for m in own['members']);bands.append((lo,hi));audit.qualified(prefix+f'/N{N}/cluster{i} uncertainty-qualified full diameter',hi-lo,GATES['cluster_diameter_max'])
        audit.qualified(prefix+f'/N{N} uncertainty-qualified cluster separation',max([GATES['cluster_separation_min']-(b[0]-a[1]) for a,b in zip(bands,bands[1:])]+[0.]),0.)
    need([r['N'] for r in c['grounds']]==list(range(2,7)) and all(type(r['N']) is int for r in c['grounds']),'Five full global ground ensembles')
    derived=[];saved_ground={}
    for row in c['grounds']:
        N=row['N'];name=prefix+f'/N{N}';members=spectra[N][0]['members'];G=np.column_stack([vectors[m['block_key'],m['eigen_index0']] for m in members]);dg=G.shape[1]
        need(type(row['dimension']) is int and row['dimension']==dg and typed_members(row['members'])==members,'Every global spinblock ground member weighted uniformly')
        keys=[block_name(*x) for x in grid if x[0]==N];minima={key:min(r['energy'] for r in records[N] if r['block_key']==key) for key in keys};E=min(minima.values())
        lo=min(minima[k]-budgets[k] for k in keys);hi=min(minima[k]+budgets[k] for k in keys);phys=physical_interval(case['mu'],N,E,lo,hi)
        source_lo=min(minima[k]-saved_budgets[k] for k in keys);source_hi=min(minima[k]+saved_budgets[k] for k in keys);sourcephys=physical_interval(case['mu'],N,E,source_lo,source_hi)
        formula=dict(energy_mu=E,energy_mu_lower=source_lo,energy_mu_upper=source_hi,energy_physical=sourcephys['energy'],energy_lower=sourcephys['lower'],energy_upper=sourcephys['upper'],physical_shift_roundoff=sourcephys['physical_shift_roundoff'])
        for k,want in formula.items():scalar_formula(audit,name+' saved own '+k,row[k],want,max(1.,abs(E),abs(case['mu']*N)) if k in ('energy_physical','energy_lower','energy_upper') else None)
        audit.check(name+' own saved individual ground energies',close(numeric_list(row['individual_energies']),[values[m['block_key'],m['eigen_index0']] for m in members],GATES['matrix_abs']))
        SG=matrix(row['vectors'],256,dg);P=G@G.conj().T
        audit.check(name+' saved global basis orthonormal',close(SG.conj().T@SG,np.eye(dg),GATES['orthogonality_abs']))
        audit.check(name+' saved same complete global projector',close(SG@SG.conj().T,P,GATES['matrix_abs']))
        audit.check(name+' saved explicit full projector',close(sparse(row['projector']),P,GATES['matrix_abs']))
        sources=[p@SG for p in pairs];need(len(row['pair_sources'])==28,'All28 pair source matrices')
        for i,(saved,S,(a,b)) in enumerate(zip(row['pair_sources'],sources,pair_grid())):
            need(all(type(saved[k]) is int for k in ('pair_index','a','b','final_N')) and (saved['pair_index'],saved['a'],saved['b'],saved['final_N'])==(i,a,b,N-2),'Exact pair-source descriptor/finalN')
            audit.check(name+f'/pair{i} source from own ground basis',close(sparse(saved['vectors'],256,dg),S,GATES['pair_operator_abs']))
        gamma=np.array([[np.vdot(a,b)/dg for b in sources] for a in sources]);saved_gamma=matrix(row['pair_density'],28,28)
        audit.check(name+' saved whole784complex pairdensity',close(saved_gamma,gamma,GATES['matrix_abs']))
        audit.qualified(name+' pairdensity Hermiticity',maxabs(saved_gamma-saved_gamma.conj().T),GATES['pair_density_hermitian_abs'])
        audit.qualified(name+' pairdensity PSD',-float(np.linalg.eigvalsh((saved_gamma+saved_gamma.conj().T)/2).min()),-GATES['pair_density_eigenvalue_min'])
        audit.qualified(name+' pairtrace chooseN2',abs(np.trace(saved_gamma)-N*(N-1)/2),GATES['pair_trace_abs'])
        probability=np.sum(abs(SG)**2,axis=1)/dg;bits=np.array([[float((u if i<4 else d)>>(i%4)&1) for u,d in coords] for i in range(8)])
        occupations=bits@probability;diagonal=np.array([np.dot(bits[a]*bits[b],probability) for a,b in pair_grid()]);onsite=np.array([np.dot(bits[j]*bits[j+4],probability) for j in range(4)])
        audit.check(name+' saved all8 occupations',close(numeric_list(row['occupations']),occupations,GATES['matrix_abs']))
        audit.qualified(name+' total fixed particleN',abs(sum(occupations)-N),GATES['particle_number_abs'])
        audit.check(name+' saved all28 direct numberpair expectations',close(numeric_list(row['number_pair_diagonal']),diagonal,GATES['matrix_abs']))
        audit.qualified(name+' all28 pairdiagonal number identities',maxabs(np.diag(saved_gamma)-diagonal),GATES['pair_diagonal_identity_abs'])
        audit.check(name+' saved all4 onsite double occupations',close(numeric_list(row['onsite_double']),onsite,GATES['onsite_identity_abs']))
        audit.qualified(name+' onsite pair identities',maxabs(np.array([saved_gamma[pair_grid().index((j,j+4)),pair_grid().index((j,j+4))] for j in range(4)])-onsite),GATES['onsite_identity_abs'])
        shift=case['mu']*N;pad=phys['physical_shift_roundoff'];arithmetic=EPS64*max(1.,abs(E),abs(shift))
        lower_range=[min(minima[k]-budget_ranges[k][1] for k in keys)+shift-pad-arithmetic,min(minima[k]-budget_ranges[k][0] for k in keys)+shift-pad+arithmetic]
        upper_range=[min(minima[k]+budget_ranges[k][0] for k in keys)+shift+pad-arithmetic,min(minima[k]+budget_ranges[k][1] for k in keys)+shift+pad+arithmetic]
        saved_ground[N]=row;derived.append(dict(N=N,dimension=dg,energy_mu=E,energy_physical=phys['energy'],energy_lower=phys['lower'],energy_upper=phys['upper'],projector=P,pair_density=gamma,diagnostic_lower_range=lower_range,diagnostic_upper_range=upper_range))
    gm={r['N']:r for r in derived};bindings=[binding_interval('B_h',gm[3],gm[2],gm[4]),binding_interval('B_e',gm[5],gm[4],gm[6])]
    need([r['name'] for r in c['bindings']]==['B_h','B_e'],'Both fixed binding combinations')
    for saved,own in zip(c['bindings'],bindings):binding_receipt(audit,prefix+'/'+own['name'],saved,own,saved_ground)
    return dict(case_id=case['case_id'],nominal_H=nominal,all_eigenvalues={N:sorted(r['energy'] for r in records[N]) for N in range(2,7)},grounds=derived,bindings=bindings)


def protocol(cfg):
    need(cfg['schema']=='X16_PAIR_BINDING_PROTOCOL_V1' and cfg['catalogue_id']=='X16' and cfg['release']==RELEASE,'Fixed X16 protocol identity')
    need(cfg['L']==4 and cfg['t']==1. and cfg['total_N']==[2,3,4,5,6] and cfg['gates']==GATES,'Fixed model/sectors/all numerical and physical gates')
    need([(c['case_id'],c['U'],c['boundary'],c['mu']) for c in cfg['cases']]==CASE_GRID,'Exact four finite cases')
    need(cfg['uncertainty']['roundoff_guard_factor']==64. and cfg['uncertainty']['unit_roundoff_budget']==np.finfo(float).eps,'Fixed explicit binary64 uncertainty assumptions')
    need(cfg['native_sha256']==NATIVE_PINS['analysis/compute_s2_ed.jl'] and cfg['native_relative_path']=='analysis/compute_s2_ed.jl' and cfg['native_project_sha256']==NATIVE_PINS['Project.toml'] and cfg['native_manifest_sha256']==NATIVE_PINS['Manifest.toml'],'Unchanged source API/environment')
    need(cfg['science_cutoff_unix']==CUTOFF and cfg['pilot_wall_cap_seconds']==300 and cfg['production_wall_cap_seconds']==1200 and cfg['pilot_serialized_suffix_seconds']==1500 and cfg['termination_guard_seconds']==6 and cfg['outer_KILL_grace_seconds']==5 and cfg['cpus']==1 and cfg['memory_mb']==3500,'Fixed pilot/production cap, actual suffix and TERM6/KILL5 roles')
    need(tuple(cfg['required_previous_job_ids'])==PRIOR_JOBS and cfg['isaac_user']=='ywang287','Original twelve complete registered predecessors')
    need(cfg['initial_ensemble']=='Uniform complete global fixedN ground projector Pg/dg in everyN2..6 over allspinblocks' and cfg['fock_order']=='index0=up_mask*16+down_mask; sitej bitj; modes up0..up3,down0..down3','Complete canonical uniform ensemble and fermion ordering')
    return cfg


def audit_backend(audit,raw,label,cfg,freeze_sha,source,execution):
    meta=raw['metadata'];mode=execution['attempt_mode'];backend='native_ed_new_pair_adapter' if label=='target' else 'independent_quspin'
    need(raw['schema']==meta['schema']==RAW_SCHEMA and meta['catalogue_id']=='X16' and meta['release']==RELEASE and meta['backend']==backend and meta['attempt_mode']==mode and meta['pilot_excluded_from_production'] is (mode=='pilot') and meta['scientific_sampling'] is False,'Exact complete raw backend schema and role')
    need(meta['initial_ensemble']==cfg['initial_ensemble'] and meta['fock_order']==cfg['fock_order'],'Raw uniform global ground-ensemble/mode ordering')
    need(meta['source_freeze_sha256']==freeze_sha and meta['manifest_sha256']==sha(source/'manifest.toml') and meta['source_sha256']=={name:sha(source/name) for name in SOURCE_ROLES} and meta['native_sha256']==NATIVE_PINS['analysis/compute_s2_ed.jl'],'Raw frozen/native source identity')
    versions={'julia':cfg['julia_version']} if label=='target' else {p:cfg[p+'_version'] for p in ('quspin','numpy','scipy')};need(meta['versions']==versions,'Pinned backend runtime versions')
    if label=='target':need(meta['native_source_sha256']==NATIVE_PINS,'Actual native source-three pins')
    finished=raw['completion']['finished_unix'];process=next(r for r in execution['processes'] if r['name']==label)
    need(raw['completion']['status']=='COMPLETE' and type(finished) in (int,float) and math.isfinite(finished) and process['start_unix']-1<=finished<=process['end_unix']+1 and finished<CUTOFF,'Complete raw receipt inside actual owned process time')
    need(len(raw['site_annihilators'])==8 and len(raw['pair_annihilators'])==28,'Every site and wedge operator retained')
    own=[annihilator(4,m) for m in range(8)];saved=[sparse(c) for c in raw['site_annihilators']]
    for i,(a,b) in enumerate(zip(saved,own)):audit.qualified(label+'/C'+str(i)+' EXACT canonical integer operator',0. if np.array_equal(a,b) else 1.,0.)
    audit.qualified(label+' all same/crossspin CAR',car_error(saved),GATES['car_abs'])
    pairs=[]
    for i,(row,(a,b)) in enumerate(zip(raw['pair_annihilators'],pair_grid())):
        need(all(type(row[k]) is int for k in ('pair_index','a','b')) and (row['pair_index'],row['a'],row['b'])==(i,a,b),'All28 ordered wedge descriptors')
        P=sparse_product(own[b],own[a]);pairs.append(P);audit.qualified(label+'/pair'+str(i)+' EXACT wedge composition',0. if np.array_equal(sparse(row['matrix']),P) else 1.,0.)
    bridge=audit_quspin_bridge(audit,raw,own) if label=='reference' else None
    need([r['case_id'] for r in raw['cases']]==[c['case_id'] for c in cfg['cases']],'Complete ordered four-case grid')
    derived=[audit_case(audit,c,case,own,pairs,label,bridge) for c,case in zip(raw['cases'],cfg['cases'])]
    return derived


def cross_backends(audit,target,reference):
    rows=[];grounds=[]
    for a,b in zip(target,reference):
        case=a['case_id'];need(case==b['case_id'],'Matching full backend case order')
        for N in range(2,7):audit.compare(case+f'/N{N} all sector eigenvalues',maxabs(np.array(a['all_eigenvalues'][N])-np.array(b['all_eigenvalues'][N])),GATES['cross_energy_abs'])
        for x,y in zip(a['grounds'],b['grounds']):
            N=x['N'];prefix=case+f'/N{N}'
            audit.compare(prefix+' global multiplicity',0. if x['dimension']==y['dimension'] else 1.,0.)
            audit.compare(prefix+' full global ground projector',maxabs(x['projector']-y['projector']),GATES['cross_projector_abs'])
            audit.compare(prefix+' all784 complex pairdensity entries',maxabs(x['pair_density']-y['pair_density']),GATES['cross_pair_density_abs'])
            audit.compare(prefix+' complete physical ground interval',max(abs(x[k]-y[k]) for k in ('energy_physical','energy_lower','energy_upper')),GATES['cross_energy_abs'])
            audit.compare(prefix+' conditional ground interval overlap',max(0.,max(x['energy_lower'],y['energy_lower'])-min(x['energy_upper'],y['energy_upper'])),0.)
            grounds.append(dict(case_id=case,N=N,dimension=x['dimension'],energy_physical=x['energy_physical'],energy_interval=[x['energy_lower'],x['energy_upper']],pair_density_max_abs_error=maxabs(x['pair_density']-y['pair_density']),pairtrace=float(np.trace(x['pair_density']).real)))
        for x,y in zip(a['bindings'],b['bindings']):
            prefix=case+'/'+x['name'];audit.compare(prefix+' complete binding interval',max(abs(x[k]-y[k]) for k in ('value','lower','upper')),GATES['cross_binding_abs'])
            audit.compare(prefix+' identical binding sign classification',0. if x['classification']==y['classification'] else 1.,0.)
            audit.compare(prefix+' conditional binding interval overlap',max(0.,max(x['lower'],y['lower'])-min(x['upper'],y['upper'])),0.)
            rows.append(dict(case_id=case,reference_classification=y['classification'],**x))
    return rows,grounds


def scientific_status(audit,bindings):
    qualified=bool(audit.numerical) and all(r['passed'] for r in audit.numerical);matched=qualified and bool(audit.comparisons) and all(r['passed'] for r in audit.comparisons);resolved=qualified and len(bindings)==8 and all(r['classification'] in ('BOUND','NONBINDING') and r.get('reference_classification',r['classification']) in ('BOUND','NONBINDING') for r in bindings)
    decision='UNRESOLVED' if not audit.integrity() or not qualified or not resolved else 'PASS' if matched else 'FAIL'
    return decision,dict(schema='X16_NUMERICAL_PILOT_READINESS_V1',status='READY' if audit.integrity() and qualified and matched else 'NOT_READY',numerical_qualification=qualified,full_grid_cross_backend_match=matched,physical_signs_resolved=resolved,physical_sign_uncertainty_only_exempt=True,production_gates_unchanged=True)


def saved_comparison(audit,saved,mode,decision,ready,target,bindings,grounds):
    need(saved['schema']=='X16_PAIR_BINDING_COMPARISON_V1' and saved['catalogue_id']=='X16' and saved['attempt_mode']==mode and saved['scientific_sampling'] is False and saved['scientific_completions']==0 and saved['catalogue_status']==('PILOT_EXCLUDED' if mode=='pilot' else 'PENDING_INDEPENDENT_AUDIT'),'Separate excluded-pilot/pending-audit comparator scope')
    need('administrative_error' not in saved,'Error-only comparison is incomplete execution, never qualified scientific FAIL')
    audit.check('saved comparator matches independent whole-grid scientific decision',saved['scientific_status']==decision)
    for key in ('status','numerical_qualification','full_grid_cross_backend_match','physical_signs_resolved'):audit.check('saved comparator independent readiness '+key,saved['pilot_readiness'][key]==ready[key])
    if not ready['numerical_qualification']:
        need(saved['scientific_status']=='UNRESOLVED' and saved['cross_checks']==[] and saved['bindings']==[] and saved['ground_receipts']==[],'Numerically unqualified full raw result remains UNRESOLVED before cross-backend comparisons')
        return dict(qualified_cross_comparison_skipped=True,reason='Original numerical qualification failed; raw full-grid audit remains primary')
    need(saved['pilot_readiness']['production_gates_unchanged'] is True,'Production scientific gates remain fixed')
    need(saved['complete_cases']==4 and saved['spinblocks_per_backend']==76 and saved['fixedN_ground_ensembles']==20 and saved['pair_modes']==28 and saved['complex_pairdensity_entries']==15680,'Complete original raw scientific scope, no filtering')
    table=saved['ground_receipts']
    need(len(table)==20 and len({(r['case_id'],r['N']) for r in table})==20,'Unique complete20 comparator ground diagnostics')
    by={(r['case_id'],r['N']):r for r in table};summary_ground={}
    for case in target:
        summary_ground[case['case_id']]={}
        for own in case['grounds']:
            key=(case['case_id'],own['N']);need(key in by,'Every comparator ground coordinate');r=by[key];need(type(r['N']) is int and type(r['dimension']) is int and r['dimension']==own['dimension'],'Exact comparator ground dimension')
            numeric_list([r['energy_physical'],r['pair_density_max_abs_error'],r['pairtrace']]);numeric_list(r['energy_interval']);need(len(r['energy_interval'])==2,'Complete comparator conditional interval')
            scalar_formula(audit,str(key)+' saved physicalE',r['energy_physical'],own['energy_physical'],max(1.,abs(own['energy_mu']),abs(own['energy_physical'])))
            for i,bounds in enumerate((own['diagnostic_lower_range'],own['diagnostic_upper_range'])):audit.check(str(key)+' propagated norm-conditioning endpoint '+str(i),bounds[0]<=r['energy_interval'][i]<=bounds[1],dict(saved=r['energy_interval'][i],allowed_interval=bounds,scientific_gate_unchanged=True))
            need(r['energy_interval'][0]<=r['energy_physical']<=r['energy_interval'][1],'Saved ground interval ordered around its estimate')
            scalar_formula(audit,str(key)+' saved exact pairtrace',r['pairtrace'],own['N']*(own['N']-1)/2)
            audit.check(str(key)+' saved crossGamma diagnostic bounded by declared roundoff',0<=r['pair_density_max_abs_error'] and abs(r['pair_density_max_abs_error']-next(x['pair_density_max_abs_error'] for x in grounds if (x['case_id'],x['N'])==key))<=2*EPS64*256)
            summary_ground[case['case_id']][own['N']]=dict(N=own['N'],energy_physical=r['energy_physical'],energy_lower=r['energy_interval'][0],energy_upper=r['energy_interval'][1])
    need(len(saved['bindings'])==8 and [(r['case_id'],r['name']) for r in saved['bindings']]==[(r['case_id'],r['name']) for r in bindings],'All eight comparator binding coordinates')
    for r,own in zip(saved['bindings'],bindings):
        gm=summary_ground[r['case_id']];odd=own['odd_N'];lo,hi=own['even_N'];formula=binding_interval(own['name'],gm[odd],gm[lo],gm[hi]);scale=max(1.,2*abs(gm[odd]['energy_physical'])+abs(gm[lo]['energy_physical'])+abs(gm[hi]['energy_physical']))
        for key in ('value','lower','upper'):scalar_formula(audit,r['case_id']+'/'+r['name']+' exact comparator-own formula '+key,r[key],formula[key],scale)
        scalar_formula(audit,r['case_id']+'/'+r['name']+' exact comparator-own halfwidth',r['halfwidth'],(r['upper']-r['lower'])/2)
        scalar_formula(audit,r['case_id']+'/'+r['name']+' exact comparator-own roundoff',r['roundoff'],formula['roundoff'])
        audit.check(r['case_id']+'/'+r['name']+' saved independent/sign labels exact',r['classification']==formula['classification']==own['classification'])
        need(r['lower']<=r['value']<=r['upper'] and r['halfwidth']>=0,'Saved binding diagnostic ordered interval')
    return dict(primitive_and_own_formula_checks='strict64epsilon; independent matrix primitives use declared frozen norm-roundoff bounds',conditioned_endpoint_policy='Monotone propagation of two implementations primitive norm allowances through the fixed polar/Weyl formula; scientific point budgets/gates unchanged',all_labels_and_grid_exact=True)


def source_closure(audit,a):
    source=a.release_dir;freeze=audit.read(a.freeze)
    need(re.fullmatch('[0-9a-f]{64}',a.freeze_sha256) and sha(a.freeze)==a.freeze_sha256,'Explicit root-approved source freeze hash')
    need(freeze['release']==RELEASE and freeze['catalogue_id']=='X16' and set(freeze['files'])==set(SOURCE_ROLES)|{'SHA256SUMS'},'Exact canonical X16 source inventory')
    files=freeze['files'];actual={p.relative_to(source).as_posix() for p in source.rglob('*') if p.is_file()}
    need(actual==set(files)|{'source_freeze.json'} and not any(p.is_symlink() for p in source.rglob('*')),'Complete local source without extras/symlinks')
    for name,want in files.items():need(re.fullmatch('[0-9a-f]{64}',want) and sha(audit.file(source/safe_relative(name)))==want,'Exact frozen source role '+name)
    need(audit.file(source/'source_freeze.json').read_bytes()==a.freeze.read_bytes(),'Executed exact administrative source freeze')
    sums={}
    for line in (source/'SHA256SUMS').read_text(encoding='utf-8').splitlines():
        m=re.fullmatch(r'([0-9a-f]{64})  (.+)',line);need(m is not None and m[2] not in sums,'Canonical complete checksum manifest');sums[m[2]]=m[1]
    need(sums=={k:v for k,v in files.items() if k!='SHA256SUMS'},'Every frozen source byte covered once')
    review=audit.read(a.source_review);review_sha='6b2e80125097a45117e7df47a8833a6e9d3ff542e8722e0528ca2df208446b10'
    need(sha(a.source_review)==a.source_review_sha256==review_sha==freeze['independent_source_review_sha256'] and review['status']=='APPROVED_FOR_SOURCE_FREEZE_ONLY' and review['files']==sums,'Exact independently reviewed pre-freeze source13 snapshot')
    need(PurePosixPath(freeze['independent_source_review']).name==a.source_review.name,'Approved source peer identity')
    archive=audit.file(a.source_archive);need(sha(archive)==freeze['archive_sha256'] and archive.stat().st_size==freeze['archive_bytes'] and archive.name==PurePosixPath(freeze['archive']).name,'Source archive identity/length')
    verify_source_archive(archive,files)
    for name,want in NATIVE_PINS.items():need(sha(audit.file(a.native_source_dir/safe_relative(name)))==want,'Unchanged original native source/API environment '+name)
    cfg=protocol(audit.read(source/'manifest.toml'));catalogue=audit.read(a.catalogue);need(sha(a.catalogue)==CATALOGUE_SHA and isinstance(catalogue,list) and len(catalogue)==100,'Original100-question catalogue exact hash')
    literal=[r for r in catalogue if r.get('id')=='X16'];need(len(literal)==1 and literal[0]['acceptance']=='Propagate all energy uncertainties; binding only if full interval has the chosen sign.' and audit.read(source/'catalogue_X16.json')==literal[0],'Literal X16 uncertainty/sign question')
    provenance=audit.read(source/'PROVENANCE.json');need(provenance['schema']=='X16_PAIR_BINDING_ADAPTER_PROVENANCE_V1' and provenance['catalogue_id']=='X16' and provenance['catalogue_sha256']==CATALOGUE_SHA and provenance['native_source_sha256']==NATIVE_PINS and provenance['public_baseline_commit']=='196bb12570ed45b9d23926059dca2ef601aa5d07','Unchanged native/catalogue/public-baseline source provenance')
    action=provenance['validated_action_source'];need(action['path']=='x13_charged_spectral_v1/target.jl' and action['sha256']=='050479d10713b7c0832f7e8aff2cca9aeb72c702244e9dce92fd614823fc2cb6' and sha(audit.file(ROOT/safe_relative(action['path'])))==action['sha256'],'Previously reviewed source primitive closure; prior raw values excluded')
    prior=ROOT/'x13_charged_spectral_v1/source_freeze.json';need(sha(audit.file(prior))=='894ac44abd92d077d5dcc8fe07b205f78f440dc0cfda78150975e7f712ff7fb2','Prior source primitive freeze unchanged')
    for name in files:need(b'\r' not in (source/name).read_bytes(),'Raw LF-only source/config/checksum '+name)
    audit.check('complete independently approved source/archive/native/catalogue closure before raw math',True)
    return cfg,freeze


def transfer_closure(audit,a,freeze):
    directory=a.job_dir;record=audit.read(directory/'transfer.json')
    need(str(record['job_id'])==a.job_id and record.get('transfer_status')=='COMPLETE' and record.get('terminal_snapshot') is True and record.get('two_pass_inventory_and_bytes_verified') is True and record.get('all_local_payload_bytes_and_hashes_verified') is True and record.get('remote_writes') is False and record.get('binaries_copied') is False and record.get('missing_paths')==[],'Complete stable original-byte terminal transfer')
    roles=set();aliases=set()
    for row in record['files']:
        role=row['local'];safe_relative(role);need(role.casefold() not in aliases and row['status']=='RETRIEVED','Unique retrieved original transfer role');aliases.add(role.casefold());p=audit.file(directory/role);need(p.is_relative_to(directory) and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'Exact original transferred bytes '+role);roles.add(role)
    need({'results/'+n for n in OUTPUT_ROLES}<=roles,'All21 original complete result roles mandatory')
    for folder,expected in [('source',dict(freeze['files'],**{'source_freeze.json':a.freeze_sha256})),('native_source',NATIVE_PINS)]:
        actual={n[len(folder)+1:] for n in roles if n.startswith(folder+'/')};need(actual==set(expected),'Exact transferred '+folder+' inventory')
        for name,want in expected.items():need(sha(audit.file(directory/folder/name))==want,'Original copied approved '+folder+' role '+name)
    need(record['source_root']==REMOTE_ROOT+'/source/manual_pc_continuation/'+RELEASE and record['result_root']==REMOTE_ROOT+'/results/'+RELEASE+'_'+a.attempt_mode+'_'+a.job_id,'Exact executed source/result namespace')
    scheduler=audit.file(directory/'scheduler.txt');need(scheduler.read_text(encoding='utf-8-sig')==record['sacct'],'Original scheduler/transfer exact text binding')
    rows=list(csv.DictReader(io.StringIO(record['sacct']),delimiter='|'));selected=[]
    for key in (a.job_id,a.job_id+'.batch'):
        candidates=[r for r in rows if r['JobIDRaw']==key];need(len(candidates)==1,'Unique own main/batch scheduler receipt');selected.append(candidates[0])
    main,batch=selected;cap=300 if a.attempt_mode=='pilot' else 1200
    need(main['JobName']=='x16_'+a.attempt_mode and main['ReqCPUS']==main['AllocCPUS']=='1' and main['ReqMem']=='3500M' and duration_seconds(main['Timelimit'])==cap,'Fixed owned attempt name/allocation/wall')
    need(all(r['State'].split()[0].rstrip('+') in TERMINAL and scheduler_epoch(r['End'])<=CUTOFF and 0<duration_seconds(r['Elapsed'])<=cap+1 for r in selected),'Owned terminal timing/cutoff')
    text=audit.file(directory/'results/scontrol.txt').read_text(encoding='utf-8');fields={}
    for m in re.finditer(r'(?:^|\s)([A-Za-z][A-Za-z0-9_]*)=([^\s]+)',text):
        k,v=m.groups()
        if k in ('JobId','JobName','UserId','Command','WorkDir','StdOut','Restarts','Requeue'):need(k not in fields,'Duplicate executed scontrol identity');fields[k]=v
    source=REMOTE_ROOT+'/source/manual_pc_continuation/'+RELEASE
    need(fields.get('JobId')==a.job_id and fields.get('JobName')=='x16_'+a.attempt_mode and re.fullmatch(r'ywang287(?:\([0-9]+\))?',fields.get('UserId','')) is not None and fields.get('Command')==source+'/'+a.attempt_mode+'.sbatch' and fields.get('WorkDir')==source and fields.get('StdOut')==REMOTE_ROOT+'/logs/'+RELEASE+'_'+a.attempt_mode+'_'+a.job_id+'.log' and fields.get('Requeue')=='0' and fields.get('Restarts')=='0','Executed original wrapper/owner/workdir/log identity')
    return record,selected


def active_queue(text,job,mode):
    need(isinstance(text,str) and text.strip() and '\0' not in text,'Nonempty original active visibility')
    rows=[];seen=set()
    for line in text.splitlines():
        fields=[p.strip() for p in line.split('|')];need(len(fields)==4,'Complete active queue fields');identity,state,cpus,name=fields
        need(re.fullmatch('[0-9]+(?:_[0-9]+)?',identity) and identity not in seen and state in ACTIVE and re.fullmatch('[1-9][0-9]*',cpus) and name and all(ord(c)>=32 for c in name),'Strict unique active queue row');seen.add(identity)
        need(identity==job or name not in ('x16_pilot','x16_production'),'No other X16 lane active');rows.append(dict(job_id=identity,state=state,cpus=int(cpus),job_name=name))
    own=[r for r in rows if r['job_id']==job];need(1<=len(rows)<=2 and sum(r['cpus'] for r in rows)<=2 and len(own)==1 and own[0]['cpus']==1 and own[0]['job_name']=='x16_'+mode,'Whole-user two-job/two-CPU and exact own oneCPU identity')
    return dict(schema='X16_ACTIVE_QUEUE_GUARD_V1',status='PASS',rows=rows,active_jobs=len(rows),active_cpus=sum(r['cpus'] for r in rows))


def previous_terminal(text,batch_start,cancelled_proof=None):
    reader=csv.DictReader(io.StringIO(text),delimiter='|');needed={'JobIDRaw','State','User','AllocCPUS','Start','End','ExitCode'};need(reader.fieldnames and len(set(reader.fieldnames))==len(reader.fieldnames) and needed<=set(reader.fieldnames),'Original full prior scheduler fields')
    rows=list(reader);selected=[]
    for job in PRIOR_JOBS:
        found=[r for r in rows if r['JobIDRaw']==job];need(len(found)==1,'Unique original predecessor main row');r=found[0]
        need(r['User']=='ywang287' and r['State'].split()[0].rstrip('+') in TERMINAL and re.fullmatch('[0-9]+',r['AllocCPUS']) and re.fullmatch('[0-9]+:[0-9]+',r['ExitCode']),'Original predecessor owned and terminal')
        if r['Start']=='None':
            need(job in CANCELLED_BEFORE_START and cancelled_proof is not None,'Only two explicitly root-proved never-allocated cancelled predecessors may omit Start')
            expected=cancelled_proof['terminal_stage3_jobs'][job]
            need(r['State']=='CANCELLED by 12640' and r['AllocCPUS']=='0' and r['ExitCode']=='0:0' and r['End']==CANCELLED_BEFORE_START[job] and expected['ElapsedRaw']=='0' and all(r[k]==expected[k]for k in ('JobIDRaw','State','User','AllocCPUS','Start','End','ExitCode')),'Exact original zero-allocation cancelled predecessor matches pinned root closure')
            need(scheduler_epoch(r['End'])<=batch_start,'Original cancelled predecessor End before actual X16 batch')
        else:need(scheduler_epoch(r['Start'])<=scheduler_epoch(r['End'])<=batch_start+1,'Allocated original predecessor strict start/end clock')
        selected.append(r)
    need(all(None not in r and set(r)==set(reader.fieldnames) and (r['JobIDRaw'] in PRIOR_JOBS or r['JobIDRaw'].split('.')[0] in PRIOR_JOBS and re.fullmatch(r'[0-9]+\.(?:batch|extern|[0-9]+)',r['JobIDRaw'])) for r in rows),'No malformed or unexpected predecessor query rows')
    return dict(schema='X16_PRIOR_TWELVE_TERMINAL_GUARD_V1',status='PASS',rows=selected)


def runtime_guard(guard,cfg,mode,job,queue,previous,cancelled_proof=None):
    cap=300 if mode=='pilot' else 1200;suffix=1500 if mode=='pilot' else 1200
    start,checked=guard['actual_batch_start_unix'],guard['checked_unix'];need(all(type(v) in (int,float) and math.isfinite(v) for v in (start,checked)) and 0<start<=checked<start+cap-6 and start+suffix<=CUTOFF,'Full actual-start cap/suffix before science cutoff')
    need(guard['schema']=='X16_ACTUAL_START_ADMISSION_V1' and guard['status']=='PASS' and guard['attempt_mode']==mode and guard['wall_cap_seconds']==cap and guard['serialized_suffix_seconds']==suffix and guard['cap_end_unix']==start+cap and guard['science_cutoff_unix']==CUTOFF and guard['termination_guard_seconds']==6 and guard['outer_KILL_grace_seconds']==5,'Exact TERM6/KILL5/full cap/suffix guard schema')
    need(guard['active_queue_guard']==active_queue(queue,job,mode) and guard['previous_terminal_guard']==previous_terminal(previous,start,cancelled_proof),'Independently reconstructed original queue/predecessor guard')
    return start+cap-6


def resource_receipt(audit,path,process,cap):
    text=audit.file(path).read_text(encoding='utf-8');fields={}
    for line in text.splitlines():
        if '): ' in line:k,v=line.strip().rsplit('): ',1);fields[k+')']=v
        elif ': ' in line:k,v=line.strip().split(': ',1);fields[k]=v
    wall=duration_seconds(fields['Elapsed (wall clock) time (h:mm:ss or m:ss)']);rss=int(fields['Maximum resident set size (kbytes)']);code=int(fields['Exit status'])
    need(0<=wall<=cap and 0<=rss<=3500*1024 and code==process['exit_code'] and abs(wall-(process['end_unix']-process['start_unix']))<=2,'Original GNU time wall/RSS/exit linked to actual process')
    return dict(wall_seconds=wall,maximum_rss_kbytes=rss,exit_code=code)


def preliminary_execution(audit,a,cfg,freeze,scheduler_rows):
    result=a.job_dir/'results';pre=audit.read(result/'preflight.json');ex=audit.read(result/'execution.json');mode=a.attempt_mode;job=a.job_id
    need(ex['schema']==pre['schema']=='X16_EXECUTION_RECEIPT_V1' and ex['catalogue_id']==pre['catalogue_id']=='X16' and ex['job_id']==pre['job_id']==job and ex['release']==pre['release']==RELEASE and ex['attempt_mode']==pre['attempt_mode']==mode and ex['pilot_excluded_from_production'] is (mode=='pilot'),'Owned exact pilot/production execution receipt')
    immutable=('source_freeze_sha256','source_sha256','native_source_sha256','actual_batch_start_unix','worker_start_unix','wall_cap_seconds','serialized_suffix_seconds','cap_end_unix','science_cutoff_unix','termination_guard_seconds','outer_KILL_grace_seconds','cpus','memory_mb','scientific_sampling','scientific_completions','new_adapter','process_group_id','actual_start_admission','threads')
    need(all(pre[k]==ex[k] for k in immutable),'Preflight immutable admission/source fields retained')
    need(ex['status']=='PROCESSES_COMPLETE' and pre['status']=='STARTED' and pre['processes']==[] and pre['process_exit_codes']=={} and ex['source_freeze_sha256']==a.freeze_sha256 and ex['source_sha256']==freeze['files'] and ex['native_source_sha256']==NATIVE_PINS and ex['source_unchanged'] is True and ex['native_unchanged'] is True,'Complete processes with unchanged full source/native closure')
    need(ex['cpus']==1 and ex['memory_mb']==3500 and ex['scientific_sampling'] is False and ex['scientific_completions']==0 and ex['new_adapter'] is True and ex['threads']=={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','JULIA_NUM_THREADS')},'Fixed process resource/scope/thread controls')
    queue=audit.file(result/'active_queue.txt').read_text(encoding='utf-8');prior=audit.file(result/'previous_twelve_terminal.txt').read_text(encoding='utf-8');limit=runtime_guard(ex['actual_start_admission'],cfg,mode,job,queue,prior,getattr(a,'_cancelled_predecessor_proof',None))
    need(ex['actual_batch_start_unix']==ex['actual_start_admission']['actual_batch_start_unix'] and ex['cap_end_unix']==ex['actual_start_admission']['cap_end_unix'],'Batch epoch/cap linked to actual guard')
    need([r['name'] for r in ex['processes']]==['target','reference','comparator'] and ex['process_exit_codes']['target']==ex['process_exit_codes']['reference']==0 and ex['process_exit_codes']['comparator'] in (0,1,2),'All native/reference/comparator processes completed with interpretable exit codes')
    need(type(ex['process_group_id']) is int and ex['process_group_id']>0 and ex['process_exit_codes']=={r['name']:r['exit_code'] for r in ex['processes']},'Exact original inherited process group and process exit vector')
    bind_execution_clock(ex,scheduler_rows)
    need(all(ex['actual_batch_start_unix']<=p['start_unix']<p['end_unix']<=limit and p['process_group_inherited'] is True for p in ex['processes']) and all(a['end_unix']<=b['start_unix'] for a,b in zip(ex['processes'],ex['processes'][1:])) and ex['processes'][-1]['end_unix']<=ex['finished_unix']<=limit,'Serial process/receipt completion before TERM reserve inside the inherited group')
    resources={}
    source=REMOTE_ROOT+'/source/manual_pc_continuation/'+RELEASE;remote_result=REMOTE_ROOT+'/results/'+RELEASE+'_'+mode+'_'+job
    for p in ex['processes']:
        name=p['name'];need(audit.read(result/(name+'.process.json'))==p,'Exact original per-process receipt '+name)
        command=p['command'];need(isinstance(command,list) and all(isinstance(x,str) for x in command),'Original exact process argv list')
        need(source+'/'+('target.jl' if name=='target' else 'reference.py' if name=='reference' else 'compare.py') in command and command.count('--manifest')==command.count('--freeze-sha256')==command.count('--output')==1 and command[command.index('--manifest')+1]==source+'/manifest.toml' and command[command.index('--freeze-sha256')+1]==a.freeze_sha256 and command[command.index('--output')+1]==remote_result+'/'+('target.toml' if name=='target' else 'reference.json' if name=='reference' else 'comparison.json'),'Executed source/freeze/output argv identity '+name)
        if name=='target':need(command[:4]==[REMOTE_ROOT+'/runtime/julia-1.12.1/bin/julia','--startup-file=no','--history-file=no','--project='+REMOTE_ROOT+'/source/renyi_disjoint_corrected'],'Pinned Julia native runtime/project argv')
        else:need(command[:2]==[REMOTE_ROOT+'/quspin-venv/bin/python','-B'],'Pinned independent reference/comparator Python runtime')
        if name=='comparator':need(command[command.index('--target')+1]==remote_result+'/target.toml' and command[command.index('--reference')+1]==remote_result+'/reference.json','Exact original raw inputs to frozen comparator')
        resources[name]=resource_receipt(audit,result/(name+'.resources.txt'),p,300 if mode=='pilot' else 1200)
    exits={}
    for line in audit.file(result/'process_exit_codes.txt').read_text(encoding='ascii').splitlines():
        k,v=line.split('=',1);need(k not in exits,'Unique process exit roles');exits[k]=v
    need(exits=={k+'_exit':str(v) for k,v in ex['process_exit_codes'].items()},'Original process exit ledger exact')
    outputs={name:sha(audit.file(result/name)) for name in OUTPUT_ROLES};need(ex['output_sha256']=={k:v for k,v in outputs.items() if k!='execution.json'},'Complete20 immutable pre-execution-output hashes and final receipt')
    audit.check('complete source-bound native/reference/process/resource/clock proof',True)
    return ex,outputs,resources


def registry_object_fragments(raw):
    text=raw.decode('utf-8-sig');whole=json.loads(text,object_pairs_hook=json_object);need(isinstance(whole,dict) and isinstance(whole.get('jobs'),list),'Whole original registry object')
    decoder=json.JSONDecoder();index=next(i for i,c in enumerate(text) if not c.isspace())+1
    def space(i):
        while i<len(text) and text[i] in ' \r\n\t':i+=1
        return i
    while True:
        index=space(index)
        if text[index]==',':index=space(index+1)
        key,index=decoder.raw_decode(text,index);index=space(index);need(text[index]==':','Registry root key separator');index=space(index+1)
        if key=='jobs':break
        _,index=decoder.raw_decode(text,index)
    need(text[index]=='[','Literal top-level jobs array');index+=1;out=[]
    while True:
        index=space(index)
        if text[index]==',':index=space(index+1)
        if text[index]==']':break
        begin=index;row,index=decoder.raw_decode(text,index);need(isinstance(row,dict),'Registry job object');out.append((row,text[begin:index].encode('utf-8')))
    need([r[0] for r in out]==whole['jobs'],'Exact original literal registry objects')
    return out


def prelaunch_queue(text,known):
    need(isinstance(text,str) and '\0' not in text,'Original complete pre-submission visibility');active=[];seen=set()
    for line in text.splitlines():
        if not line.strip():continue
        row=[p.strip() for p in line.split('|')];need(len(row)==5 and re.fullmatch('[0-9]+(?:_[0-9]+)?',row[0]) and row[0] not in seen and row[1] in ACTIVE|{'PENDING'} and re.fullmatch('[1-9][0-9]*',row[2]),'Strict unique prelaunch queue row');seen.add(row[0]);need(row[0] not in known,'Required terminal predecessor still queued')
        if row[1] in ACTIVE:need(row[3] not in ('x16_pilot','x16_production'),'Another X16 science lane already active');active.append(row)
    need(len(active)<=1 and sum(int(r[2]) for r in active)<=1,'One prospective oneCPU slot before sbatch')
    return len(active),sum(int(r[2]) for r in active)


def remote_submission_closure(audit,a,config,receipt,known):
    folder=a.job_dir/'submissions';intent=audit.read(folder/'intent.json');durable=audit.read(folder/'durable_before_sbatch.json');batch=audit.read(folder/'sbatch_result.json');original=audit.read(folder/'receipt.json');plan=dict(intent['plan']);initial=plan.pop('additional_live_guards')
    need(plan==config and original==receipt and original['job_id']==a.job_id,'Original transferred remote plan/acknowledgement exactly matches exclusive local phase')
    need(durable['intent_sha256']==sha(folder/'intent.json') and durable['maximum_sbatch_calls']==1 and durable['argv']==batch['argv']==receipt['argv'] and all(batch[k]==receipt[k] for k in ('returncode','stdout','stderr','utc')) and intent['queue_before']==receipt['queue_before'],'Durable original remote intent linked to exactly one original sbatch result')
    prelaunch_queue(intent['queue_before'],known);submitted=scheduler_epoch(receipt['utc'])
    for label,proof in [('initial',initial),('durable',durable['live_guard'])]:
        counts=prelaunch_queue(proof['queue'],known);reader=csv.DictReader(io.StringIO(proof['initial_accounting']),delimiter='|');rows=list(reader)
        for job in known:
            selected=[r for r in rows if r['JobIDRaw']==job];need(len(selected)==1 and selected[0]['User']=='ywang287' and selected[0]['State'].split()[0].rstrip('+') in TERMINAL,'Fresh exact owned predecessor terminal before original sbatch')
            if job not in PRIOR_JOBS:need(selected[0]['JobName']=='x16_pilot' and selected[0]['State']=='COMPLETED' and selected[0]['ExitCode']=='0:0' and selected[0]['AllocCPUS']=='1','Production freshly qualified its exact pilot COMPLETED0')
        need(type(proof['checked_unix']) in (int,float) and math.isfinite(proof['checked_unix']) and 0<proof['checked_unix']<=submitted+1 and proof['checked_unix']+config['full_suffix_seconds']<=CUTOFF and proof['original_twelve_terminal'] is True and (proof['strict_active_jobs'],proof['strict_active_cpus'])==counts and proof['prospective_jobs']==counts[0]+1 and proof['prospective_cpus']==counts[1]+1,'Original fresh whole-user/predecessor/full-suffix live gate '+label)
    need(initial['checked_unix']<=durable['live_guard']['checked_unix']<=submitted+1,'Prelaunch guards precede original submission')
    argv=receipt['argv'];dependencies=config['dependencies'];need(argv[0]=='sbatch' and '--no-requeue' in argv and '--cpus-per-task=1' in argv and '--mem=3500M' in argv and '--time='+str(config['wall_minutes']) in argv and argv[-1]==config['source']+'/'+config['wrapper'] and [s for s in argv if s.startswith('--dependency=')]==(['--dependency=afterok:'+dependencies[0]] if dependencies else []),'One exact allocated/dependent original sbatch argv')
    for row in config['input_files']:
        name=str(safe_relative(row['name']));need('/' not in name and audit.file(folder/name).read_bytes()==row['content'].encode('utf-8'),'Exact original staged qualification alias bytes '+name)
    return submitted


def checkpoint_roster(checkpoint,freeze_sha,checksum_sha):
    """The original nine releases remain, and X16 is explicitly the tenth."""
    expected={'continuum_production_v1b':'69981840a5b9c4a2280f180655abb3e2fe9710381830440220d28b3f09dda0de','p21_sampling_pilot_v1':'1dc2c0daa3dc273bb03758a44b13f0dbf4f68dbf335c85654c44fb94437e4879','continuum_production_v1c':'5caa7b1e27b793909e449d55c7ae4e91cc2d868a5124ddfef1cbd87929298e06','x01_continuum_stage2_v1':'d61f31f2849da116c03b46e1098ce3f245d490d34045f71d534671e175e316a0','x01_continuum_stage2_v2':'72f01b0ab2b53c18a9de6b3079b96954191f1eaeb82004e7b3f62608e90632e3','p21_sampled_production_v1':'1c11152366e1802ba2a650550d52a97e0e758ab08b2c2a1498f3bdf89a2bb48a','x13_charged_spectral_v1':'3b5d399c35940dc3bf7571b59376d62c42e96615f341181ae189776174cf2bde','x01_stage2_analysis_repair_v1':'b3cd596752091780ba240e0914ef766a998fb0be2e04f9b911379f4322d22182','x01_continuum_stage3_v1':'b851d498e72e8a3ba22e0d7f7b27af73136fd3e9923a19997a6349b8ee86dab0',RELEASE:checksum_sha}
    rows=keyed(checkpoint['frozen_releases'],'directory',set(expected))
    need(all(rows[name]['manifest_sha256']==want for name,want in expected.items()),'Original nine checksum-manifest pins and exact X16 tenth checksum manifest preserved')
    need(rows[RELEASE]['freeze_record']==RELEASE+'/source_freeze.json' and rows[RELEASE]['freeze_record_sha256']==freeze_sha,'Initial checkpoint explicitly binds active X16 source freeze before pilot admission')
    return True


def acknowledged_transition(before,after,phase,document):
    """Only the named reserved row may change, and only to its exact receipt."""
    old,new=registry_object_fragments(before),registry_object_fragments(after);index=12+(phase=='production');other=25-index
    need(len(old)==len(new)==14 and all(old[i]==new[i] for i in range(14) if i!=index),'Acknowledgement preserves every other registry object AND literal byte fragment')
    receipt=document['receipt'];config=document['configuration'];job=str(document['complete']['job_id'])
    expected=dict(old[index][0],job_id=job,state='SUBMITTED',submitted_utc=receipt['utc'],dependencies=config['dependencies'],result=receipt['result'],result_prefix=config['result_prefix'],receipt='x16_two_phase_transaction_v1/'+phase+'/remote_receipt.json',log_path=REMOTE_ROOT+'/logs/'+config['result_prefix']+job+'.log')
    need(new[index][0]==expected and receipt['job_id']==job,'Acknowledged reserved row equals exact phase receipt and no other field changes')
    if phase=='production':need(new[other][0]['job_id']!=job and new[other][0]['state']=='SUBMITTED','Production preserves a distinct already acknowledged pilot')
    return True


def mechanical_authority(audit,a):
    """Root-pinned postpilot metadata repair, preserving original preapproval."""
    record=audit.read(a.mechanical_amendment);need(re.fullmatch('[0-9a-f]{64}',a.mechanical_amendment_sha256) and sha(a.mechanical_amendment)==a.mechanical_amendment_sha256,'Explicit external root mechanical-amendment SHA')
    need(record['schema']=='ISAAC_X16_POSTPILOT_MECHANICAL_AUDITOR_AMENDMENT_V1' and record['status']=='ROOT_APPROVED_FOR_METADATA_ONLY_POSTPILOT_AUDIT' and record['pilot_job_id']==ORIGINAL_PILOT and record['source_freeze_sha256']==a.freeze_sha256,'Explicit audit-only metadata repair for this original pilot/source')
    need(all(record[k]is True for k in ('original_pilot_preapproval_unchanged','no_new_compute','scientific_gates_unchanged')) and record['new_source_preapproved_before_pilot']is False,'Replacement authority honestly postpilot; original compute/preapproval unchanged')
    def artifact(row):
        path=ROOT/safe_relative(row['path']);need(sha(audit.file(path))==row['sha256'],'Exact root amendment artifact '+row['path']);return audit.read(path) if path.suffix=='.json' else path
    original=record['original_auditor'];need(original['path']=='audit_x16_pair_binding_v2.py' and original['sha256']==ORIGINAL_AUDITOR_SHA,'Preserved original approved5c source')
    artifact(original)
    for role in ('peer_review','validation','tests','root_approval'):artifact(original[role])
    workflow=artifact(record['original_admission_root_approval']);need(record['original_admission_root_approval']['sha256']==ORIGINAL_WORKFLOW_SHA and workflow['independent_auditor']==original and workflow['source_freeze']['sha256']==a.freeze_sha256,'Preserved originala064 workflow and full old auditor approval')
    replacement=record['replacement_auditor'];need(replacement['path']==Path(__file__).name and replacement['sha256']==sha(__file__),'Only exact newly reviewed replacement verifier')
    artifact(replacement);peer=artifact(replacement['peer_review']);validation=artifact(replacement['validation']);artifact(replacement['tests'])
    need(peer['status']=='APPROVED_FOR_LOCAL_RECEIPT_AUDITOR_SOURCE_ONLY' and peer['auditor_sha256']==replacement['sha256'] and peer['test_sha256']==replacement['tests']['sha256'],'Independent source/test peer for exact replacement')
    need(validation['status']=='PASS' and validation['returncode']==0 and validation['files'][replacement['path']]==replacement['sha256'] and validation['files'][Path(replacement['tests']['path']).name]==replacement['tests']['sha256'] and validation['all_scientific_function_class_ASTs_identical']is True and validation['scientific_gate_changes']is False and validation['source_freeze_changes']is False,'Replacement validated as administrative-only with exact science AST/gates/source')
    failed=artifact(record['original_failed_audit']);need(record['original_failed_audit']['sha256']==ORIGINAL_FAILED_AUDIT_SHA and failed['schema']==SCHEMA and failed['job_id']==ORIGINAL_PILOT and failed['attempt_mode']=='pilot' and failed['audit_tool_sha256']==ORIGINAL_AUDITOR_SHA and failed['integrity_status']=='UNRESOLVED' and failed['check_count']==159 and failed['numerical_qualification_checks']==[] and failed['cross_backend_checks']==[],'Original failed pre-math159-check audit remains exact evidence')
    for path,want in failed['input_sha256'].items():need(sha(audit.file(path))==want,'Every original failed audit input remains unchanged')
    transfer=artifact(record['pilot_transfer']);scheduler=artifact(record['pilot_scheduler']);execution=artifact(record['pilot_execution']);registry=artifact(record['pilot_acknowledged_registry'])
    need(record['pilot_transfer']['sha256']=='591d248a1bf894e8fe345c197ce2a074d3cd9df80e8f89c7663cdb1180f562e5' and record['pilot_scheduler']['sha256']=='757d9b0efe86e7a2ca46d21a104312975ef30e0a9259c8a99e40f07e133bb321' and record['pilot_acknowledged_registry']['sha256']=='63bf643a87eac304c82f94d01bb94ecaa3559db2cde29ac3f915212f322fb6a4' and transfer['job_id']==execution['job_id']==ORIGINAL_PILOT and execution['source_freeze_sha256']==a.freeze_sha256 and registry['jobs'][12]['job_id']==ORIGINAL_PILOT,'Exact original pilot receipts and acknowledged registry bound')
    pilot_parent=(ROOT/safe_relative(record['pilot_transfer']['path'])).parent
    need(ROOT/safe_relative(record['pilot_scheduler']['path'])==pilot_parent/'scheduler.txt' and ROOT/safe_relative(record['pilot_execution']['path'])==pilot_parent/'results/execution.json' and ROOT/safe_relative(record['original_failed_audit']['path'])==pilot_parent/'independent_audit_v1.json','All original failed-audit/pilot receipt roles derive from the exact original transfer parent')
    execution_rows=[r for r in transfer['files']if r['local']=='results/execution.json'];need(len(execution_rows)==1 and execution_rows[0]['sha256']==record['pilot_execution']['sha256'],'Mechanical authority uses the exact original transferred pilot execution receipt')
    closure=artifact(record['stage3_terminal_closure']);need(record['stage3_terminal_closure']['sha256']==STAGE3_CLOSURE_SHA and closure['schema']=='ISAAC_X01_STAGE3_ROOT_TERMINAL_CLOSURE_V1' and closure['status']=='ROOT_CLOSED_EXECUTION_INCOMPLETE_WITHOUT_LEDGER_REVISION' and closure['blocked_dependents_cancelled']==list(CANCELLED_BEFORE_START) and closure['cancelled_before_start']is True and closure['cancelled_allocation_seconds']==0 and closure['new_catalogue_completions']==0 and closure['scientific_ledger_changed']is False and closure['fourth_look_authorized']is False,'Exact root proof of original zero-allocation cancellations and no promotion')
    created=dt.datetime.fromisoformat(record['created_utc'].replace('Z','+00:00'));failed_time=dt.datetime.fromisoformat(failed['created_utc'].replace('Z','+00:00'));need(created.tzinfo is not None and failed_time.tzinfo is not None and max(execution['finished_unix'],failed_time.timestamp())<=created.timestamp()<=dt.datetime.now(dt.timezone.utc).timestamp(),'New metadata-only authority was approved after original pilot/failed audit and before this local audit')
    a._original_auditor_tool=original;a._cancelled_predecessor_proof=closure;a._mechanical_amendment_record=record
    audit.check('postpilot metadata repair has separate explicit root authority; original preapproved5c/a064 retained',True)
    need(a.mechanical_amendment.is_relative_to(ROOT),'Mechanical root amendment remains BASE-local immutable evidence')
    proof=dict(schema='X16_POSTPILOT_MECHANICAL_AUDITOR_AUTHORITY_V1',amendment=dict(path=a.mechanical_amendment.relative_to(ROOT).as_posix(),sha256=a.mechanical_amendment_sha256),original_auditor_sha256=ORIGINAL_AUDITOR_SHA,current_auditor_sha256=sha(__file__),original_pilot_preapproval_unchanged=True,replacement_preapproved_before_pilot=False,no_new_compute=True,original_failed_audit_sha256=ORIGINAL_FAILED_AUDIT_SHA,source_freeze_sha256=a.freeze_sha256,pilot_job_id=ORIGINAL_PILOT,postpilot_metadata_repair=True,scientific_gates_unchanged=True)
    a._mechanical_authority_proof=proof;return proof


def production_controller_authority(audit,a,policy,config,found):
    """Explicit new outer-controller authority; original d39 remains delegated."""
    need(a.production_controller_amendment is not None and re.fullmatch('[0-9a-f]{64}',a.production_controller_amendment_sha256 or ''),'Explicit external production-controller amendment required')
    need(sha(audit.file(a.production_controller_amendment))==a.production_controller_amendment_sha256,'External root controller amendment hash is exact')
    readiness=policy['x16_production_readiness'];mechanical=a._mechanical_authority_proof['amendment'];controller_pin=dict(path=a.production_controller_amendment.relative_to(ROOT).as_posix(),sha256=a.production_controller_amendment_sha256)
    need(readiness['auditor_mechanical_amendment']==mechanical and readiness['production_controller_amendment']==controller_pin,'Only explicit readiness addition licenses both separate postpilot amendments')
    def artifact(row):
        original=str((ROOT/safe_relative(row['path'])).resolve());need(original in found,'Required controller artifact archived before production');path=audit.file(found[original]);need(sha(path)==row['sha256'],'Exact original archived controller/amendment bytes '+row['path']);return audit.read(path) if Path(row['path']).suffix=='.json' else path
    record=artifact(controller_pin);need(record['schema']=='ISAAC_X16_POSTPILOT_PRODUCTION_CONTROLLER_AMENDMENT_V1' and record['status']=='ROOT_APPROVED_FOR_ONE_FIXED_PRODUCTION_AFTER_METADATA_REPAIR' and record['original_admission_root_approval']['sha256']==ORIGINAL_WORKFLOW_SHA and record['original_auditor']==a._original_auditor_tool and record['pilot_job_id']==ORIGINAL_PILOT,'Separate root amendment preserves original workflow/pilot/verifier identity')
    need(record['auditor_mechanical_amendment']==mechanical and record['replacement_auditor']==a._mechanical_amendment_record['replacement_auditor'] and record['source_freeze']['sha256']==a.freeze_sha256 and record['stage3_terminal_closure']['sha256']==STAGE3_CLOSURE_SHA,'Outer production controller binds exact approved audit-only repair/source/closure')
    need(all(record[k]is True for k in ('original_admission_preapproval_unchanged','scientific_source_and_gates_unchanged','no_pilot_rerun','reuse_existing_production_reservation')) and record['cpus']==1 and record['memory_mb']==3500 and record['wall_minutes']==20 and record['full_suffix_seconds']==1200 and record['science_cutoff_unix']==CUTOFF and record['maximum_sbatch_calls']==1,'Original production cap/reservation/suffix/no-rerun/single-sbatch controls preserved')
    need(record['ready_pilot_audit']==readiness['artifacts']['independent_audit'] and record['current_checkpoint']==readiness['current_checkpoint'],'Exact approved repaired pilot audit and current production checkpoint')
    artifact(record['original_admission_root_approval']);artifact(record['original_admission_helper']);artifact(record['pilot_acknowledged_registry']);artifact(record['pilot_admission_intent']);artifact(record['auditor_mechanical_amendment']);artifact(record['ready_pilot_audit']);artifact(record['current_checkpoint']);artifact(record['stage3_terminal_closure']);artifact(record['source_freeze'])
    need(record['pilot_acknowledged_registry']['sha256']=='63bf643a87eac304c82f94d01bb94ecaa3559db2cde29ac3f915212f322fb6a4' and record['pilot_admission_intent']['sha256']==sha(a.admission_dir/'pilot/intent.json'),'Original pilot ACK/intent remain immutable')
    original_helper=record['original_admission_helper'];need(original_helper['path']=='overnight_20261022/submit_x16_after_queue_candidate_v3.py' and original_helper['sha256']==policy['x16_admission']['submission_helper_sha256'],'Exact original d39 controller remains delegated, no fiction new preapproval')
    controller=record['controller'];artifact(controller);peer=artifact(controller['peer_review']);validation=artifact(controller['validation']);artifact(controller['tests'])
    need(peer['status']=='APPROVED_FOR_X16_POSTPILOT_PRODUCTION_CONTROLLER_SOURCE_ONLY' and peer['controller_sha256']==controller['sha256'] and peer['test_sha256']==controller['tests']['sha256'] and validation['status']=='PASS' and validation['returncode']==0 and validation['files'][Path(controller['path']).name]==controller['sha256'] and validation['files'][Path(controller['tests']['path']).name]==controller['tests']['sha256'],'Independent source/fixture approval for exact new outer production controller')
    expected=dict(schema='ISAAC_X16_POSTPILOT_PRODUCTION_CONTROLLER_ANCESTRY_V1',path=controller['path'],sha256=controller['sha256'],delegated_helper=original_helper,auditor_mechanical_amendment=mechanical,production_controller_amendment=controller_pin,original_pilot_preapproval_unchanged=True)
    need(config['actual_production_controller']==expected,'Actual production configuration truthfully records outer controller ancestry')
    root=artifact(readiness['readiness_root_approval']);need(root['auditor_mechanical_amendment']==mechanical and root['production_controller_amendment']==controller_pin,'Root readiness mirrors both explicit postpilot amendments')
    for pin in (mechanical,controller_pin):
        matches=[row for row in config['input_files']if hashlib.sha256(row['content'].encode('utf-8')).hexdigest()==pin['sha256']];need(len(matches)==1 and matches[0]['content'].encode('utf-8')==(ROOT/safe_relative(pin['path'])).read_bytes(),'Both exact root amendment JSON bytes staged as original production aliases')
    audit.check('separately approved production controller delegates original one-shot workflow with explicit ancestry',True)
    return expected


def admission_closure(audit,a,cfg,freeze):
    """Independent original archived policy/budget/one-shot transaction proof."""
    tx=a.admission_dir;phase=a.attempt_mode;folder=tx/phase
    need(not (tx/'outcome_unknown_or_partial.json').exists(),'Unknown or partial submission is reconcile-only, never accepted')
    root=audit.read(a.admission_root_approval);need(re.fullmatch('[0-9a-f]{64}',a.admission_root_approval_sha256) and sha(a.admission_root_approval)==a.admission_root_approval_sha256,'Explicit external root admission approval SHA')
    need(root['schema']=='ISAAC_X16_TWO_PHASE_ROOT_APPROVAL_V1' and root['status']=='ROOT_APPROVED_FOR_EXACT_TWO_PHASE_WORKFLOW' and root['release']==RELEASE and root['source_freeze']['sha256']==a.freeze_sha256 and root['initial_registry_sha256']=='580cbda04dca1f6d995e3dde70cf38d123fcef79006a02b0f4298375c99aea08' and tuple(root['initial_job_ids'])==PRIOR_JOBS and root['wall_minutes']==[5,20] and root['cpus']==1 and root['memory_mb']==3500,'Root-approved exact two-phase source/predecessor/resource scope')
    need(root['source_peer_review']['sha256']=='6b2e80125097a45117e7df47a8833a6e9d3ff542e8722e0528ca2df208446b10' and root['source_archive']['sha256']==freeze['archive_sha256'],'Root approval binds exact independent science source and frozen archive')
    need(all(root[k] is True for k in ('reserve_both_caps_before_pilot','one_sbatch_per_phase','separate_production_readiness_required','no_retry_after_intent')) and root['pilot_suffix_seconds']==1500 and root['production_suffix_seconds']==1200 and root['science_cutoff_unix']==CUTOFF,'Root approved full suffix, reservation and no-retry guards')
    intents={};archives={};documents={}
    phases=['pilot'] if phase=='pilot' else ['pilot','production']
    for label in phases:
        sub=tx/label;intent=audit.read(sub/'intent.json');need(intent['schema']=='ISAAC_X16_DURABLE_ATTEMPT_INTENT_V1' and intent['phase']==label and intent['maximum_sbatch_calls']==1 and intent['any_existing_attempt_blocks_retry'] is True and abs(intent['full_reserved_cpu_hours']-25/60)<1e-12,'One durable exclusive original phase intent')
        need(set(intent['archived_inputs'])==set(intent['inputs']),'Every qualified original input has exact immutable archive')
        found={};aliases=set()
        for original,row in intent['archived_inputs'].items():
            rel=str(safe_relative(row['archive']));need(rel.startswith('inputs/') and rel.casefold() not in aliases,'Unique archive role inside phase inputs');aliases.add(rel.casefold());p=audit.file(sub/rel)
            need(p.stat().st_size==row['bytes'] and sha(p)==row['sha256']==intent['inputs'][original],'Exact original archived input bytes '+original);found[original]=p
        intents[label]=intent;archives[label]=found
        complete=audit.read(sub/'completed.json');transport=audit.read(sub/'transport.json');receipt=audit.read(sub/'remote_receipt.json');remote=audit.read(sub/'remote_configuration.json')
        need(complete['schema']=='ISAAC_X16_ATTEMPT_ACKNOWLEDGED_V1' and complete['phase']==label and complete['status']=='ACKNOWLEDGED' and complete['one_sbatch_acknowledged'] is True and complete['retry_permitted'] is False and complete['automatic_production'] is False and complete['reserved_production_reused'] is (label=='production') and abs(complete['full_reserved_cpu_hours']-25/60)<1e-12,'Unambiguous original acknowledgement and retained full reservation')
        need(transport['returncode']==0 and json.loads(transport['stdout'],object_pairs_hook=json_object)==receipt and receipt['returncode']==0 and str(receipt['job_id'])==str(complete['job_id']),'Original transport stdout and unique phase receipt agree')
        policy_hash=sha(tx/'reserve_both_policy.json') if label=='pilot' else intent['inputs'][str((ROOT/'overnight_20261022/policy.json').resolve())]
        expected=dict(intent['configuration'],local_attempt_intent_sha256=sha(sub/'intent.json'),current_registry_sha256=sha(tx/('reserve_both_registry.json' if label=='pilot' else 'pilot_acknowledged_registry.json')),current_policy_sha256=policy_hash,fully_reserved_cpu_hours=25/60)
        need(remote==expected,'Original remote configuration is exactly this archived exclusive intent')
        need(complete['current_registry_sha256']==sha(tx/(label+'_acknowledged_registry.json')) and complete['current_policy_sha256']==sha(tx/(label+'_acknowledged_policy.json')),'Original phase completion binds exact paired acknowledgement')
        documents[label]=dict(complete=complete,receipt=receipt,configuration=remote)
    original=(tx/'registry_before.json').read_bytes();audit.file(tx/'registry_before.json');audit.file(tx/'policy_before.json')
    need(hashlib.sha256(original).hexdigest()==root['initial_registry_sha256'],'Exact initial twelve registry bytes')
    original_rows=registry_object_fragments(original);need(tuple(str(r[0]['job_id']) for r in original_rows)==PRIOR_JOBS,'Exact original twelve registry IDs/order')
    previous_reg,previous_pol=original,(tx/'policy_before.json').read_bytes()
    for label in ['reserve_both']+['pilot_acknowledged']+(['production_acknowledged'] if phase=='production' else []):
        regp=tx/(label+'_registry.json');polp=tx/(label+'_policy.json');reg=audit.read(regp);pol=audit.read(polp);intent=audit.read(tx/(label+'_commit_intent.json'));complete=audit.read(tx/(label+'_commit_complete.json'));regbytes=regp.read_bytes();polbytes=polp.read_bytes()
        if label=='production_acknowledged':
            archived_policy=archives['production'][str((ROOT/'overnight_20261022/policy.json').resolve())];previous_pol=archived_policy.read_bytes()
            baseline=json.loads((tx/'pilot_acknowledged_policy.json').read_bytes());current=json.loads(previous_pol);current.pop('x16_production_readiness',None);need(current==baseline,'Only explicit production readiness added between original paired phases')
        need(intent['before_registry_sha256']==hashlib.sha256(previous_reg).hexdigest() and intent['before_policy_sha256']==hashlib.sha256(previous_pol).hexdigest() and intent['after_registry_sha256']==sha(regp)==complete['registry_sha256'] and intent['after_policy_sha256']==sha(polp)==complete['policy_sha256'],'Durable paired transaction before/after hash chain')
        fragments=registry_object_fragments(regbytes);need(len(fragments)==14 and fragments[:12]==original_rows,'Every original twelve object and literal UTF8/CRLF fragment preserved')
        ours=reg['jobs'][12:];need([r['reservation_id'] for r in ours]==['x16_v1_pilot','x16_v1_production'] and [r['wall_minutes'] for r in ours]==[5,20] and all(r['catalogue_id']=='X16' and r['cpus']==1 and r['memory_mb']==3500 and r['freeze_sha256']==a.freeze_sha256 and r['scientific_compute'] is True and r['scientific_sampling'] is False for r in ours) and abs(sum(r['reserved_cpu_hours'] for r in ours)-25/60)<1e-12,'Both original caps reserved once before any science submission')
        if label=='reserve_both':need(all(r['job_id'] is None for r in ours) and ours[0]['state'] in ('RESERVED_NOT_SUBMITTED','PREPARED_UNKNOWN_UNTIL_RECONCILED') and ours[1]['state']=='RESERVED_NOT_SUBMITTED','Both phases fully reserved before first remote intent; pilot attempt cannot be silently retried')
        else:acknowledged_transition(previous_reg,regbytes,label.removesuffix('_acknowledged'),documents[label.removesuffix('_acknowledged')])
        changed=json.loads(previous_pol);changed['submitted_registry_sha256']=sha(regp);need(pol==changed,'Only policy registry digest changed in paired commit')
        previous_reg,previous_pol=regbytes,polbytes
    need(str(documents[phase]['complete']['job_id'])==a.job_id and str(json.loads(previous_reg)['jobs'][12+(phase=='production')]['job_id'])==a.job_id,'Exact canonical registered phase ID, no replacement')
    if phase=='production':need(str(json.loads(previous_reg)['jobs'][12]['job_id'])==str(documents['pilot']['complete']['job_id']),'Final production registry retains exact independently acknowledged pilot ID')
    for label in phases:
        found=archives[label];intent=intents[label];config=documents[label]['configuration'];receipt=documents[label]['receipt'];sub=tx/label
        def archived(path):
            original=str(Path(path).resolve());need(original in found,'Required pre-admission original input archive '+original);return audit.read(found[original])
        policy=archived(ROOT/'overnight_20261022/policy.json');registry=archived(ROOT/'overnight_20261022/submissions.json');resource=archived(ROOT/'overnight_20261022/resource_ledger.json');checkpoint=archived(ROOT/'latest_checkpoint.json');approval=policy['x16_admission']
        mirror=('release','initial_registry_sha256','initial_job_ids','source_freeze','source_peer_review','source_archive','staging_receipt','source_preflight','submission_helper_sha256','submission_peer_review','submission_validation','submission_tests','selected_base_ledger','selected_checkpoint','stage3_terminal_closure','required_consumer_helper_pins','independent_auditor','wall_minutes','cpus','memory_mb')
        need(approval['schema']=='ISAAC_X16_TWO_PHASE_ADMISSION_V1' and approval['status']=='ROOT_APPROVED' and all(approval[k]==root[k] for k in mirror) and approval['root_approval']['sha256']==a.admission_root_approval_sha256,'Archived policy mirrors exact external root source/admission approval')
        root_original=ROOT/safe_relative(approval['root_approval']['path']);need(str(root_original.resolve()) in found and sha(found[str(root_original.resolve())])==a.admission_root_approval_sha256,'Trusted root approval also retained inside original phase archives')
        def artifact(row):
            original=str((ROOT/safe_relative(row['path'])).resolve());need(original in found and sha(found[original])==row['sha256'],'Every root-approved original artifact retained exactly '+row['path']);p=found[original];return audit.read(p) if Path(row['path']).suffix=='.json' else p
        helper_roles=[name for name,want in policy['reviewed_local_helper_sha256'].items() if want==root['submission_helper_sha256'] and name in ('overnight_20261022/submit_x16_after_queue_candidate.py','overnight_20261022/submit_x16_after_queue_candidate_v2.py','overnight_20261022/submit_x16_after_queue_candidate_v3.py')];need(len(helper_roles)==1,'One exactly approved two-phase admission source role');helper_role=helper_roles[0];helper_original=str((ROOT/helper_role).resolve());need(helper_original in found and sha(found[helper_original])==root['submission_helper_sha256'],'Exact root-reviewed admission source retained before execution')
        for name,want in [('submit_one.py','77352ab312dbaaccb294e889db734e34a8ae86af455b3603e3948669af4d8bb0'),('submit_x13_after_queue.py','209c94bda9c90850d0d7205fd8efd23b8b6dcae4189af57eda82857a2d737660')]:
            original=str((ROOT/'overnight_20261022'/name).resolve());need(original in found and sha(found[original])==want,'Preserved independently reviewed transaction/single-sbatch primitive')
        submission_peer=artifact(root['submission_peer_review']);submission_validation=artifact(root['submission_validation']);artifact(root['submission_tests'])
        need(submission_peer['status']=='APPROVED_FOR_X16_ADMISSION_CANDIDATE_SOURCE_ONLY' and submission_peer['submission_helper_sha256']==root['submission_helper_sha256'] and submission_peer['test_sha256']==root['submission_tests']['sha256'] and submission_validation['status']=='PASS' and submission_validation['files'][Path(helper_role).name]==root['submission_helper_sha256'],'Independent source/fixture approval of exact one-shot admission helper')
        need(set(root['required_consumer_helper_pins'])=={'cancellation','checkpoint_refresh','retrieval','accounting'},'All four policy consumers explicitly source-reviewed')
        for role,row in root['required_consumer_helper_pins'].items():
            artifact(row);review=artifact(row['peer_review']);need(review['status']=='APPROVED_FOR_X16_CONSUMER_USE' and review['helper_sha256']==row['sha256'] and review['helper_path'] in (row['path'],ROOT.name+'/'+row['path']) and policy['reviewed_local_helper_sha256'][row['path']]==row['sha256'],'Exact reviewed consumer source/peer policy binding '+role)
        stage=artifact(root['staging_receipt']);preflight=artifact(root['source_preflight']);need(stage['all_remote_bytes_verified'] is True and stage['source_freeze_sha256']==a.freeze_sha256 and stage['remote_source']==REMOTE_ROOT+'/source/manual_pc_continuation/'+RELEASE and preflight['status']=='PASS' and preflight['source_freeze_sha256']==a.freeze_sha256 and preflight['native_or_model_execution'] is False and all(r['passed'] is True for r in preflight['checks']),'Root verified exact staged bytes and source-only runtime syntax/native preflight')
        tool=approval['independent_auditor'];need(tool==a._original_auditor_tool and tool['path']=='audit_x16_pair_binding_v2.py' and tool['sha256']==ORIGINAL_AUDITOR_SHA,'Original exact5c auditor remains the historically preapproved pilot verifier; replacement has separate postpilot authority')
        artifact(tool)
        for role in ('peer_review','validation','tests','root_approval'):
            original=str((ROOT/safe_relative(tool[role]['path'])).resolve());need(original in found and sha(found[original])==tool[role]['sha256'],'Original independently approved auditor artifact '+role)
        peer=archived(ROOT/tool['peer_review']['path']);validation=archived(ROOT/tool['validation']['path']);auditor_root=archived(ROOT/tool['root_approval']['path'])
        need(peer['status']=='APPROVED_FOR_LOCAL_RECEIPT_AUDITOR_SOURCE_ONLY' and peer['auditor_sha256']==tool['sha256'] and peer['test_sha256']==tool['tests']['sha256'] and validation['status']=='PASS' and validation['files'][tool['path']]==tool['sha256'] and validation['files'][Path(tool['tests']['path']).name]==tool['tests']['sha256'],'Independent pre-run source/fixture approval for original historical auditor')
        need(auditor_root['schema']=='ISAAC_X16_LOCAL_RECEIPT_AUDITOR_APPROVAL_V1' and auditor_root['status']=='ROOT_APPROVED_FOR_LOCAL_RECEIPT_AUDIT' and all(auditor_root[k]==tool[k] for k in ('path','sha256','peer_review','validation','tests')) and auditor_root['source_freeze_sha256']==a.freeze_sha256,'Separate root approval for exact independent raw auditor')
        need(policy['submitted_registry_sha256']==intent['inputs'][str((ROOT/'overnight_20261022/submissions.json').resolve())] and registry['jobs'][:12]==[r[0] for r in original_rows],'Archived canonical pair and original twelve identities')
        ledger=artifact(approval['selected_base_ledger']);catalogue=audit.read(ROOT/'benchmark_ideas_100.json');need(sha(ROOT/'benchmark_ideas_100.json')==CATALOGUE_SHA and approval['selected_base_ledger']['sha256']=='e44c10b49f987aa052772b9316228c04a93f7d477c2dd0f3f5188885085ead54' and checkpoint['ledger']==approval['selected_base_ledger']['path'] and checkpoint['ledger_sha256']==approval['selected_base_ledger']['sha256'] and checkpoint['counts']==ledger['counts'] and len(ledger['rows'])==100 and [r['id'] for r in ledger['rows']]==[r['id'] for r in catalogue],'Explicit root-selected exact archived current100-ID ledger and checkpoint')
        need(approval['selected_checkpoint']['path']=='latest_checkpoint.json','Root-selected canonical initial checkpoint role');checkpoint_roster(checkpoint,a.freeze_sha256,freeze['files']['SHA256SUMS'])
        initial_cp_path=audit.file(tx/'checkpoint_before.json');need(sha(initial_cp_path)==approval['selected_checkpoint']['sha256'],'Initial exact checkpoint preserved separately before reservations')
        initial_alias=archives['pilot'][str((ROOT/'latest_checkpoint.json').resolve())];need(initial_alias.read_bytes()==initial_cp_path.read_bytes(),'Original pilot checkpoint archival alias agrees with deterministic preservation')
        if label=='pilot':need(intent['inputs'][str((ROOT/'latest_checkpoint.json').resolve())]==approval['selected_checkpoint']['sha256'],'Exact root-selected initial checkpoint bytes at first admission')
        else:
            readiness_record=policy['x16_production_readiness'];readiness_root=artifact(readiness_record['readiness_root_approval']);pin=dict(path='latest_checkpoint.json',sha256=intent['inputs'][str((ROOT/'latest_checkpoint.json').resolve())]);need(readiness_record['current_checkpoint']==readiness_root['current_checkpoint']==pin,'Separate explicit root pin for any current production checkpoint refresh')
            initial_cp=load(initial_cp_path);need(initial_cp['ledger']==checkpoint['ledger'] and initial_cp['ledger_sha256']==checkpoint['ledger_sha256'] and initial_cp['counts']==checkpoint['counts'] and initial_cp['frozen_releases']==checkpoint['frozen_releases'] and initial_cp['x01_stage3_terminal_closure']==checkpoint['x01_stage3_terminal_closure']==approval['stage3_terminal_closure'],'Production checkpoint refresh changes no selected science, frozen roster or original incomplete-stage closure')
            production_controller_authority(audit,a,policy,config,found)
        closed=artifact(approval['stage3_terminal_closure']);need(closed['schema']=='ISAAC_X01_STAGE3_ROOT_TERMINAL_CLOSURE_V1' and closed['status']=='ROOT_CLOSED_EXECUTION_INCOMPLETE_WITHOUT_LEDGER_REVISION' and closed['registry_sha256']==root['initial_registry_sha256'] and closed['registered_job_ids']==list(PRIOR_JOBS) and closed['selected_ledger']==approval['selected_base_ledger'] and closed['all_twelve_registered_jobs_terminal'] is True and closed['new_catalogue_completions']==0 and closed['final_reserved_statistical_look_closed'] is True,'Original exact-hash incomplete stage3 remains execution history, no synthetic scientific ledger adoption')
        need(all(closed[k] is False for k in ('scientific_ledger_changed','complete_scientific_audit_PASS','ledger_updater_executed','comparator_replay_executed','primary_mathematical_reconstruction_executed','fourth_look_authorized','retry_extension_pooling_seed_gate_changes_authorized','X16_execution_authorized_by_this_record')),'No fourth look or partial promotion licensed by closure')
        snapshot_path=resource['inputs']['snapshot']['path'];need(snapshot_path in found,'Original fresh visibility snapshot archived');snapshot=audit.read(found[snapshot_path]);created=dt.datetime.fromisoformat(intent['created_utc']).timestamp();finished=dt.datetime.fromisoformat(snapshot['finished_utc']).timestamp()
        need(snapshot['status']=='VISIBLE' and snapshot['visibility_loss'] is False and not snapshot.get('registry_changed_during_query') and snapshot['registry_sha256']==intent['inputs'][str((ROOT/'overnight_20261022/submissions.json').resolve())] and snapshot['policy_sha256']==intent['inputs'][str((ROOT/'overnight_20261022/policy.json').resolve())] and -5<=created-finished<=120,'Original fresh exact policy/registry-bound whole-user visibility')
        for role,original in [('registry',ROOT/'overnight_20261022/submissions.json'),('policy',ROOT/'overnight_20261022/policy.json'),('snapshot',Path(snapshot_path))]:need(resource['inputs'][role]['sha256']==intent['inputs'][str(original.resolve())],'Original accounting input hash binding '+role)
        need(resource['fresh_visibility_for_releases'] is True and resource['budget_and_visibility_within_policy'] is True and resource['accounting_issues']==[] and resource['exceeds_cap'] is False and resource['new_overnight_cpu_hour_cap']==policy['agent_chosen_new_overnight_cpu_hour_cap']==10,'Original fresh ten-hour accounting and visibility qualification')
        known=PRIOR_JOBS+((str(documents['pilot']['complete']['job_id']),) if label=='production' else ())
        actual=0.
        for job in known:
            rows=[r for r in snapshot['remote']['sacct_rows'] if r['JobIDRaw']==job];need(len(rows)==1 and rows[0]['User']=='ywang287' and rows[0]['State'].split()[0].rstrip('+') in TERMINAL,'Original required predecessor terminal/owned');row=rows[0];hours=int(row['ElapsedRaw'])*int(row['AllocCPUS'])/3600;actual+=hours;accounted=[r for r in resource['jobs'] if str(r.get('job_id'))==job];need(len(accounted)==1 and accounted[0]['accounting_state']=='VERIFIED_TERMINAL' and accounted[0]['reserved_cpu_hours']==0 and abs(accounted[0]['verified_actual_cpu_hours']-hours)<1e-12,'Original main allocation fully reconciled once')
        governed=resource['governed_actual_plus_open_reservations'];need(governed>=actual-1e-12 and governed+(25/60 if label=='pilot' else 0)<=10+1e-12,'Complete originally reserved workflow fits governed budget')
        if label=='pilot':need(resource['reserved_new_cpu_hours']==0,'Both caps first reserved after zero previous open reservations')
        else:need(any(r.get('reservation_id')=='x16_v1_production' and r['accounting_state']=='OPEN_FULL_RESERVATION' and r['reserved_cpu_hours']>=20/60 for r in resource['jobs']),'Production reuses existing full20-minute reservation')
        expected_job=str(documents[label]['complete']['job_id']);pilot=str(documents['pilot']['complete']['job_id']) if label=='production' else None
        need(config['phase']==label and config['key']=='x16_v1_'+label and config['release']==RELEASE and config['wrapper']==label+'.sbatch' and config['wall_minutes']==(5 if label=='pilot' else 20) and config['full_suffix_seconds']==(1500 if label=='pilot' else 1200) and config['freeze_sha256']==a.freeze_sha256 and config['dependencies']==([pilot] if pilot else []) and config['dependency_kind']=='afterok' and config['pilot_job_id']==pilot and config['environment']['X16_EXPECTED_FREEZE_SHA256']==a.freeze_sha256,'Exact originally frozen phase/dependency/full-suffix configuration')
        need(config['single_remote_sha256']=='0c53c774dc78985e9765c8400aa3bd6c5431fe1dd0afa34603e06d27f3867257' and config['source']==REMOTE_ROOT+'/source/manual_pc_continuation/'+RELEASE and config['result_prefix']==RELEASE+'_'+label+'_' and receipt['result']==REMOTE_ROOT+'/results/'+config['result_prefix']+expected_job and receipt['source']==config['source'] and receipt['submission_dir']==REMOTE_ROOT+'/submissions/overnight_20261023/x16_v1_'+label,'Exact source/result/submission namespace and preserved original single-sbatch code')
        if label==phase:remote_submission_closure(audit,a,config,receipt,known)
    audit.check('complete original two-phase admission/current-ledger/root/source/auditor/budget/paired-transaction proof',True)
    return dict(schema='X16_INDEPENDENT_TWO_PHASE_ADMISSION_PROOF_V1',phase=phase,job_id=a.job_id,root_approval_sha256=a.admission_root_approval_sha256,phase_intent_sha256=sha(folder/'intent.json'),initial_registry_sha256=root['initial_registry_sha256'],acknowledged_registry_sha256=sha(tx/(phase+'_acknowledged_registry.json')),reserved_total_cpu_hours=25/60,original_twelve_literal_rows_preserved=True,no_retry=True,production_separate_readiness_required=True)


def production_pilot_closure(audit,a,cfg,freeze,execution,record):
    """Repeat the full independent pilot mathematics; never pool its values."""
    qualification=execution['actual_start_admission']['pilot_qualification'];need(qualification['schema']=='X16_PRODUCTION_PILOT_GUARD_V1' and qualification['status']=='PASS' and qualification['source_freeze_sha256']==a.freeze_sha256 and qualification['read_only_qualification'] is True and qualification['production_physics_pooling'] is False,'Production preserves independent pilot qualification as metadata only')
    pilot_id=qualification['job_id'];need(re.fullmatch('[0-9]+',pilot_id) and pilot_id!=a.job_id and pilot_id not in PRIOR_JOBS,'One exact distinct excluded pilot identity')
    tx=a.admission_dir;intent=audit.read(tx/'production/intent.json');aliases=intent['archived_inputs'];original_policy=str((ROOT/'overnight_20261022/policy.json').resolve());policy=audit.read(tx/'production'/aliases[original_policy]['archive']);selected=policy['x16_production_readiness']
    need(selected['schema']=='ISAAC_X16_PRODUCTION_READINESS_ADMISSION_V1' and selected['status']=='ROOT_APPROVED_FOR_ONE_FIXED_PRODUCTION' and selected['pilot_job_id']==pilot_id,'Separate root-approved production readiness admission')
    by_remote={}
    for row in record['files']:
        need(row['remote'] not in by_remote,'Unique original production transferred remote path');by_remote[row['remote']]=a.job_dir/row['local']
    readiness_path=qualification['readiness_path'];need(readiness_path in by_remote,'Original staged root readiness alias retained')
    readiness=audit.read(by_remote[readiness_path]);need(sha(by_remote[readiness_path])==qualification['readiness_sha256']==selected['readiness_root_approval']['sha256'],'Exact staged root approval bytes equal original approved hash')
    if a.pilot_readiness is not None:need(a.pilot_readiness==by_remote[readiness_path].resolve(),'Explicit pilot readiness path matches original staged alias')
    need(readiness['schema']=='X16_PILOT_READINESS_APPROVAL_V1' and readiness['status']=='ROOT_APPROVED_FOR_FIXED_PRODUCTION' and readiness['job_id']==pilot_id and readiness['release']==RELEASE and readiness['source_freeze_sha256']==a.freeze_sha256 and readiness['pilot_wall_cap_seconds']==300 and readiness['scientific_completions']==0,'Original root-qualified excluded300-second pilot and identical source')
    roles={'independent_audit','auditor_source','auditor_peer_review','auditor_validation','transfer','scheduler'}
    need(set(readiness['artifacts'])==set(selected['artifacts'])==set(qualification['artifacts'])==roles and readiness['artifacts']==qualification['artifacts'],'Complete six immutable pilot qualification aliases')
    artifacts={}
    for role,row in readiness['artifacts'].items():
        need(row['path'] in by_remote and row['sha256']==selected['artifacts'][role]['sha256'],'Original remote alias corresponds to root-approved local artifact');p=audit.file(by_remote[row['path']]);need(sha(p)==row['sha256'],'Original staged pilot qualification artifact bytes '+role);artifacts[role]=p
        original=str((ROOT/safe_relative(selected['artifacts'][role]['path'])).resolve());need(original in aliases and aliases[original]['sha256']==row['sha256'] and sha(tx/'production'/aliases[original]['archive'])==row['sha256'],'Exact original local qualification bytes permanently archived before production')
    need(readiness['artifacts']['auditor_source']['sha256']==sha(__file__),'Pilot used this exact independently reviewed auditor source')
    old=audit.read(artifacts['independent_audit']);need(old['schema']==SCHEMA and old['integrity_status']=='PASS' and old['attempt_mode']=='pilot' and old['job_id']==pilot_id and old['source_freeze_sha256']==a.freeze_sha256 and old['audit_tool_sha256']==sha(__file__) and old['scientific_sampling'] is False and old['remote_calls'] is False and old['failed_checks']==0 and old['check_count']==len(old['checks'])>0 and all(row['passed'] is True for row in old['checks']),'Complete original independent pilot integrity report with approved source')
    need(old['postpilot_mechanical_auditor_amendment']==a._mechanical_authority_proof,'Repaired pilot audit truthfully carries this same separately approved postpilot authority')
    ready=old['pilot_readiness'];need(ready['status']=='READY' and all(ready[k] is True for k in ('numerical_qualification','full_grid_cross_backend_match','timing_precision_ready')),'Original excluded pilot independently numerically and operationally READY')
    outputs=readiness['pilot_output_sha256'];need(set(outputs)==OUTPUT_ROLES and outputs==old['pilot_output_sha256']==selected['pilot_output_sha256']==qualification['pilot_output_sha256'],'All21 complete immutable pilot outputs bound by every readiness layer')
    pilot_dir=(ROOT/safe_relative(selected['artifacts']['transfer']['path'])).resolve().parent
    need((ROOT/safe_relative(selected['pilot_job_dir'])).resolve()==pilot_dir/'results','Selected pilot output directory corresponds to its original transfer parent')
    if a.pilot_job_dir is not None:need(a.pilot_job_dir==pilot_dir,'Explicit original pilot transfer parent differs')
    original_audit=(ROOT/safe_relative(selected['artifacts']['independent_audit']['path'])).resolve()
    if a.pilot_audit is not None:need(a.pilot_audit==original_audit,'Explicit original pilot audit path differs')
    need(sha(audit.file(original_audit))==readiness['artifacts']['independent_audit']['sha256'],'Retain original pilot audit as well as staged exact alias')
    for name,want in outputs.items():need(sha(audit.file(pilot_dir/'results'/safe_relative(name)))==want and want in old['input_sha256'].values(),'Original complete excluded pilot output unchanged and fully audited')
    for name,want in old['input_sha256'].items():need(sha(audit.file(name))==want,'Original pilot audit input closure remains unchanged')
    for role in ('transfer','scheduler'):need(sha(audit.file(pilot_dir/('transfer.json' if role=='transfer' else 'scheduler.txt')))==readiness['artifacts'][role]['sha256'],'Original pilot scheduler/transfer match staged qualification')
    live=qualification['live_terminal_row'];need(live['JobIDRaw']==pilot_id and live['JobName']=='x16_pilot' and live['State']=='COMPLETED' and live['User']=='ywang287' and live['AllocCPUS']=='1' and live['ExitCode']=='0:0' and scheduler_epoch(live['End'])<=execution['actual_batch_start_unix']+1,'Fresh owned pilot COMPLETED0 before actual production')
    need(all(row['job_id']!=pilot_id for row in execution['actual_start_admission']['active_queue_guard']['rows']),'Expected pilot absent under every active JobName')
    pilot_args=argparse.Namespace(**vars(a));pilot_args.job_dir=pilot_dir;pilot_args.job_id=pilot_id;pilot_args.attempt_mode='pilot';pilot_args.pilot_job_dir=pilot_args.pilot_audit=pilot_args.pilot_readiness=None;pilot_args.replay_frozen=False
    child=Audit();child.file(__file__);repeated=reconstruct_job(child,pilot_args,cfg,freeze)
    for row in child.checks:audit.checks.append(dict(row,name='pilot_reaudit/'+row['name']))
    for name,want in child.inputs.items():
        need(name not in audit.inputs or audit.inputs[name]==want,'Pilot and production bind same persistent bytes');audit.inputs[name]=want
    need(child.integrity() and repeated['pilot_readiness']['status']=='READY' and all(repeated['pilot_readiness'][k] is True for k in ('numerical_qualification','full_grid_cross_backend_match','timing_precision_ready')),'Full independent repeat of pilot raw mathematics must qualify READY; stored READY alone is insufficient')
    audit.check('production independent repeat of complete excluded pilot source/process/raw mathematics',True)
    return dict(schema='X16_INDEPENDENT_PILOT_REAUDIT_PROOF_V1',job_id=pilot_id,original_audit_sha256=sha(original_audit),readiness_root_approval_sha256=sha(by_remote[readiness_path]),pilot_output_sha256=outputs,independent_scientific_status=repeated['scientific_status'],pilot_readiness=repeated['pilot_readiness'],primary_raw_reconstruction_repeated=True,physics_samples_pooled=False,child_integrity_checks=len(child.checks),child_numerical_checks=len(child.numerical),child_cross_checks=len(child.comparisons))


def reconstruct_job(audit,a,cfg,freeze):
    record,scheduler_rows=transfer_closure(audit,a,freeze)
    admission=admission_closure(audit,a,cfg,freeze)
    execution,outputs,resources=preliminary_execution(audit,a,cfg,freeze,scheduler_rows)
    if a.attempt_mode=='production':pilot=production_pilot_closure(audit,a,cfg,freeze,execution,record)
    else:
        need(execution['actual_start_admission']['pilot_qualification']==dict(status='EXCLUDED_PILOT_NO_PRIOR_PRODUCTION_INPUT'),'Excluded pilot consumes no production or prior pilot numerical values')
        pilot=None
    need(audit.integrity(),'Complete original source/admission/transfer/process closure required before independent raw mathematics')
    result=a.job_dir/'results';target=audit.read(result/'target.toml');reference=audit.read(result/'reference.json');saved=audit.read(result/'comparison.json')
    expected={REMOTE_ROOT+'/source/manual_pc_continuation/'+RELEASE+'/manifest.toml':sha(a.release_dir/'manifest.toml'),record['result_root']+'/target.toml':sha(result/'target.toml'),record['result_root']+'/reference.json':sha(result/'reference.json')}
    need(saved['input_sha256']==expected,'Original frozen comparator binds exact original manifest and both raw receipts')
    t=audit_backend(audit,target,'target',cfg,a.freeze_sha256,a.release_dir,execution);r=audit_backend(audit,reference,'reference',cfg,a.freeze_sha256,a.release_dir,execution)
    if all(row['passed'] for row in audit.numerical):bindings,grounds=cross_backends(audit,t,r)
    else:
        bindings=[dict(case_id=x['case_id'],reference_classification=b['classification'],**aa) for x,y in zip(t,r) for aa,b in zip(x['bindings'],y['bindings'])];grounds=[]
    decision,ready=scientific_status(audit,bindings)
    diagnostics=saved_comparison(audit,saved,a.attempt_mode,decision,ready,t,bindings,grounds)
    comparator_exit={'PASS':0,'FAIL':1,'UNRESOLVED':2}[decision];audit.check('independently qualified comparator exact exit',execution['process_exit_codes']['comparator']==comparator_exit)
    wrapper_exit=(0 if ready['status']=='READY' and comparator_exit in (0,2) else 2) if a.attempt_mode=='pilot' else comparator_exit
    if a.attempt_mode=='pilot':audit.check('pilot process readiness exit distinct from physical sign/comparator exit',execution['pilot_readiness_exit_code']==wrapper_exit)
    audit.check('final scheduler main/batch status matches exact actual wrapper exit',all(row['State']==('COMPLETED' if wrapper_exit==0 else 'FAILED') and row['ExitCode']==str(wrapper_exit)+':0' for row in scheduler_rows))
    decision,ready=scientific_status(audit,bindings);ready['timing_precision_ready']=ready['status']=='READY' and audit.integrity()
    return dict(scientific_status=decision,pilot_readiness=ready,bindings=bindings,ground_receipts=grounds,backends=dict(target=t,reference=r),original_output_sha256=outputs,execution=execution,resources=resources,admission=admission,production_pilot=pilot,comparison_diagnostics=diagnostics)


def compact_backend(cases):
    result=[]
    for c in cases:
        grounds=[]
        for g in c['grounds']:
            gamma=g['pair_density'];grounds.append(dict(N=g['N'],dimension=g['dimension'],energy_mu=g['energy_mu'],energy_physical=g['energy_physical'],energy_lower=g['energy_lower'],energy_upper=g['energy_upper'],pair_density=dict(rows=28,columns=28,real=gamma.real.ravel().tolist(),imag=gamma.imag.ravel().tolist()),diagnostic_lower_range=g['diagnostic_lower_range'],diagnostic_upper_range=g['diagnostic_upper_range']))
        result.append(dict(case_id=c['case_id'],all_eigenvalues=c['all_eigenvalues'],grounds=grounds,bindings=c['bindings']))
    return result


def plain(value):
    if isinstance(value,np.ndarray):return plain(value.tolist())
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {str(k):plain(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [plain(v) for v in value]
    return value


def replay_frozen(audit,a,saved):
    """Optional supplemental deterministic receipt comparator, no native solve."""
    need(audit.integrity(),'No candidate code before independent full integrity closure')
    need(all(sha(p)==want for p,want in audit.inputs.items()),'All original bound inputs unchanged immediately before optional replay')
    import subprocess
    with tempfile.TemporaryDirectory(prefix='x16_receipt_replay_',dir=a.output.parent) as d:
        output=Path(d)/'comparison.json';env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',X16_ATTEMPT_MODE=a.attempt_mode,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        command=[sys.executable,'-B',str(a.release_dir/'compare.py'),'--manifest',str(a.release_dir/'manifest.toml'),'--freeze-sha256',a.freeze_sha256,'--target',str(a.job_dir/'results/target.toml'),'--reference',str(a.job_dir/'results/reference.json'),'--output',str(output)]
        completed=subprocess.run(command,cwd=a.release_dir,env=env,capture_output=True,timeout=120,check=False);need(output.is_file(),'Frozen receipt replay must produce complete comparison')
        replay=load(output);keys=('scientific_status','attempt_mode','catalogue_status','pilot_readiness','complete_cases','spinblocks_per_backend','fixedN_ground_ensembles','pair_modes','complex_pairdensity_entries')
        audit.check('supplemental frozen comparator verdict/readiness/grid exact',all(replay.get(k)==saved.get(k) for k in keys) and completed.returncode=={'PASS':0,'FAIL':1,'UNRESOLVED':2}[saved['scientific_status']])
        return dict(executed=True,exit_code=completed.returncode,saved_sha256=sha(a.job_dir/'results/comparison.json'),replayed_sha256=sha(output),decision_and_readiness_exact=all(replay.get(k)==saved.get(k) for k in keys),byte_identical=False,numeric_summary_policy='Primary raw formula checks and explicit conditioned diagnostics are authoritative; this optional replay compares exact decisions and scope.')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('job-dir','release-dir','freeze','source-review','source-archive','native-source-dir','catalogue','admission-dir','admission-root-approval','mechanical-amendment','output'):parser.add_argument('--'+name,type=Path,required=True)
    for name in ('job-id','freeze-sha256','source-review-sha256','admission-root-approval-sha256','mechanical-amendment-sha256'):parser.add_argument('--'+name,required=True)
    parser.add_argument('--production-controller-amendment',type=Path);parser.add_argument('--production-controller-amendment-sha256')
    parser.add_argument('--attempt-mode',choices=('pilot','production'),required=True);parser.add_argument('--pilot-job-dir',type=Path);parser.add_argument('--pilot-audit',type=Path);parser.add_argument('--pilot-readiness',type=Path);parser.add_argument('--replay-frozen',action='store_true');a=parser.parse_args(argv)
    for k,v in vars(a).items():
        if isinstance(v,Path):setattr(a,k,v.resolve())
    need(re.fullmatch('[0-9]+',a.job_id) and not a.output.exists(),'Numeric original job ID and fresh output required')
    protected=[a.release_dir,a.native_source_dir,a.admission_dir]+[a.job_dir/role for role in ('results','source','native_source','submissions')]
    if a.pilot_job_dir:protected.extend(a.pilot_job_dir/role for role in ('results','source','native_source','submissions'))
    need(all(not a.output.is_relative_to(p) for p in protected),'Audit output outside all preserved source/raw/admission trees')
    a.output.parent.mkdir(parents=True,exist_ok=True);audit=Audit();audit.file(__file__)
    report=dict(schema=SCHEMA,created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),audit_tool_sha256=sha(__file__),auditor_sha256=sha(__file__),job_id=a.job_id,attempt_mode=a.attempt_mode,catalogue_id='X16',source_freeze_sha256=a.freeze_sha256,scientific_sampling=False,scientific_completions=0,remote_calls=False,scientific_status='UNRESOLVED',pilot_excluded_from_production=a.attempt_mode=='pilot',catalogue_status='PILOT_EXCLUDED' if a.attempt_mode=='pilot' else 'PENDING_SEPARATE_APPEND_ONLY_LEDGER_REVIEW')
    try:
        cfg,freeze=source_closure(audit,a);report['postpilot_mechanical_auditor_amendment']=mechanical_authority(audit,a);result=reconstruct_job(audit,a,cfg,freeze)
        report.update(scientific_status=result['scientific_status'],pilot_readiness=result['pilot_readiness'],independent_reconstruction=dict(bindings=result['bindings'],ground_receipts=result['ground_receipts'],backends={key:compact_backend(value) for key,value in result['backends'].items()}),original_output_sha256=result['original_output_sha256'],resources=result['resources'],admission=result['admission'],production_pilot=result['production_pilot'],comparison_diagnostics=result['comparison_diagnostics'])
        if a.attempt_mode=='pilot':report['pilot_output_sha256']=result['original_output_sha256']
        elif result['production_pilot']:report['pilot_output_sha256']=result['production_pilot']['pilot_output_sha256']
        report['required_grid_counts']=dict(cases=4,spinblocks_per_backend=76,eigenvalues_per_backend=952,fixedN_ground_ensembles_per_backend=20,site_operators=8,pair_modes=28,pair_sources_per_backend=560,complex_pairdensity_entries_per_backend=15680,real_pairdensity_components_per_backend=31360,binding_intervals_per_backend=8)
        report['frozen_comparator_replay']=replay_frozen(audit,a,audit.read(a.job_dir/'results/comparison.json')) if a.replay_frozen else dict(executed=False,reason='Optional supplemental replay not requested; primary independent complete raw reconstruction remains authoritative')
    except Exception as error:audit.check('complete original X16 independent receipt audit',False,type(error).__name__+': '+str(error))
    for name,want in list(audit.inputs.items()):
        try:audit.check('original bound input unchanged at audit completion '+name,sha(name)==want)
        except Exception as error:audit.check('original bound input unchanged at audit completion '+name,False,type(error).__name__+': '+str(error))
    failed=[r for r in audit.checks if not r['passed']];report.update(integrity_status='PASS' if not failed else 'UNRESOLVED',audit_status='PASS' if not failed else 'UNRESOLVED',check_count=len(audit.checks),failed_checks=len(failed),failure_count=len(failed),failures=failed,checks=audit.checks,numerical_qualification_checks=audit.numerical,cross_backend_checks=audit.comparisons,input_sha256=audit.inputs)
    if failed:
        report['scientific_status']='UNRESOLVED';report['pilot_readiness']=dict(status='NOT_READY',numerical_qualification=False,full_grid_cross_backend_match=False,timing_precision_ready=False,reason='Incomplete or failed integrity cannot qualify production or a catalogue result')
    report['limitations']=['Conditional binary64 and BLAS/LAPACK numerical budgets with explicit64epsilon assumptions, not directed-rounding interval arithmetic or formal proof.','Finite four-case4site global-ground question through unchanged nativeH and NEW pair observer; no thermodynamic or preexisting native pair-solver claim.','All imaginary/zero/fullmatrix components and all spinblock/ground members retained; no spatial/spin filtering or averaging.','Pilot numerical readiness can exempt only physical sign uncertainty; no pilot physics is pooled into production or counted.','Original receipt/admission paths may be absolute Windows/ISAAC paths; relocation needs a separately reviewed procedure without rewriting original evidence.']
    raw=(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False)+'\n').encode('utf-8')
    with a.output.open('xb') as f:f.write(raw)
    print(json.dumps(dict(integrity_status=report['integrity_status'],scientific_status=report['scientific_status'],attempt_mode=a.attempt_mode,checks=report['check_count'],failed_checks=len(failed),output=str(a.output)),indent=2))
    return 0 if not failed else 2


if __name__=='__main__':raise SystemExit(main())


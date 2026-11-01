#!/usr/bin/env python3
"""Independent TeNPy finite-MPS oracle; run only in the allocated pilot/job."""
import argparse
import hashlib
import importlib.metadata
import itertools
import json
import os
from pathlib import Path
import time
import tomllib
import traceback

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sector_ed(c):
    """Small, independent sparse ED for energy/gap gates, not a new catalogue idea."""
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import eigsh
    L, nu, nd = c['L'], c['Nup'], c['Ndown']
    masks = []
    for up in itertools.combinations(range(L), nu):
        for dn in itertools.combinations(range(L), nd):
            masks.append(sum(1 << (2*i) for i in up) + sum(1 << (2*i+1) for i in dn))
    assert len(masks) <= 10000, 'Frozen small-sector resource guard'
    index = {m: i for i, m in enumerate(masks)}
    rows, cols, values = [], [], []
    for k, m in enumerate(masks):
        double = sum(((m >> (2*i)) & 3) == 3 for i in range(L))
        rows.append(k); cols.append(k); values.append(c['U'] * double)
        for i in range(L-1):
            for spin in (0, 1):
                a, b = 2*i+spin, 2*(i+1)+spin
                for src, dst in ((a,b), (b,a)):
                    if (m >> src) & 1 and not ((m >> dst) & 1):
                        lo, hi = sorted((src, dst))
                        between = ((1 << hi)-1) ^ ((1 << (lo+1))-1)
                        sign = -1 if (m & between).bit_count() % 2 else 1
                        rows.append(index[m ^ (1 << src) ^ (1 << dst)])
                        cols.append(k); values.append(-c['t']*sign)
    H = coo_matrix((values, (rows, cols)), shape=(len(masks),)*2).tocsr()
    rng = np.random.default_rng(417071)
    vals, vecs = eigsh(H, k=2, which='SA', tol=1e-12, maxiter=30000,
                      v0=rng.normal(size=len(masks)))
    order = np.argsort(vals); vals, vecs = vals[order], vecs[:,order]
    return {'energy': float(vals[0]), 'gap': float(vals[1]-vals[0]),
            'residuals': [float(np.linalg.norm(H@vecs[:,j]-vals[j]*vecs[:,j])) for j in range(2)],
            'dimension': len(masks), 'gap_scope': 'numerically resolved same-N,Sz sector; not a theorem certificate'}


def physical_rdm_sectors(psi, cut, floor):
    import numpy as np
    rho = psi.get_rho_segment(list(range(cut)))
    labels = list(rho.get_leg_labels())
    desired = [f'p{i}' for i in range(cut)] + [f'p{i}*' for i in range(cut)]
    a = rho.to_ndarray().transpose([labels.index(label) for label in desired]).reshape(4**cut, 4**cut)
    a = (a+a.conj().T)/2
    total = float(np.trace(a).real)
    n, sz = [], []
    for site in psi.sites[:cut]:
        n.append(np.rint(np.diag(site.get_op('Ntot').to_ndarray()).real).astype(int))
        sz.append(np.rint(2*np.diag(site.get_op('Sz').to_ndarray()).real).astype(int))
    groups = {}
    for row, states in enumerate(itertools.product(range(4), repeat=cut)):
        q = (sum(n[i][s] for i,s in enumerate(states)), sum(sz[i][s] for i,s in enumerate(states)))
        groups.setdefault(q, []).append(row)
    result = []
    negative = 0.
    for (N,S), rows in sorted(groups.items()):
        p = np.linalg.eigvalsh(a[np.ix_(rows,rows)])/total
        negative = min(negative, float(p.min()))
        p = np.sort(np.maximum(p,0.))[::-1]
        result.append({'N': int(N), 'twoSz': int(S), 'weight': float(p.sum()),
                       'probabilities': [float(x) for x in p if x > floor]})
    native = []
    for q, xi in psi.entanglement_spectrum(by_charge=True)[cut-1]:
        native.append({'charge': [int(x) for x in q],
                       'probabilities': sorted([float(np.exp(-x)) for x in xi if np.isfinite(x) and np.exp(-x)>floor], reverse=True)})
    # Charge conservation: qtotal_i = qphysical_i + Q_i - Q_{i+1}.
    # Hence physical prefix charge = Q_cut - Q_left + sum_{i<cut} qtotal_i.
    # This is a fixed tensor-conservation transformation, never a fit to levels.
    names = list(psi.chinfo.names)
    n_index = names.index('N')
    s_index = next(i for i,name in enumerate(names) if 'Sz' in name)
    left_charge = psi.get_B(0,form=None).get_leg('vL').get_charge(0)
    tensor_charge = np.sum([psi.get_B(i,form=None).qtotal for i in range(cut)],axis=0)
    offset = tensor_charge-left_charge
    assert psi.get_B(cut,form=None).get_leg('vL').qconj == 1
    for s in native:
        s['physical_charge'] = [int(x) for x in np.asarray(s['charge'])+offset]
    native_map = {(s['physical_charge'][n_index],s['physical_charge'][s_index]):s['probabilities'] for s in native}
    error = 0.
    for s in result:
        p, q = s['probabilities'], native_map.get((s['N'],s['twoSz']), [])
        d = max(len(p),len(q))
        if d:
            error = max(error, float(np.max(np.abs(np.pad(p,(0,d-len(p)))-np.pad(q,(0,d-len(q)))))))
    extra = set(native_map)-{(s['N'],s['twoSz']) for s in result}
    error = max(error, max((sum(native_map[q]) for q in extra), default=0.))
    return result, {'trace': total, 'minimum_eigenvalue': negative, 'native_charge_names': names,
                    'native_sectors': native, 'native_physical_probability_error': error,
                    'native_prefix_offset': [int(x) for x in offset],
                    'gauge_policy': 'Qphysical_prefix=Qnative_cut-Qleft+sum(prefix tensor qtotal); mismatch with physical RDM is UNRESOLVED; no fitted relabelling'}


def solve(c, arm, floor):
    import numpy as np
    from tenpy.models.hubbard import FermiHubbardChain
    from tenpy.networks.mps import MPS
    from tenpy.algorithms import dmrg
    start = time.monotonic()
    model = FermiHubbardChain({'L':c['L'], 't':c['t'], 'U':c['U'], 'mu':0., 'V':0.,
                              'cons_N':'N', 'cons_Sz':'Sz', 'bc_MPS':'finite', 'bc_x':'open'})
    assert c['Nup'] == c['Ndown']
    state = ['up','down']*c['Nup'] + ['empty']*(c['L']-c['Nup']-c['Ndown'])
    np.random.default_rng(arm['seed']).shuffle(state)
    psi = MPS.from_product_state(model.lat.mps_sites(), state, bc='finite')
    chi = arm['chi_max']
    ramp = {0:16, 2:32, 4:64, 6:128}
    if chi >= 256:
        ramp[8] = 256
    ramp = {k:min(v,chi) for k,v in ramp.items()}
    info = dmrg.run(psi, model, {'mixer':True, 'max_sweeps':arm['nsweeps'],
                    'min_sweeps':arm['nsweeps']-2, 'max_E_err':1e-12, 'max_S_err':1e-10,
                    'chi_list':ramp, 'trunc_params':{'chi_max':chi, 'svd_min':1e-14,
                                                  'trunc_cut':float(np.sqrt(arm['cutoff']))}})
    psi.canonical_form()
    energy = float(np.real(model.H_MPO.expectation_value(psi)))
    variance = float(np.real(model.H_MPO.variance(psi, exp_val=energy)))
    density = np.asarray(psi.expectation_value('Ntot')).real
    spin = np.asarray(psi.expectation_value('Sz')).real
    strings = []
    if 'P15' in c['catalogue_ids']:
        for i in range(c['L']):
            for j in range(i,c['L']):
                v = (-1.)**(j-i+1)*psi.expectation_value_multi_sites(['JW']*(j-i+1),i)
                strings.append({'i':i+1,'j':j+1,'value':float(np.real(v))})
    sectors, rdm = physical_rdm_sectors(psi,c['cut'],floor) if 'P11' in c['catalogue_ids'] else ([],{})
    return {'energy':energy, 'reported_energy':float(np.real(info['E'])), 'variance':variance,
            'norm':float(np.real(psi.overlap(psi))), 'N':float(density.sum()),
            'twoSz':float(2*spin.sum()), 'density':density.tolist(), 'sectors':sectors,
            'strings':strings, 'rdm_diagnostics':rdm, 'max_bond':int(max(psi.chi)),
            'elapsed_seconds':time.monotonic()-start}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')))
    p.add_argument('--output',required=True)
    p.add_argument('--pilot',action='store_true')
    args = p.parse_args()
    output = Path(args.output); output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    cfg = tomllib.loads(Path(args.manifest).read_text())
    result = {'schema':'tenpy_itensor_distinct_v1','backend':'TeNPy','status':'RUNNING',
              'manifest_sha256':sha(args.manifest),'oracle_source_sha256':sha(__file__),
              'versions':{x:importlib.metadata.version(x) for x in ('physics-tenpy','numpy','scipy')},
              'pilot':args.pilot,'cases':{}}
    for name,expected in {'physics-tenpy':cfg['tenpy_version'],'numpy':'2.4.6','scipy':'1.17.1'}.items():
        assert result['versions'][name] == expected, (name,result['versions'][name],expected)
    start=time.monotonic()
    try:
        for name,c in cfg['fixtures'].items():
            if args.pilot and name != 'doped': continue
            case={'ed':sector_ed(c),'arms':{}}
            result['cases'][name]=case
            for arm,a in cfg['arms'].items():
                if args.pilot and arm != 'coarse': continue
                case['arms'][arm]=solve(c,a,cfg['schmidt_probability_floor'])
        result['status']='PILOT_PARTIAL' if args.pilot else 'COMPLETED'
    except Exception:
        result['status']='EXECUTION_ERROR'
        result['error']=traceback.format_exc()
    result['elapsed_seconds']=time.monotonic()-start
    with output.open('x') as f: json.dump(result,f,indent=2,allow_nan=False)
    if result['status']=='EXECUTION_ERROR': raise SystemExit(2)

if __name__ == '__main__':
    main()

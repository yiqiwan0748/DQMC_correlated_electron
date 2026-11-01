#!/usr/bin/env python3
"""Standard-library complete-grid/provenance gate; never executes physics."""
import argparse,hashlib,json,math,pathlib,tomllib
HERE=pathlib.Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def finite(x):
    if isinstance(x,float) and not math.isfinite(x):raise ValueError('Nonfinite receipt')
    if isinstance(x,dict):
        for v in x.values():finite(v)
    if isinstance(x,list):
        for v in x:finite(v)

def error(a,b):
    if isinstance(a,list) or isinstance(b,list):
        if not isinstance(a,list) or not isinstance(b,list) or len(a)!=len(b):raise ValueError('Array shape mismatch')
        return max((error(x,y) for x,y in zip(a,b)),default=0.)
    return abs(a-b)

class Gates:
    def __init__(self):self.checks=[]
    def yes(self,name,value,eligibility=False):self.checks.append(dict(name=name,value=bool(value),passed=bool(value),eligibility=eligibility))
    def le(self,name,value,bound,eligibility=False):self.checks.append(dict(name=name,value=value,bound=bound,passed=value<=bound,eligibility=eligibility))
    def ge(self,name,value,bound,eligibility=False):self.checks.append(dict(name=name,value=value,bound=bound,passed=value>=bound,eligibility=eligibility))
    def result(self):
        status='UNRESOLVED' if any(x['eligibility'] and not x['passed'] for x in self.checks) else 'PASS' if all(x['passed'] for x in self.checks) else 'FAIL'
        return dict(status=status,checks=self.checks)

def keyed(rows,fields,expected):
    out={}
    for row in rows:
        key=tuple(row[f] for f in fields)
        if key in out:raise ValueError('Duplicate fixture '+str(key))
        out[key]=row
    if set(out)!=set(expected):raise ValueError('Incomplete or extra fixture grid')
    return out

def a15(m,t,o):
    p=m['A15'];g=Gates();oracle=o['A15']
    target=keyed(t['A15'],['U','Nup','Ndown'],[(u,*s) for u in p['interactions'] for s in p['spin_sectors']])
    for name in ['residual','orthogonality']:g.le('oracle_'+name,oracle[name],p['eigensystem_tolerance'],True)
    for key,row in target.items():
        label=str(key);g.yes(label+'_dimension',row['dimension']==8,True)
        g.le(label+'_centered_shift',abs(row['centered_shift']-1.5*row['U']),1e-12)
        g.le(label+'_physical_energies',error(row['energies'],oracle['energies']),p['energy_tolerance'])
        for name in ['residual','orthogonality','native_ground_residual']:g.le(label+'_'+name,row[name],p['eigensystem_tolerance'],True)
        g.yes(label+'_native_converged',row['native_converged'],True)
        g.le(label+'_native_ground_energy',abs(row['native_ground_energy_physical']-oracle['energies'][0]),p['energy_tolerance'])
        g.le(label+'_one_electron_double_occupancy',abs(row['native_double_occupancy']),1e-12)
        if len(row['groups'])!=len(oracle['groups']):raise ValueError('Degeneracy group mismatch')
        for j,(a,b) in enumerate(zip(row['groups'],oracle['groups'])):
            g.yes(label+f'_group{j}_multiplicity',a['multiplicity']==b['multiplicity'])
            for name in ['energies','projector','site_probabilities','orbital_weights']:
                g.le(label+f'_group{j}_'+name,error(a[name],b[name]),p['projector_tolerance'])
    return g.result()

def a16(m,t,o):
    p=m['A16'];g=Gates();oracle=o['A16'];meshes=p['meshes'];widths=p['smoothing_widths']
    rows=keyed(t['A16'],['Nk'],[(n,) for n in meshes]);ref=keyed(oracle['curves'],['sigma'],[(s,) for s in widths])
    g.le('oracle_bin_identity',oracle['elliptic_vs_reduced_BZ_bin_error'],p['oracle_identity_tolerance'],True)
    g.le('oracle_quad_error',oracle['maximum_reported_quad_error'],p['oracle_identity_tolerance'],True)
    g.le('oracle_total_mass',abs(oracle['normalization']-1),p['normalization_tolerance'],True)
    g.yes('held_out_mesh_is_final',meshes[-1]==p['held_out_mesh'],True)
    for (nk,),row in rows.items():
        g.yes(f'{nk}_counts',row['total_counts']==2*nk*nk,True)
        g.le(f'{nk}_energy_grid',error(row['energy'],oracle['energy']),1e-12)
        for field,value in [('minimum_band_energy',p['energy_min']),('maximum_band_energy',p['energy_max'])]:
            g.le(f'{nk}_'+field,abs(row[field]-value),p['spectrum_edge_tolerance'])
        curves=keyed(row['curves'],['sigma'],[(s,) for s in widths])
        for (sigma,),curve in curves.items():
            if len(curve['density'])!=p['bins']:raise ValueError('Incomplete DOS curve')
            g.le(f'{nk}_{sigma}_normalization',abs(curve['normalization']-1),p['normalization_tolerance'])
            g.ge(f'{nk}_{sigma}_positivity',min(curve['density']),-1e-14)
    final=rows[(meshes[-1],)];previous=rows[(meshes[-2],)]
    fc=keyed(final['curves'],['sigma'],[(s,) for s in widths]);pc=keyed(previous['curves'],['sigma'],[(s,) for s in widths])
    for s in widths:
        g.le(f'{s}_full_curve',error(fc[(s,)]['density'],ref[(s,)]['density']),p['curve_tolerance'])
        g.le(f'{s}_last_mesh_change',error(fc[(s,)]['density'],pc[(s,)]['density']),p['last_mesh_curve_tolerance'],True)
        g.le(f'{s}_peak_location',abs(fc[(s,)]['peak_energy']),p['peak_location_bins']*final['bin_width'])
    g.ge('narrower_broadening_peak_growth',fc[(widths[0],)]['peak']-fc[(widths[1],)]['peak'],p['minimum_peak_growth'])
    for name in ['lower_edge_mass','upper_edge_mass']:
        g.le(name+'_continuum',abs(final[name]-oracle[name]),p['edge_mass_tolerance'])
        g.le(name+'_last_mesh',abs(final[name]-previous[name]),p['last_mesh_edge_mass_tolerance'],True)
    return g.result()

def a17(m,t,o):
    p=m['A17'];g=Gates();grid=[(u,mu,b) for u in p['interactions'] for mu in p['centered_chemical_potentials'] for b in p['inverse_temperatures']]
    target=keyed(t['A17'],['U','mu','beta'],grid);oracle=keyed(o['A17'],['U','mu','beta'],grid)
    for key,a in target.items():
        b=oracle[key];label=str(key)
        g.yes(label+'_complete_basis',a['dimension']==16 and a['sector_count']==9,True)
        g.le(label+'_eigensystem',a['residual'],p['eigensystem_tolerance'],True)
        g.le(label+'_trace',abs(a['rho_trace']-1),p['factorization_tolerance'],True)
        g.ge(label+'_positive_rho',a['rho_minimum_eigenvalue'],-p['factorization_tolerance'],True)
        g.le(label+'_rho_offdiagonal',a['rho_offdiagonal'],p['factorization_tolerance'])
        g.le(label+'_product_state_probabilities',error(a['probabilities'],b['probabilities']),p['factorization_tolerance'])
        for name in ['density','double_occupancy','density_connected','spin_connected','local_spin_squared','local_density_variance','S_alm_q0','delta_alm']:
            g.le(label+'_'+name,error(a[name],b[name]),p['observable_tolerance'])
        g.ge(label+'_nonzero_local_spin_moment',min(a['local_spin_squared']),p['minimum_local_moment'])
        for name in ['density','double_occupancy','S_alm_q0','delta_alm']:
            g.le(label+'_native_'+name,error(a['native_'+name],b[name]),p['observable_tolerance'])
    return g.result()

def compare(target,oracle):
    m=tomllib.loads((HERE/'manifest.toml').read_text());t=tomllib.loads(pathlib.Path(target).read_text());o=json.loads(pathlib.Path(oracle).read_text())
    finite(t);finite(o)
    if t['release']!=m['release'] or o['release']!=m['release']:raise ValueError('Release mismatch')
    if t.get('mode')!='full':raise ValueError('Pilot cannot produce scientific completions')
    for r in [t,o]:
        if r['manifest_sha256']!=sha(HERE/'manifest.toml'):raise ValueError('Manifest mismatch')
    if t['driver_sha256']!=sha(HERE/'target.jl') or o['driver_sha256']!=sha(HERE/'oracle.py'):raise ValueError('Driver mismatch')
    if t['source_hashes']!=m['source_hashes']:raise ValueError('Recorded source mismatch')
    for p,h in m['source_hashes'].items():
        if sha(HERE/m['source_root']/p)!=h:raise ValueError('Local source mismatch')
    if t['julia_version']!=m['julia_version'] or t['blas_threads']!=1 or t['julia_threads']!=1:raise ValueError('Target runtime mismatch')
    if o['numpy_version']!=m['numpy_version'] or o['scipy_version']!=m['scipy_version']:raise ValueError('Oracle runtime mismatch')
    if tuple(map(int,o['python_version'].split('.')[:2]))<tuple(map(int,m['python_minimum'].split('.'))):raise ValueError('Python version mismatch')
    rows=[]
    for name,fn in [('A15',a15),('A16',a16),('A17',a17)]:
        try:r=fn(m,t,o)
        except Exception as ex:r=dict(status='UNRESOLVED',reason=str(ex))
        rows.append(dict(id=name,**r))
    return dict(release=m['release'],manifest_sha256=sha(HERE/'manifest.toml'),
                target_sha256=sha(target),oracle_sha256=sha(oracle),comparator_sha256=sha(__file__),
                scope=m['new_observable_scope'],benchmarks=rows,
                unsupported=[dict(id='A25',status='NOT_RUN',reason=m['A25']['reason'])])

def main():
    p=argparse.ArgumentParser();p.add_argument('--target',required=True);p.add_argument('--oracle',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    out=pathlib.Path(a.output)
    if out.exists():raise FileExistsError(out)
    try:r=compare(a.target,a.oracle)
    except Exception as ex:r=dict(benchmarks=[dict(id=i,status='UNRESOLVED',reason=str(ex)) for i in ['A15','A16','A17']])
    out.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    for row in r['benchmarks']:print(row['id'],row['status'])
    if not all(r['status']=='PASS' for r in r['benchmarks']):raise SystemExit(1)

if __name__=='__main__':main()

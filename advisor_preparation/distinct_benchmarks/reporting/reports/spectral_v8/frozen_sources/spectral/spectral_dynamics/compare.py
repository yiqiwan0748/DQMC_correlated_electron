#!/usr/bin/env python3
"""Frozen complete-grid spectral comparisons; stdlib postprocessing only."""
import argparse,hashlib,json,math,tomllib
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def finite(x):
    if isinstance(x,dict):return all(finite(y) for y in x.values())
    if isinstance(x,list):return all(finite(y) for y in x)
    return not isinstance(x,float) or math.isfinite(x)
def require(ok,label):
    if not ok:raise ValueError(label)
def lookup(rows,keys,expected):
    d={tuple(r[k] for k in keys):r for r in rows}
    require(len(d)==len(rows) and set(d)==set(expected),'Missing/extra/duplicate grid point')
    return d

def evaluate(cfg,o,t):
    g=cfg['gates'];require(finite(o) and finite(t),'Nonfinite receipt');result=[]
    global_reasons=[]
    for name,data in [('QuSpin',o),('user_ED',t)]:
        if data['status']!='COMPLETED' or data['schema']!=cfg['schema']:global_reasons.append(name+': incomplete execution/schema')
        if max(data['diagnostics'].values())>g['normalization_abs']:global_reasons.append(name+': Hamiltonian Hermiticity/number conservation')
    for idea in ('P08','P09','P10'):
        reasons=list(global_reasons);metrics={};different=False;s=cfg[idea]
        if idea=='P08':
            times=[s['time_start']+(s['time_stop']-s['time_start'])*i/(s['time_count']-1) for i in range(s['time_count'])]
            for name,data in [('QuSpin',o),('user_ED',t)]:
                rows=data[idea]['rows'];require([r['time'] for r in rows]==times,'Incomplete real-time grid')
                initial=rows[0]
                if max(abs(initial['density'][j]-(2. if j+1 in s['initial_doublon_sites_one_based'] else 0.)) for j in range(cfg['L']))>g['normalization_abs']:
                    reasons.append(name+': initial CDW profile')
                if abs(initial['imbalance']-1)>g['normalization_abs'] or abs(initial['double_occupancy_total']-2)>g['normalization_abs']:
                    reasons.append(name+': initial normalization')
                if abs(data[idea]['initial_energy']-cfg['U']*len(s['initial_doublon_sites_one_based']))>g['conservation_abs']:
                    reasons.append(name+': initial uncentered physical energy')
                if max(r['imbalance'] for r in rows)-min(r['imbalance'] for r in rows)<s['minimum_dynamic_range']:reasons.append(name+': nontrivial quench not resolved')
                for r in rows:
                    require(len(r['density'])==cfg['L'],'Incomplete density profile')
                    if abs(r['norm']-1)>g['normalization_abs'] or abs(r['energy']-data[idea]['initial_energy'])>g['conservation_abs']:
                        reasons.append(name+': unitary norm/energy conservation')
                    if abs(r['Nup']-s['Nup'])>g['conservation_abs'] or abs(r['Ndown']-s['Ndown'])>g['conservation_abs'] or abs(sum(r['density'])-s['Nup']-s['Ndown'])>g['conservation_abs']:
                        reasons.append(name+': conserved electron counts')
                    if min(r['density']) < -g['positivity_roundoff'] or max(r['density'])>2+g['positivity_roundoff'] or not -g['positivity_roundoff']<=r['double_occupancy_total']<=2+g['positivity_roundoff']:
                        reasons.append(name+': density/doublon positivity')
            if o[idea]['exp_action_vs_dense_spectral_error']>g['exp_action_spectral_abs']:reasons.append('QuSpin exponential-action accuracy')
            if t[idea]['Hamiltonian_eigen_residual']>g['eigen_residual_max']:reasons.append('User-sector eigensystem residual')
            for k in ('imbalance','double_occupancy_total','norm','energy','Nup','Ndown'):
                metrics[k]=max(abs(a[k]-b[k]) for a,b in zip(o[idea]['rows'],t[idea]['rows']))
            metrics['density']=max(abs(x-y) for a,b in zip(o[idea]['rows'],t[idea]['rows']) for x,y in zip(a['density'],b['density']))
            different=any(x>g['cross_observable_abs'] for x in metrics.values())
        elif idea=='P09':
            expected=[(spin,m) for spin in s['spin'] for m in s['momentum_indices']]
            maps=[]
            for name,data in [('QuSpin',o),('user_ED',t)]:
                z=data[idea]['ground'];rows=lookup(data[idea]['rows'],('spin','momentum_index'),expected);maps.append(rows)
                if max(data[idea]['operator_checks'].values())>g['sum_rule_abs']:
                    reasons.append(name+': cross-spin fermionic anticommutation')
                if z['residual']>g['eigen_residual_max'] or z['gap_same_sector']<=g['ground_gap_min'] or abs(z['norm']-1)>g['normalization_abs']:
                    reasons.append(name+': unique ground-state/residual gate')
                if abs(z['Nup']-s['Nup'])>g['normalization_abs'] or abs(z['Ndown']-s['Ndown'])>g['normalization_abs']:reasons.append(name+': canonical ground sector')
                for r in rows.values():
                    a,b=r['addition'],r['removal']
                    require(a and b,'Empty spectral energy groups')
                    for branch in (a,b):
                        require([z['omega'] for z in branch]==sorted(z['omega'] for z in branch),'Unsorted spectral groups')
                        require(sum(z['multiplicity'] for z in branch)==24,'Incomplete N±1 eigenspectrum')
                        if any(z['span']>s['degeneracy_cluster_abs'] or z['weight']< -g['positivity_roundoff'] or z['multiplicity']<1 for z in branch):reasons.append(name+': spectral cluster eligibility')
                    if r['addition_dimension']!=24 or r['removal_dimension']!=24 or max(r['addition_eigen_residual'],r['removal_eigen_residual'])>g['eigen_residual_max']:
                        reasons.append(name+': full spin-resolved N±1 sectors')
                    if min(z['omega'] for z in a)< -g['sum_rule_abs'] or max(z['omega'] for z in b)>g['sum_rule_abs']:reasons.append(name+': addition/removal frequency convention')
                    aw=sum(z['weight'] for z in a);rw=sum(z['weight'] for z in b)
                    first=sum(z['weight']*z['omega'] for z in a+b);second=sum(z['weight']*z['omega']**2 for z in a+b)
                    errors=[abs(aw+rw-1),abs(r['sum_weight']-1),abs(aw-r['addition_norm']),abs(rw-r['occupation']),
                            abs(first-r['first_moment']),abs(second-r['second_moment']),
                            abs(r['first_moment']-r['first_moment_direct']),abs(r['first_moment']-r['first_moment_expected']),
                            abs(r['second_moment']-r['second_moment_direct']),r['anticommutator_matrix_error']]
                    if max(errors)>g['sum_rule_abs']:reasons.append(name+': fermionic normalization/occupation/moment sum rule')
            z=t[idea]['ground']
            if max(z['dense_ground_overlap_error'],z['dense_ground_energy_error'])>g['sum_rule_abs']:reasons.append('User Lanczos/dense ground-state consistency')
            metrics.update({'ground_'+k:abs(o[idea]['ground'][k]-t[idea]['ground'][k]) for k in ('energy','gap_same_sector','Nup','Ndown')})
            pole_error=weight_error=0.
            for key in expected:
                a,b=maps[0][key],maps[1][key]
                for branch in ('addition','removal'):
                    if len(a[branch])!=len(b[branch]):different=True;continue
                    for x,y in zip(a[branch],b[branch]):
                        pole_error=max(pole_error,abs(x['omega']-y['omega']));weight_error=max(weight_error,abs(x['weight']-y['weight']))
                        if x['multiplicity']!=y['multiplicity']:different=True
                for k in ('occupation','addition_norm','sum_weight','first_moment','second_moment'):
                    metrics[k]=max(metrics.get(k,0.),abs(a[k]-b[k]))
            metrics.update(pole=pole_error,residue=weight_error)
            different|=pole_error>g['pole_abs'] or weight_error>g['residue_abs'] or any(v>g['cross_observable_abs'] for k,v in metrics.items() if k not in ('pole','residue'))
        else:
            expected=[(b,m,f) for b in s['beta'] for m in s['momentum_indices'] for f in s['tau_fractions']];maps=[]
            for name,data in [('QuSpin',o),('user_ED',t)]:
                rows=lookup(data[idea]['rows'],('beta','momentum_index','tau_fraction'),expected);maps.append(rows)
                thermal=lookup(data[idea]['thermal'],('beta',),[(b,) for b in s['beta']])
                dynamic=0.
                for beta in s['beta']:
                    z=thermal[beta,]
                    if abs(z['trace']-1)>g['normalization_abs'] or z['rho_expm_error']>g['normalization_abs']:reasons.append(name+': Gibbs trace/exponential consistency')
                    if abs(z['mu_physical']-s['physical_mu'])>g['normalization_abs'] or abs(z['mu_centered']-(s['physical_mu']-cfg['U']/2))>g['normalization_abs']:reasons.append(name+': physical/centered chemical potential')
                    if z['spin_mean_abs_max']>g['normalization_abs']:reasons.append(name+': zero-field spin symmetry')
                    total=sum(rows[beta,m,0.]['value'] for m in s['momentum_indices'])
                    if abs(total-z['local_moment_sum'])>g['sum_rule_abs'] or abs(z['local_moment_sum']-(z['Nmean']/4-z['double_occupancy_total']/2))>g['sum_rule_abs']:
                        reasons.append(name+': physical-Sz local-moment sum rule')
                    for m in s['momentum_indices']:
                        seq=[rows[beta,m,f]['value'] for f in s['tau_fractions']]
                        if m:dynamic=max(dynamic,max(seq)-min(seq))
                        for fraction in s['tau_fractions']:
                            r=rows[beta,m,fraction];mirror=rows[beta,(-m)%cfg['L'],1-fraction]
                            if abs(r['tau']-beta*fraction)>1e-12:reasons.append(name+': imaginary-time coordinate')
                            if min(r['value'],r['spectral_weight_min'])< -g['positivity_roundoff']:reasons.append(name+': response/spectral positivity')
                            if abs(r['unconnected_value']-r['spin_mean_abs']**2-r['value'])>g['sum_rule_abs']:reasons.append(name+': connected subtraction')
                            if abs(r['value']-mirror['value'])>g['endpoint_kms_abs']:reasons.append(name+': bosonic KMS relation')
                            if fraction in (0.,1.) and abs(r['value']-r['equal_time_direct'])>g['endpoint_kms_abs']:reasons.append(name+': equal-time endpoint')
                            if m==0 and abs(r['value']-seq[0])>g['endpoint_kms_abs']:reasons.append(name+': conserved uniform spin response')
                if dynamic<s['minimum_dynamic_range']:reasons.append(name+': nontrivial imaginary-time response unresolved')
            for key in expected:
                for field in ('value','equal_time_direct','spin_mean_abs','unconnected_value'):
                    metrics[field]=max(metrics.get(field,0.),abs(maps[0][key][field]-maps[1][key][field]))
            a=lookup(o[idea]['thermal'],('beta',),[(b,) for b in s['beta']]);b=lookup(t[idea]['thermal'],('beta',),[(b,) for b in s['beta']])
            for key in a:
                for field in ('Nmean','double_occupancy_total','local_moment_sum'):
                    metrics[field]=max(metrics.get(field,0.),abs(a[key][field]-b[key][field]))
            different=any(v>g['cross_observable_abs'] for v in metrics.values())
        status='UNRESOLVED' if reasons else ('FAIL' if different else 'PASS')
        result.append({'id':idea,'status':status,'max_abs_errors':metrics,'eligibility_reasons':sorted(set(reasons)),
                       'scope':'Frozen L4 fixture and new observable adapters only; native DQMC dynamics not tested.'})
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(Path(__file__).with_name('manifest.toml')))
    for key in ('oracle','target','output'):p.add_argument('--'+key,required=True)
    args=p.parse_args();cfg=tomllib.loads(Path(args.manifest).read_text());digest=sha(args.manifest)
    try:
        o=json.loads(Path(args.oracle).read_text());t=tomllib.loads(Path(args.target).read_text())
        require(o['manifest_sha256']==t['manifest_sha256']==digest,'Manifest receipt differs')
        require(o['oracle_source_sha256']==sha(Path(__file__).with_name('oracle_quspin.py')),'Oracle adapter differs')
        require(t['adapter_sha256']==sha(Path(__file__).with_name('target_ed.jl')),'Target adapter differs')
        require(o['versions']=={x:cfg[x+'_version'] for x in ('quspin','numpy','scipy')},'Oracle versions differ')
        require(t['julia_version']==cfg['julia_version'],'Julia version differs')
        for key in ('thermal_source_sha256','ed_source_sha256'):require(t[key]==cfg[key],'User source closure differs')
        rows=evaluate(cfg,o,t)
    except Exception as exc:rows=[{'id':x,'status':'UNRESOLVED','reason':f'{type(exc).__name__}: {exc}'} for x in ('P08','P09','P10')]
    result={'schema':cfg['schema'],'manifest_sha256':digest,'rows':rows,'distinct_ideas_passed':sum(r['status']=='PASS' for r in rows)}
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Independent stdlib-only audit of frozen spectral receipts; no physics runs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import tarfile
import tomllib

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE_SHA='3fb6ce6accac61140ee21584dddbb79418b121f7d62ac9391b31996f1cf2239b'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text()) if p.suffix=='.json' else tomllib.loads(p.read_text())
def require(ok,msg):
    if not ok:raise ValueError(msg)
def finite(x):
    if isinstance(x,dict):return all(finite(v) for v in x.values())
    if isinstance(x,list):return all(finite(v) for v in x)
    return not isinstance(x,float) or math.isfinite(x)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();require(not a.output.exists(),'Use a new audit output')
    own=ROOT/'spectral_dynamics';result=a.result_dir
    freeze_path=ROOT/'source_releases/spectral-dynamics-v1-freeze.json'
    archive=ROOT/'source_releases/spectral-dynamics-v1-source.tar.gz'
    freeze=load(freeze_path);require(sha(archive)==ARCHIVE_SHA==freeze['archive_sha256'],'Frozen archive identity')
    cfg=load(own/'manifest.toml');g=cfg['gates'];seen={}
    with tarfile.open(archive) as z:
        for m in z:
            if m.isdir():continue
            require(m.isfile(),'Non-regular frozen member')
            path=m.name if m.name in freeze['files'] else '/'.join(Path(m.name).parts[1:]);data=z.extractfile(m).read()
            digest=hashlib.sha256(data).hexdigest()
            require(path in freeze['files'] and path not in seen and digest==freeze['files'][path],'Frozen member identity')
            require((ROOT/path).read_bytes()==data,'Authored source differs from frozen archive')
            seen[path]=digest
    require(seen==freeze['files'],'Incomplete frozen closure')
    paths=[result/n for n in ('oracle.json','target.toml','comparison.json','comparison_independent_replay_attempt01.json','retrieval_sha256.txt')]
    o,t,reported,replay=[load(p) for p in paths[:4]]
    require(finite(o) and finite(t),'Nonfinite raw data')
    require(reported==replay and replay['distinct_ideas_passed']==3,'Comparator replay differs')
    require([(r['id'],r['status'],r['eligibility_reasons']) for r in replay['rows']]==[(i,'PASS',[]) for i in ('P08','P09','P10')], 'Comparator eligibility or ID set')
    require(o['manifest_sha256']==t['manifest_sha256']==sha(own/'manifest.toml'),'Manifest receipt')
    require(o['oracle_source_sha256']==sha(own/'oracle_quspin.py') and t['adapter_sha256']==sha(own/'target_ed.jl'),'Adapter receipts')
    require(o['versions']=={k:cfg[k+'_version'] for k in ('quspin','numpy','scipy')} and t['julia_version']==cfg['julia_version'],'Runtime versions')
    sources=[]
    for key,name in [('ed_source_sha256','compute_s2_ed.jl'),('thermal_source_sha256','audit_u4_l4_thermal_ed.jl')]:
        p=ROOT.parent/'validated_package/HubbardRenyiDQMC/analysis'/name;sources.append(p)
        require(sha(p)==cfg[key]==t[key],'Actual user source identity')
    remote={Path(line.split(None,1)[1].strip()).name:line.split(None,1)[0] for line in paths[4].read_text().splitlines()}
    for p in paths[:3]+sources+[result/'oracle.resources.txt',result/'target.resources.txt',ROOT/'logs/spectral_dynamics_6458468.log']:
        require(remote[p.name]==sha(p),'ISAAC transfer or actual source mismatch: '+p.name)
    diagnostics={}
    def bound(name,value,tol):
        require(math.isfinite(value) and value<=tol,name)
        diagnostics[name]=max(diagnostics.get(name,0.),value)
    for name,d in [('QuSpin',o),('user_ED',t)]:
        require(d['status']=='COMPLETED' and d['schema']==cfg['schema'],'Incomplete execution')
        for k,v in d['diagnostics'].items():bound(name+'/'+k,abs(v),g['normalization_abs'])
        rows=d['P08']['rows'];s=cfg['P08'];times=[s['time_start']+(s['time_stop']-s['time_start'])*i/(s['time_count']-1) for i in range(s['time_count'])]
        require(len(rows)==25 and [r['time'] for r in rows]==times,'P08 complete time grid')
        bound(name+'/initial_energy',abs(d['P08']['initial_energy']-2*cfg['U']),g['conservation_abs'])
        for r in rows:
            require(len(r['density'])==4,'P08 density shape')
            for k,want in [('norm',1.),('Nup',2.),('Ndown',2.),('energy',d['P08']['initial_energy'])]:bound(name+'/P08_'+k,abs(r[k]-want),g['normalization_abs'] if k=='norm' else g['conservation_abs'])
            bound(name+'/P08_density_number',abs(sum(r['density'])-4),g['conservation_abs'])
            bound(name+'/P08_imbalance_identity',abs(r['imbalance']-sum((-1 if j%2==0 else 1)*v for j,v in enumerate(r['density']))/4),g['conservation_abs'])
            require(min(r['density'])>=-g['positivity_roundoff'] and max(r['density'])<=2+g['positivity_roundoff'],'P08 density positivity')
            require(-g['positivity_roundoff']<=r['double_occupancy_total']<=2+g['positivity_roundoff'],'P08 doublon bound')
        require(max(r['imbalance'] for r in rows)-min(r['imbalance'] for r in rows)>1e-5,'P08 nontrivial quench')
        for actual,want in zip(rows[0]['density'],[0.,2.,0.,2.]):bound(name+'/P08_initial_density',abs(actual-want),g['normalization_abs'])
        for k,want in [('imbalance',1.),('double_occupancy_total',2.)]:bound(name+'/P08_initial_'+k,abs(rows[0][k]-want),g['normalization_abs'])
        if name=='QuSpin':bound('QuSpin/exp_action_spectral',d['P08']['exp_action_vs_dense_spectral_error'],g['exp_action_spectral_abs'])
        else:bound('user_ED/P08_eigen_residual',d['P08']['Hamiltonian_eigen_residual'],g['eigen_residual_max'])
        ground=d['P09']['ground'];bound(name+'/P09_residual',ground['residual'],g['eigen_residual_max'])
        require(ground['gap_same_sector']>g['ground_gap_min'],'P09 unique ground gap')
        for k,want in [('norm',1.),('Nup',2.),('Ndown',2.)]:bound(name+'/P09_ground_'+k,abs(ground[k]-want),g['normalization_abs'])
        if name=='user_ED':
            for k in ('dense_ground_energy_error','dense_ground_overlap_error'):bound(name+'/'+k,ground[k],g['sum_rule_abs'])
        require(set(d['P09']['operator_checks'])=={'cross_spin_annihilation_anticommutator','cross_spin_mixed_anticommutator'},'P09 missing cross-spin CAR')
        for k,v in d['P09']['operator_checks'].items():bound(name+'/'+k,v,g['sum_rule_abs'])
        spectra=d['P09']['rows'];require(len(spectra)==8 and {(r['spin'],r['momentum_index']) for r in spectra}=={(s,m) for s in ('up','down') for m in range(4)},'P09 full channels')
        for r in spectra:
            for side in ('addition','removal'):
                groups=r[side];require(groups and sum(x['multiplicity'] for x in groups)==24 and r[side+'_dimension']==24,'P09 complete sector')
                require([x['omega'] for x in groups]==sorted(x['omega'] for x in groups),'P09 sorted clusters')
                require(all(x['multiplicity']>0 and x['span']<=cfg['P09']['degeneracy_cluster_abs'] and x['weight']>=-g['positivity_roundoff'] for x in groups),'P09 cluster eligibility')
                bound(name+'/P09_'+side+'_eigen_residual',r[side+'_eigen_residual'],g['eigen_residual_max'])
            add,rem=r['addition'],r['removal'];wa=sum(x['weight'] for x in add);wr=sum(x['weight'] for x in rem)
            require(min(x['omega'] for x in add)>=-g['sum_rule_abs'] and max(x['omega'] for x in rem)<=g['sum_rule_abs'],'P09 frequency branches')
            moment=sum(x['omega']*x['weight'] for x in add+rem);second=sum(x['omega']**2*x['weight'] for x in add+rem)
            expected=-2*cfg['t']*math.cos(2*math.pi*r['momentum_index']/4)+cfg['U']*.5-cfg['P09']['physical_mu_frequency']
            checks={'total_weight':abs(wa+wr-1),'reported_weight':abs(r['sum_weight']-1),'addition_norm':abs(wa-r['addition_norm']),'occupation':abs(wr-r['occupation']),
                'first_cluster':abs(moment-r['first_moment']),'first_direct':abs(moment-r['first_moment_direct']),'first_Hubbard':abs(moment-expected),
                'second_cluster':abs(second-r['second_moment']),'second_direct':abs(second-r['second_moment_direct']),'same_spin_CAR':r['anticommutator_matrix_error']}
            for k,v in checks.items():bound(name+'/P09_'+k,v,g['sum_rule_abs'])
        rr=d['P10']['rows'];s=cfg['P10'];grid={(b,m,f) for b in s['beta'] for m in range(4) for f in s['tau_fractions']}
        lookup={(r['beta'],r['momentum_index'],r['tau_fraction']):r for r in rr}
        require(len(rr)==len(lookup)==len(grid)==108 and set(lookup)==grid,'P10 full thermal grid')
        thermal={r['beta']:r for r in d['P10']['thermal']};require(len(d['P10']['thermal'])==len(thermal)==3 and set(thermal)==set(s['beta']),'P10 thermal cases')
        dynamic=0.
        for beta,z in thermal.items():
            for k,v in [('trace',abs(z['trace']-1)),('rho_expm',z['rho_expm_error']),('physical_mu',abs(z['mu_physical']-s['physical_mu'])),('centered_mu',abs(z['mu_centered']-s['physical_mu']+cfg['U']/2)),('spin_mean',z['spin_mean_abs_max'])]:bound(name+'/P10_'+k,v,g['normalization_abs'])
            total=sum(lookup[beta,m,0.]['value'] for m in range(4));local=z['Nmean']/4-z['double_occupancy_total']/2
            bound(name+'/P10_local_moment',max(abs(total-local),abs(z['local_moment_sum']-local)),g['sum_rule_abs'])
            for m in range(4):
                seq=[lookup[beta,m,f]['value'] for f in s['tau_fractions']]
                if m:dynamic=max(dynamic,max(seq)-min(seq))
                for f in s['tau_fractions']:
                    r=lookup[beta,m,f];mirror=lookup[beta,(-m)%4,1-f]
                    bound(name+'/P10_tau',abs(r['tau']-beta*f),1e-12)
                    require(min(r['value'],r['spectral_weight_min'])>=-g['positivity_roundoff'],'P10 positivity')
                    bound(name+'/P10_connected',abs(r['unconnected_value']-r['spin_mean_abs']**2-r['value']),g['sum_rule_abs'])
                    bound(name+'/P10_KMS',abs(r['value']-mirror['value']),g['endpoint_kms_abs'])
                    if f in (0.,1.):bound(name+'/P10_endpoint',abs(r['value']-r['equal_time_direct']),g['endpoint_kms_abs'])
                    if m==0:bound(name+'/P10_uniform_conservation',abs(r['value']-seq[0]),g['endpoint_kms_abs'])
        require(dynamic>=s['minimum_dynamic_range'],'P10 nontrivial dynamics')
    paths += [freeze_path,archive,*[own/p for p in ('manifest.toml','oracle_quspin.py','target_ed.jl','compare.py')],*sources,ROOT/'logs/spectral_dynamics_6458468.log',Path(__file__)]
    receipt={'schema':'spectral_receipt_independent_audit_v1','status':'PASS','catalogue_ids':['P08','P09','P10'],
        'method':'Independent stdlib provenance/completeness and physical-identity audit plus separately executed frozen comparator replay; no physics',
        'source_archive_sha256':ARCHIVE_SHA,'manifest_sha256':sha(own/'manifest.toml'),'auditor_sha256':sha(__file__),
        'input_sha256':{p.resolve().relative_to(ROOT.parent).as_posix():sha(p) for p in paths},
        'grid_counts':{'P08_time_points':25,'P09_spin_momentum_channels':8,'P09_each_Nplus_or_Nminus_states':24,'P10_beta_q_tau_points':108},
        'independent_guard_maxima':diagnostics,'cross_backend_results':replay['rows'],
        'scientific_runs_executed_by_auditor':0,'scope':'Actual user ED/thermal sources plus new dynamics/operator adapters, not native DQMC or DMRG dynamics'}
    a.output.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':receipt['status'],'grid_counts':receipt['grid_counts'],'auditor_sha256':receipt['auditor_sha256']},indent=2))

if __name__=='__main__':main()

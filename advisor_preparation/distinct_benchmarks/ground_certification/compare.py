#!/usr/bin/env python3
"""External ground-projector certification and prescribed-ensemble comparison."""
import argparse,hashlib,json,math,os,tomllib
from pathlib import Path
OWN=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,reason):
    if not ok:raise ValueError(reason)
def finite(x):
    if isinstance(x,dict):return all(finite(v) for v in x.values())
    if isinstance(x,list):return all(finite(v) for v in x)
    return not isinstance(x,float) or math.isfinite(x)
def error(a,b):
    if isinstance(a,list) or isinstance(b,list):
        require(isinstance(a,list) and isinstance(b,list) and len(a)==len(b),'Array shape mismatch')
        return max((error(x,y) for x,y in zip(a,b)),default=0.)
    return abs(a-b)
def matvec(A,v):return [sum(a*x for a,x in zip(row,v)) for row in A]
def dot(a,b):return sum(x.conjugate()*y for x,y in zip(a,b))

def x09(cfg,o,t):
    s=cfg['X09'];g=cfg['gates'];reasons=[];m={};checks=[]
    expected={(u,d) for u in range(16) for d in range(16) if u.bit_count()==2 and d.bit_count()==2}
    oq=[tuple(x) for x in o['coordinates']];tq=[tuple(x) for x in t['coordinates']]
    require(len(oq)==len(tq)==36 and set(oq)==set(tq)==expected,'X09 complete fixed-spin coordinates')
    idx={x:i for i,x in enumerate(oq)};permutation=[idx[x] for x in tq]
    for key in ('Hamiltonian','ground_projector'):
        require(len(o[key])==36 and all(len(r)==36 for r in o[key]),'X09 oracle matrix shape')
    H=[[o['Hamiltonian'][i][j] for j in permutation] for i in permutation]
    P=[[o['ground_projector'][i][j] for j in permutation] for i in permutation]
    for key in ('Hamiltonian','MPO_qn_Hamiltonian','MPO_dense_Hamiltonian'):
        m[key]=error(t[key],H)
    for key in ('MPO_qn_max_imaginary','MPO_dense_max_imaginary'):m[key]=abs(t[key])
    if o['ground_dimension']!=1 or o['gap']<=g['ground_gap_min']:reasons.append('Independent unique ground/gap unresolved')
    for key,value in [('oracle_eigen_residual',o['eigen_residual']),('source_eigen_residual',t['source_eigen_residual']),('oracle_projector_idempotency',o['projector_idempotency'])]:
        if value>g['eigen_residual_max']:reasons.append(key)
    wanted=s['starts']+s['challenge_states'];rows={r['label']:r for r in t['rows']}
    require(len(rows)==len(t['rows'])==6 and set(rows)==set(wanted),'X09 full starts/challenges')
    different=any(v>g['matrix_abs'] for v in m.values())
    for name in wanted:
        r=rows[name];kind='native_dmrg' if name in s['starts'] else 'manufactured_excited_mps'
        require(r['kind']==kind,'X09 candidate attribution')
        require(len(r['coeff_real'])==len(r['coeff_imag'])==36,'X09 full state vector')
        v=[complex(x,y) for x,y in zip(r['coeff_real'],r['coeff_imag'])];norm=dot(v,v).real
        require(norm>0,'X09 zero state')
        fidelity=dot(v,matvec(P,v)).real/norm
        hv=matvec(H,v);E=dot(v,hv).real/norm;E2=dot(hv,hv).real/norm;variance=E2-E*E
        for key,value in [('norm',abs(r['norm']-1)),('sector_norm',abs(norm-r['norm'])),('outside_sector',abs(r['outside_sector_weight'])),('Nup',abs(r['Nup']-2)),('Ndown',abs(r['Ndown']-2))]:
            if value>g['norm_abs']:reasons.append(name+': '+key)
        for key,value in [('physical_energy',abs(E-r['energy'])),('physical_second_moment',abs(E2-r['energy_second_moment'])),('physical_variance',abs(variance-r['variance']))]:
            m[key]=max(m.get(key,0.),value)
            different|=value>g['energy_abs']
        low_variance=g['variance_min']<=r['variance']<=g['variance_max']
        excess=E-o['ground_energy']
        certificate=low_variance and -g['energy_abs']<=excess<=g['energy_abs'] and fidelity>=g['ground_fidelity_min']
        if kind=='native_dmrg':
            pos=s['starts'].index(name)
            require(r['seed']==s['seeds'][pos] and len(r['truncation_errors'])==len(s['maxdim_schedule']),'X09 actual seeded sweep record')
            different|=abs(r['reported_energy']-r['energy'])>g['energy_abs']
            if not low_variance:reasons.append(name+': native DMRG variance unresolved')
            different|=not certificate
        else:
            if r['source_eigen_residual']>g['eigen_residual_max']:reasons.append(name+': challenge source residual')
            different|=not low_variance or certificate or fidelity>g['challenge_fidelity_max'] or excess<.5*o['gap']
            m['challenge_energy_vs_source']=max(m.get('challenge_energy_vs_source',0.),abs(E-r['source_eigenvalue']))
            different|=m['challenge_energy_vs_source']>g['energy_abs']
        checks.append({'label':name,'kind':kind,'external_ground_fidelity':fidelity,'external_energy_excess':excess,
                       'variance':r['variance'],'low_variance_gate_pass':low_variance,'external_ground_certificate_pass':certificate})
    return {'id':'X09','status':'UNRESOLVED' if reasons else ('FAIL' if different else 'PASS'),
        'max_abs_errors':m,'eligibility_reasons':sorted(set(reasons)),'candidate_certificates':checks,
        'scope':'Four actual native-DMRG starts must be independently ground certified. Two manufactured excited MPS diagnostic states must be rejected; they do not demonstrate natural optimizer failure.'}

def x12(cfg,o,t):
    s=cfg['X12'];g=cfg['gates'];reasons=[];errors={};different=False
    keys=[(name,-1.) for name in s['prescriptions']]+[('canonical_gibbs',b) for b in s['beta']]
    theta=math.pi/s['rotation_theta_pi_denominator'];phi=math.pi/s['rotation_phi_pi_denominator'];per_backend={}
    for backend,data in [('QuSpin',o),('user_ED_thermal_RDM',t)]:
        require(data['canonical_dimension']==4 and len(data['canonical_eigenvalues'])==4,'X12 all canonical-N3 states')
        delta=error(data['canonical_eigenvalues'],[3.,3.,5.,5.]);errors['analytic_canonical_eigenvalues']=max(errors.get('analytic_canonical_eigenvalues',0.),delta)
        different|=delta>g['energy_abs']
        require(data['ground_dimension']==s['expected_ground_dimension']==2,'X12 ground doublet')
        different|=abs(data['gap']-2)>g['energy_abs'] or abs(data['spin_lowering_norm']-1)>g['observable_abs']
        if data['eigen_residual']>g['eigen_residual_max']:reasons.append(backend+': eigensystem residual')
        if data['mixture_rotation_error']>g['matrix_abs']:reasons.append(backend+': mixed-ground basis invariance')
        rows={(r['label'],r['beta']):r for r in data['rows']}
        require(len(rows)==len(data['rows'])==len(keys) and set(rows)==set(keys),'X12 full prescription/beta grid')
        per_backend[backend]=rows
        for key,r in rows.items():
            label,beta=key;pure=label.startswith('pure_')
            eps=1/(1+math.exp(2*beta)) if beta>=0 else 0.
            expected_E=3+2*eps;expected_distance=1. if pure else 2*eps
            expected_spin=.25 if label=='pure_up' else math.cos(2*theta)/4 if label=='pure_rotated' else 0.
            connected=-expected_spin*expected_spin
            for field,want in [('energy',expected_E),('trace',1.),('Nmean',3.),('Czz',0.),('connected_Czz',connected),('distance_to_equal_ground',expected_distance)]:
                delta=abs(r[field]-want);errors['analytic_'+field]=max(errors.get('analytic_'+field,0.),delta);different|=delta>g['observable_abs']
            delta=error(r['spin_z_means'],[expected_spin,expected_spin]);errors['analytic_spin_z_means']=max(errors.get('analytic_spin_z_means',0.),delta)
            different|=delta>g['observable_abs']
            xyz=[0.,0.,.25] if label=='pure_up' else [math.sin(2*theta)*math.cos(phi)/4,math.sin(2*theta)*math.sin(phi)/4,expected_spin] if label=='pure_rotated' else [0.,0.,0.]
            delta=error(r['spin_xyz_means'],[xyz,xyz]);errors['analytic_spin_xyz_means']=max(errors.get('analytic_spin_xyz_means',0.),delta)
            different|=delta>g['observable_abs']
            if beta>=0:different|=abs(r['excited_weight']-eps)>g['observable_abs']
            rdms={x['nA']:x for x in r['rdms']};require(len(r['rdms'])==len(rdms)==2 and set(rdms)=={1,2},'X12 complete RDM cuts')
            for nA,z in rdms.items():
                require(len(z['spectrum'])==4**nA and z['spectrum']==sorted(z['spectrum'],reverse=True),'X12 complete sorted RDM spectrum')
                if abs(z['trace']-1)>g['rdm_trace_abs'] or min(z['spectrum'])<g['rdm_eigenvalue_min']:reasons.append(backend+': RDM positivity/trace')
                wanted=[.5,.5,0.,0.] if pure else [.5,.25,.25,0.]
                if nA==2:
                    wanted=([1.]+[0.]*15) if pure else [((1-eps)/2)]*2+[eps/2]*2+[0.]*12
                delta=sum(abs(x-y) for x,y in zip(z['spectrum'],wanted));errors['analytic_spectrum_L1']=max(errors.get('analytic_spectrum_L1',0.),delta);different|=delta>g['observable_abs']
                purity=sum(x*x for x in wanted)
                for field,delta in [('analytic_purity',abs(z['purity']-purity)),('analytic_S2',abs(z['S2']+math.log(purity)))]:
                    errors[field]=max(errors.get(field,0.),delta);different|=delta>g['observable_abs']
        p=rows['pure_up',-1.];mix=rows['equal_ground_mixture',-1.]
        p_cut={r['nA']:r for r in p['rdms']};m_cut={r['nA']:r for r in mix['rdms']}
        if abs(p_cut[1]['S2']-m_cut[1]['S2'])<g['pure_mixed_entropy_separation_min'] or abs(p['connected_Czz']-mix['connected_Czz'])<g['pure_mixed_connected_separation_min']:
            reasons.append(backend+': pure/mixed fixture does not distinguish prescriptions')
        if rows['canonical_gibbs',s['beta'][-1]]['distance_to_equal_ground']>g['last_beta_trace_distance_max']:reasons.append(backend+': low-T manifold limit unresolved')
    if t['pure_lanczos_residual']>g['eigen_residual_max']:reasons.append('User pure-state Lanczos residual')
    different|=abs(t['pure_lanczos_energy']-3)>g['energy_abs']
    for key in keys:
        a,b=per_backend['QuSpin'][key],per_backend['user_ED_thermal_RDM'][key]
        for field in ('energy','trace','Nmean','spin_z_means','spin_xyz_means','Czz','connected_Czz','distance_to_equal_ground'):
            errors['cross_'+field]=max(errors.get('cross_'+field,0.),error(a[field],b[field]))
        a_cuts={r['nA']:r for r in a['rdms']};b_cuts={r['nA']:r for r in b['rdms']}
        for nA in s['subsystem_sizes']:
            za,zb=a_cuts[nA],b_cuts[nA]
            for field in ('purity','S2'):
                errors['cross_'+field]=max(errors.get('cross_'+field,0.),abs(za[field]-zb[field]))
            errors['cross_spectrum_L1']=max(errors.get('cross_spectrum_L1',0.),sum(abs(x-y) for x,y in zip(za['spectrum'],zb['spectrum'])))
    different|=any(v>g['observable_abs'] for v in errors.values())
    return {'id':'X12','status':'UNRESOLVED' if reasons else ('FAIL' if different else 'PASS'),'max_abs_errors':errors,
        'eligibility_reasons':sorted(set(reasons)),'prescription_count':len(keys),
        'scope':'Canonical total-N3 interacting dimer; pure vector, coherent spin rotation, equal ground mixture, basis-rotated mixture and finite-beta trace each compared with its own oracle. No arbitrary eigenvector gauge is compared.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default=str(OWN/'manifest.toml'))
    for key in ('oracle','target','output'):p.add_argument('--'+key,required=True)
    p.add_argument('--source-root',default=os.environ.get('RENYI_ROOT',str(OWN.parents[1]/'validated_package/HubbardRenyiDQMC')))
    a=p.parse_args();cfg=tomllib.loads(Path(a.manifest).read_text());out=Path(a.output);require(not out.exists(),'Use a new comparison output')
    try:
        o=json.loads(Path(a.oracle).read_text());t=tomllib.loads(Path(a.target).read_text())
        require(finite(o) and finite(t),'Nonfinite receipt')
        require(o['status']==t['status']=='COMPLETED' and o['schema']==t['schema']==cfg['schema'],'Incomplete execution/schema')
        require(o['manifest_sha256']==t['manifest_sha256']==sha(a.manifest),'Manifest receipt')
        require(o['oracle_source_sha256']==sha(OWN/'oracle_quspin.py') and t['adapter_sha256']==sha(OWN/'target.jl'),'Adapter source receipt')
        require(o['versions']=={k:cfg[k+'_version'] for k in ('quspin','numpy','scipy')},'Oracle versions')
        for key in ('julia_version','itensors_version','itensormps_version'):require(t[key]==cfg[key],'Target runtime')
        for rel,key in [('analysis/compute_s2_ed.jl','ed_source_sha256'),('analysis/audit_u4_l4_thermal_ed.jl','thermal_source_sha256'),('analysis/compute_s2_itensor_dmrg.jl','dmrg_source_sha256'),('Project.toml','project_sha256'),('Manifest.toml','environment_manifest_sha256')]:
            require(t[key]==cfg[key]==sha(Path(a.source_root)/rel),'User source/environment receipt')
        rows=[]
        for idea,fn in [('X09',x09),('X12',x12)]:
            try:rows.append(fn(cfg,o[idea],t[idea]))
            except Exception as exc:rows.append({'id':idea,'status':'UNRESOLVED','reason':f'{type(exc).__name__}: {exc}'})
    except Exception as exc:rows=[{'id':i,'status':'UNRESOLVED','reason':f'{type(exc).__name__}: {exc}'} for i in ('X09','X12')]
    result={'schema':cfg['schema'],'manifest_sha256':sha(a.manifest),'rows':rows,'distinct_ideas_passed':sum(x['status']=='PASS' for x in rows),
            'target_sha256':sha(a.target) if Path(a.target).exists() else None,'oracle_sha256':sha(a.oracle) if Path(a.oracle).exists() else None,'comparator_sha256':sha(__file__)}
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2,allow_nan=False))
if __name__=='__main__':main()

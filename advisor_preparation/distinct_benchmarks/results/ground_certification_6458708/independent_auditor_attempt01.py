#!/usr/bin/env python3
"""Independent mixed-status audit of initial frozen X09/X12 receipts; no model run."""
import argparse,hashlib,json,math,os,pathlib,subprocess,sys,tarfile,tempfile,tomllib
BASE=pathlib.Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def need(v,s):
    if not v:raise ValueError(s)
def finite(x):
    if isinstance(x,float):need(math.isfinite(x),'Nonfinite receipt')
    elif isinstance(x,dict):
        for v in x.values():finite(v)
    elif isinstance(x,list):
        for v in x:finite(v)
def diff(a,b):
    if isinstance(a,list) or isinstance(b,list):
        need(isinstance(a,list) and isinstance(b,list) and len(a)==len(b),'Array shape')
        return max((diff(x,y) for x,y in zip(a,b)),default=0.)
    return abs(a-b)
def dot(a,b):return sum(x.conjugate()*y for x,y in zip(a,b))
def apply(H,v):return [sum(row[j]*v[j] for j in range(len(v))) for row in H]
def keyed(rows,fields,keys):
    d={tuple(r[k] for k in fields):r for r in rows};need(len(d)==len(rows) and set(d)==set(keys),'Complete unique grid');return d
def entry(p):return {'path':os.path.relpath(p,BASE),'bytes':p.stat().st_size,'sha256':sha(p)}

def audit(a):
    folder=a.results.resolve();source=a.source.resolve();native=a.renyi_source_root.resolve()
    cfg=tomllib.loads((source/'manifest.toml').read_text());g=cfg['gates']
    t=tomllib.loads((folder/'target.toml').read_text());o=json.loads((folder/'oracle.json').read_text());c=json.loads((folder/'comparison.json').read_text())
    for r in (t,o,c):finite(r)
    archive=BASE/'source_releases/ground-certification-v1-source.tar.gz';freeze=archive.with_name('ground-certification-v1-freeze.json');fr=json.loads(freeze.read_text())
    need(sha(archive)==fr['archive_sha256']=='45298d46a32e3f2d2d2af273a989cd8c10f7cd850e7aea631e8c1e6725dc2dd4','Frozen archive')
    with tarfile.open(archive) as tf:
        ms=[m for m in tf.getmembers() if m.isfile()];need(len(ms)==8 and {m.name for m in ms}==set(fr['files']),'Exact8-file source closure')
        for m in ms:
            raw=tf.extractfile(m).read();need(raw==(source.parent/m.name).read_bytes() and hashlib.sha256(raw).hexdigest()==fr['files'][m.name],'Archive/local byte identity')
    need(o['schema']==t['schema']==cfg['schema'] and o['status']==t['status']=='COMPLETED','Complete producer status')
    need(o['manifest_sha256']==t['manifest_sha256']==c['manifest_sha256']==sha(source/'manifest.toml'),'Manifest identity')
    need(t['adapter_sha256']==sha(source/'target.jl') and o['oracle_source_sha256']==sha(source/'oracle_quspin.py'),'Actual adapter identities')
    need(o['versions']=={k:cfg[k+'_version'] for k in ['quspin','numpy','scipy']},'Independent runtime versions')
    for k in ['julia_version','itensors_version','itensormps_version']:need(t[k]==cfg[k],'Native runtime version')
    native_files=[]
    for rel,k in [('analysis/compute_s2_ed.jl','ed_source_sha256'),('analysis/audit_u4_l4_thermal_ed.jl','thermal_source_sha256'),('analysis/compute_s2_itensor_dmrg.jl','dmrg_source_sha256'),('Project.toml','project_sha256'),('Manifest.toml','environment_manifest_sha256')]:
        p=native/rel;need(sha(p)==t[k]==cfg[k],'Exact user source/environment pin');native_files.append(p)
    need(c['target_sha256']==sha(folder/'target.toml') and c['oracle_sha256']==sha(folder/'oracle.json') and c['comparator_sha256']==sha(source/'compare.py'),'Comparator input identities')
    with tempfile.TemporaryDirectory(prefix='ground_postprocess_') as td:
        out=pathlib.Path(td)/'comparison.json'
        subprocess.run([sys.executable,str(source/'compare.py'),'--manifest',str(source/'manifest.toml'),'--oracle',str(folder/'oracle.json'),'--target',str(folder/'target.toml'),'--source-root',str(native),'--output',str(out)],check=True,capture_output=True,text=True)
        need(out.read_bytes()==(folder/'comparison.json').read_bytes()==(folder/'comparison_independent_replay_attempt01.json').read_bytes(),'Byte-identical comparator replay')
    reported={r['id']:r['status'] for r in c['rows']}
    need(len(c['rows'])==2 and set(reported)=={'X09','X12'} and c['distinct_ideas_passed']==sum(v=='PASS' for v in reported.values()),'Faithful complete mixed scientific status')

    # Recompute ground fidelities, moments and eigen-residuals from all36 complex
    # coefficients, never reading the comparator's derived certificate values.
    s=cfg['X09'];oo=o['X09'];tt=t['X09'];expected={(u,d) for u in range(16) for d in range(16) if u.bit_count()==d.bit_count()==2}
    oc=[tuple(v) for v in oo['coordinates']];tc=[tuple(v) for v in tt['coordinates']]
    need(len(oc)==len(tc)==36 and set(oc)==set(tc)==expected,'All36 physical fixed-spin coordinates')
    idx={v:i for i,v in enumerate(oc)};ix=[idx[v] for v in tc]
    H=[[oo['Hamiltonian'][i][j] for j in ix] for i in ix];P=[[oo['ground_projector'][i][j] for j in ix] for i in ix]
    need(oo['ground_dimension']==1 and len(oo['eigenvalues'])==36 and oo['gap']>g['ground_gap_min'],'Externally unique gapped ground projector')
    matrix_errors={k:diff(tt[k],H) for k in ['Hamiltonian','MPO_qn_Hamiltonian','MPO_dense_Hamiltonian']}
    matrix_errors.update({k:abs(tt[k]) for k in ['MPO_qn_max_imaginary','MPO_dense_max_imaginary']})
    need(max(matrix_errors.values())<=g['matrix_abs'],'Whole native Hamiltonian/MPO comparison')
    for field in ['source_eigen_residual']:need(tt[field]<=g['eigen_residual_max'],'Actual ED eigensystem residual')
    for field in ['eigen_residual','projector_idempotency']:need(oo[field]<=g['eigen_residual_max'],'External projector eligibility')
    source_rows=keyed(tt['rows'],['label'],[(x,) for x in s['starts']+s['challenge_states']]);certificates=[];unresolved=[];candidate_failure=False
    for (label,),r in source_rows.items():
        need(len(r['coeff_real'])==len(r['coeff_imag'])==36,'Complete MPS sector amplitude export')
        vector=[complex(x,y) for x,y in zip(r['coeff_real'],r['coeff_imag'])];norm=dot(vector,vector).real
        need(abs(norm-1)<=g['norm_abs'] and abs(norm-r['norm'])<=g['norm_abs'] and abs(r['outside_sector_weight'])<=g['sector_abs'],'Full normalized sector support')
        need(abs(r['Nup']-2)<=g['sector_abs'] and abs(r['Ndown']-2)<=g['sector_abs'],'Actual MPS quantum numbers')
        Hv=apply(H,vector);energy=dot(vector,Hv).real/norm;second=dot(Hv,Hv).real/norm
        variance=sum(abs(h-energy*v)**2 for h,v in zip(Hv,vector))/norm
        fidelity=dot(vector,apply(P,vector)).real/norm;excess=energy-oo['ground_energy']
        need(max(abs(energy-r['energy']),abs(second-r['energy_second_moment']),abs(variance-r['variance']))<=g['energy_abs'],'Independent Rayleigh moment replay')
        lowvar=g['variance_min']<=variance<=g['variance_max'];isground=lowvar and abs(excess)<=g['energy_abs'] and fidelity>=g['ground_fidelity_min']
        if label in s['starts']:
            pos=s['starts'].index(label);need(r['kind']=='native_dmrg' and r['seed']==s['seeds'][pos] and len(r['truncation_errors'])==len(s['maxdim_schedule']),'All actual seeded native sweeps')
            need(abs(r['reported_energy']-energy)<=g['energy_abs'],'Actual optimizer reported/matrix energy consistency')
            if not lowvar:unresolved.append(label+': native DMRG variance unresolved')
            candidate_failure|=not isground
        else:
            residual=math.sqrt(sum(abs(x-r['source_eigenvalue']*v)**2 for x,v in zip(Hv,vector))/norm)
            need(r['kind']=='manufactured_excited_mps' and lowvar and not isground and fidelity<=g['challenge_fidelity_max'] and excess>=oo['gap']/2,'Manufactured low-variance excited state rejected')
            need(r['source_eigen_residual']<=g['eigen_residual_max'] and residual<=g['eigen_residual_max'],'Manufactured source/MPS eigenstate identity')
        certificates.append({'label':label,'kind':r['kind'],'norm':norm,'external_energy':energy,'external_energy_excess':excess,'external_ground_fidelity':fidelity,'independent_variance':variance,'low_variance':lowvar,'ground_certificate':isground})

    x09_status='UNRESOLVED' if unresolved else 'FAIL' if candidate_failure else 'PASS'
    need(reported['X09']==x09_status,'Independent X09 status agrees with preserved comparator')
    # Recompute analytic dimer ensembles and spectra rather than treating a
    # comparison PASS or a single ground eigenvector as the mixed-state oracle.
    s=cfg['X12'];theta=math.pi/s['rotation_theta_pi_denominator'];phi=math.pi/s['rotation_phi_pi_denominator'];metrics={}
    wanted=[(name,-1.) for name in s['prescriptions']]+[('canonical_gibbs',b) for b in s['beta']]
    for backend,data in [('QuSpin',o['X12']),('user_ED_thermal_RDM',t['X12'])]:
        need(data['canonical_dimension']==4 and data['ground_dimension']==2 and diff(data['canonical_eigenvalues'],[3.,3.,5.,5.])<=g['energy_abs'],'Full totalN3 canonical doublets')
        need(abs(data['gap']-2)<=g['energy_abs'] and data['eigen_residual']<=g['eigen_residual_max'] and data['mixture_rotation_error']<=g['matrix_abs'],'Ground-manifold invariance and full spectrum')
        need(abs(data['spin_lowering_norm']-1)<=g['observable_abs'],'Fixed spin-doublet lowering phase')
        rr=keyed(data['rows'],['label','beta'],wanted);maxima={}
        def check(name,value,bound=g['observable_abs']):
            maxima[name]=max(maxima.get(name,0.),value);need(value<=bound,'Analytic ensemble '+backend+' '+name)
        for (name,beta),r in rr.items():
            pure=name.startswith('pure_');eps=1/(1+math.exp(2*beta)) if beta>=0 else 0.
            spin=[0.,0.,.25] if name=='pure_up' else [math.sin(2*theta)*math.cos(phi)/4,math.sin(2*theta)*math.sin(phi)/4,math.cos(2*theta)/4] if name=='pure_rotated' else [0.,0.,0.]
            for field,value in [('energy',3+2*eps),('trace',1.),('Nmean',3.),('Czz',0.),('connected_Czz',-spin[2]**2),('distance_to_equal_ground',1. if pure else 2*eps)]:check(field,abs(r[field]-value))
            check('spin_xyz',diff(r['spin_xyz_means'],[spin,spin]));check('spin_z',diff(r['spin_z_means'],[spin[2]]*2))
            if beta>=0:check('excited_weight',abs(r['excited_weight']-eps))
            cuts=keyed(r['rdms'],['nA'],[(1,),(2,)])
            for (nA,),rho in cuts.items():
                spectrum=rho['spectrum'];need(len(spectrum)==4**nA and spectrum==sorted(spectrum,reverse=True) and min(spectrum)>=g['rdm_eigenvalue_min'],'Complete positive RDM spectrum')
                expected=([.5,.5,0.,0.] if pure else [.5,.25,.25,0.]) if nA==1 else ([1.]+[0.]*15 if pure else [(1-eps)/2]*2+[eps/2]*2+[0.]*12)
                check('cut'+str(nA)+'_spectrum_L1',sum(abs(x-y) for x,y in zip(spectrum,expected)))
                purity=sum(x*x for x in spectrum);exact_purity=sum(x*x for x in expected)
                check('RDM_trace',max(abs(sum(spectrum)-1),abs(rho['trace']-1)),g['rdm_trace_abs'])
                check('purity_vs_spectrum',abs(rho['purity']-purity));check('purity_vs_analytic',abs(purity-exact_purity));check('S2_vs_analytic',abs(rho['S2']+math.log(exact_purity)))
        pure=rr['pure_up',-1.];mix=rr['equal_ground_mixture',-1.]
        purecut={r['nA']:r for r in pure['rdms']};mixcut={r['nA']:r for r in mix['rdms']}
        separation={'S2':abs(purecut[1]['S2']-mixcut[1]['S2']),'connected_Czz':abs(pure['connected_Czz']-mix['connected_Czz'])}
        need(separation['S2']>=g['pure_mixed_entropy_separation_min'] and separation['connected_Czz']>=g['pure_mixed_connected_separation_min'],'Failure-sensitive ensemble separation')
        need(rr['canonical_gibbs',s['beta'][-1]]['distance_to_equal_ground']<=g['last_beta_trace_distance_max'],'Cooling to full equal-ground mixture')
        metrics[backend]={'prescriptions':10,'RDMs':20,'max_errors':maxima,'pure_mixed_separation':separation}
    need(t['X12']['pure_lanczos_residual']<=g['eigen_residual_max'] and abs(t['X12']['pure_lanczos_energy']-3)<=g['energy_abs'],'Native fixed-spin Lanczos endpoint')
    need(reported['X12']=='PASS','Independent X12 result matches preserved comparator')
    inputs=[folder/n for n in ['oracle.json','target.toml','comparison.json','comparison_independent_replay_attempt01.json']]+[archive,freeze]+[source.parent/rel for rel in fr['files']]+native_files
    metadata=folder/'remote_transfer_metadata.json'
    if metadata.exists():
        remote=json.loads(metadata.read_text());need(remote['job_id']==str(a.job_id),'Remote scheduler/transfer job')
        for rel,v in remote['files'].items():
            p=BASE/rel;need(sha(p)==v['sha256'] and p.stat().st_size==v['bytes'],'All remote/local text transfer bytes');inputs.append(p)
        need(f'{a.job_id}|COMPLETED|0:0|' in remote['sacct'],'Completed job exit');inputs.append(metadata)
    else:raise ValueError('Missing exact remote transfer metadata')
    records=[entry(p) for p in dict.fromkeys(inputs)]
    return {'status':'PASS','audit_scope':'Successful audit of a mixed scientific outcome; top-level PASS does not certify every benchmark','scientific_completions':sum(v=='PASS' for v in reported.values()),'catalogue_ids':['X09','X12'],'source_archive_sha256':sha(archive),'independent_comparator_replay':'BYTE_IDENTICAL','source_closure_members':8,'native_source_environment_pins':5,'benchmarks':{'X09':{'status':x09_status,'eligibility_reasons':unresolved,'native_optimizer_starts':4,'manufactured_challenges':2,'matrix_errors':matrix_errors,'certificates':certificates},'X12':{'status':'PASS','backends':metrics}},'inputs':records,'input_sha256':{r['path']:r['sha256'] for r in records},'limits':'Manufactured excited MPS states test certificate failure sensitivity, not natural optimizer failure. This audit executes only scalar/vector arithmetic and frozen postprocessing, no model, eigensolver, native optimizer or Monte Carlo.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=pathlib.Path,required=True);p.add_argument('--source',type=pathlib.Path,default=BASE/'ground_certification');p.add_argument('--renyi-source-root',type=pathlib.Path,default=BASE.parent/'validated_package/HubbardRenyiDQMC');p.add_argument('--job-id',type=int,default=6458708);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args();need(not a.output.exists(),'Use a new immutable audit output')
    try:r=audit(a)
    except Exception as exc:r={'status':'UNRESOLVED','scientific_completions':0,'reason':f'{type(exc).__name__}: {exc}'}
    r.update(schema='ground_certification_independent_receipt_audit_v1',auditor_sha256=sha(__file__),no_physics_executed=True)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':r['status'],'scientific_completions':r['scientific_completions'],'reason':r.get('reason')},indent=2))
    if r['status']!='PASS':raise SystemExit(1)
if __name__=='__main__':main()

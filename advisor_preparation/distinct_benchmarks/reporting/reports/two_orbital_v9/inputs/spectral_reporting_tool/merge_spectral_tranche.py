#!/usr/bin/env python3
"""Read-only scientific postprocessing: immutable32v7 plus spectral P08/P09/P10.

No physics, jobs or edits to frozen sources. Python3.11+ standard library only.
Every output directory is new and retains all base evidence and authored source.
"""
import argparse,copy,csv,hashlib,json,math,shutil,subprocess,sys,tarfile,tomllib
from datetime import datetime,timezone
from pathlib import Path,PurePosixPath

BASE=Path(__file__).resolve().parents[1]
PREP=BASE.parent
BASE32_SHA='219abfcee07899eba7f6c6501d463b2afd32816999cb7944c58a149ddd5fabe5'
BASE32_MANIFEST_SHA='1d434cc1a38577eb9500a3818d5edf67fcf00ed83bd3ca18cfec7913493c8ee8'
SPECTRAL_ARCHIVE_SHA='3fb6ce6accac61140ee21584dddbb79418b121f7d62ac9391b31996f1cf2239b'
SPECTRAL_FREEZE_SHA='8f7ad7607c618d4ca56684ee81c3dcdee825e8901c105fe5ed53ae36146686ff'
SPECTRAL_AUDIT_SHA='1ccbd313086027ab974f37255908dcb8da64bd5bd3cd78e292abb2dfe3bd13ae'
SPECTRAL_AUDITOR_SHA='30088a735570b2f68a263652a1f68fa27c2bb07e569420f418d758def2944dbe'
IDS=('P08','P09','P10')
RAW_SHA={'oracle.json':'31c75279dacfb753ce40723f8f6c11a8060a6a4fd08da4aafba360fde311a5b8',
         'target.toml':'308f528021a63b11b2f11fbf79cde70eb5efa96814f21f058b6a171392bae7df',
         'comparison.json':'d55b83e6b809a952e41de8ea4ac8e6afebac0be6f3f32c9856914a33f01032af'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read_json(p):
    return json.loads(Path(p).read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON '+x)))
def dump(p,data):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def relative(value):
    p=PurePosixPath(value)
    require(not p.is_absolute() and '..' not in p.parts and str(p) not in ('','.'),'Unsafe relative path')
    return Path(*p.parts)
def finite(x):
    if isinstance(x,dict):return all(finite(v) for v in x.values())
    if isinstance(x,list):return all(finite(v) for v in x)
    return not isinstance(x,float) or math.isfinite(x)
def row_hash(row):return hashlib.sha256(json.dumps(row,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def counts(rows):
    out={s:sum(r['status']==s for r in rows) for s in ('PASS','FAIL','UNRESOLVED','NOT_RUN')}
    require(sum(out.values())==100,'Unrecognized or incomplete statuses')
    return out
def replay_equal(a,b):
    diffs=[]
    def visit(x,y,path):
        require(type(x) is type(y),'Replay type differs:'+path)
        if isinstance(x,dict):
            require(x.keys()==y.keys(),'Replay keys differ:'+path)
            for k in x:visit(x[k],y[k],path+'/'+k)
        elif isinstance(x,list):
            require(len(x)==len(y),'Replay length differs:'+path)
            for i,(u,v) in enumerate(zip(x,y)):visit(u,v,path+'/'+str(i))
        elif isinstance(x,float):
            error=abs(x-y);require(math.isfinite(x) and math.isfinite(y) and error<=64*sys.float_info.epsilon*max(1.,abs(x),abs(y)),'Replay float differs beyond roundoff:'+path)
            if error:diffs.append({'path':path,'absolute_difference':error,'reported':x,'replayed':y})
        else:require(x==y,'Replay nonfloat differs:'+path)
    visit(a,b,'')
    return {'nonfloat_values':'exact','float_allowance':'64*binary64_epsilon*max(1,abs(a),abs(b)); frozen scientific gates unchanged','roundoff_differences':diffs,'maximum_absolute_roundoff_difference':max((r['absolute_difference'] for r in diffs),default=0.)}

class Evidence:
    def __init__(self,out,manifest):self.out=out;self.manifest=manifest
    def add(self,role,p,expected=None):
        p=Path(p).resolve();require(p.is_file() and not p.is_symlink(),'Missing/nonregular evidence:'+str(p))
        digest=sha(p);require(expected is None or digest==expected,'Evidence SHA differs:'+role)
        require(role not in self.manifest,'Duplicate evidence role:'+role)
        dest=Path('inputs')/role/p.name;target=self.out/dest
        require(not target.exists(),'Duplicate evidence destination:'+str(dest));target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
        try:original=p.relative_to(PREP.resolve()).as_posix()
        except ValueError:original='external/'+role+'/'+p.name
        self.manifest[role]={'original_path':original,'evidence_copy':dest.as_posix(),'sha256':digest,'bytes':target.stat().st_size}
        return target

def copy_base32(path,out):
    path=Path(path).resolve();ledger=path/'ledger_100.json';manifest=path/'input_manifest.json'
    require(sha(ledger)==BASE32_SHA and sha(manifest)==BASE32_MANIFEST_SHA,'Immutable32v7 identity differs')
    base=read_json(ledger);inputs=read_json(manifest)
    require(base['input_manifest']==inputs and base['counts']=={'PASS':32,'FAIL':0,'UNRESOLVED':0,'NOT_RUN':68},'Base32 manifest/counts differ')
    expected={f'{prefix}{i:02d}' for prefix in ('A','J','P','X') for i in range(1,26)}
    require(len(base['rows'])==100 and {r['id'] for r in base['rows']}==expected,'Base catalogue is not exactly100uniqueIDs')
    require(counts(base['rows'])==base['counts'],'Base row counts differ')
    for role,r in inputs.items():
        file=path/relative(r['evidence_copy'])
        require(file.is_file() and not file.is_symlink() and sha(file)==r['sha256'] and file.stat().st_size==r['bytes'],'Base evidence differs:'+role)
    for rel,digest in base['analysis_receipts'].items():require(sha(path/relative(rel))==digest,'Base analysis receipt differs:'+rel)
    closure=[]
    for child in sorted(path.rglob('*')):
        require(not child.is_symlink(),'Base report contains symlink')
        if not child.is_file():continue
        rel=child.relative_to(path)
        # Keep all evidence/analysis/frozen source paths stable for original rows.
        dest=rel if len(rel.parts)>1 else Path('inputs/base32_metadata')/rel
        target=out/dest;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(child,target)
        closure.append({'base_relative_path':rel.as_posix(),'evidence_copy':dest.as_posix(),'sha256':sha(child),'bytes':child.stat().st_size})
    for name in sorted(p.name for p in path.iterdir() if p.is_file()):
        role='base32_metadata_'+name.replace('.','_');p=out/'inputs/base32_metadata'/name
        inputs[role]={'original_path':str((path/name).relative_to(PREP.resolve())),'evidence_copy':str(p.relative_to(out)),'sha256':sha(p),'bytes':p.stat().st_size}
    dump(out/'analysis/base32_file_closure.json',closure)
    return base,inputs,closure

def spectral_release(ev):
    freeze=read_json(ev.add('spectral_freeze',BASE/'source_releases/spectral-dynamics-v1-freeze.json',SPECTRAL_FREEZE_SHA))
    archive=ev.add('spectral_archive',BASE/'source_releases/spectral-dynamics-v1-source.tar.gz',SPECTRAL_ARCHIVE_SHA)
    require(freeze['release']=='spectral-dynamics-v1' and freeze['archive_sha256']==SPECTRAL_ARCHIVE_SHA and freeze['distinct_ids']==list(IDS),'Spectral release identity/IDs differ')
    seen={};root=ev.out/'frozen_sources/spectral'
    with tarfile.open(archive) as tar:
        for item in tar:
            if item.isdir():continue
            require(item.isfile(),'Nonregular archive member')
            rel=relative(item.name)
            require(rel.parts[0]=='spectral_dynamics','Unexpected archive root')
            key=rel.as_posix()
            require(key in freeze['files'] and key not in seen,'Extra/duplicate frozen file')
            data=tar.extractfile(item).read();digest=hashlib.sha256(data).hexdigest()
            require(digest==freeze['files'][key],'Frozen spectral source differs:'+key)
            dest=root/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data);seen[key]=digest
    require(seen==freeze['files'],'Spectral archive closure incomplete')
    return root/'spectral_dynamics',freeze

def audit_raw(cfg,o,t):
    require(finite(o) and finite(t),'Nonfinite spectral receipt')
    require(o['schema']==t['schema']=='spectral_dynamics_v1' and o['status']==t['status']=='COMPLETED','Incomplete spectral execution')
    require(cfg['catalogue_ids']==list(IDS) and cfg['L']==4 and cfg['t']==1. and cfg['U']==4. and cfg['boundary']=='pbc','Frozen fixture differs')
    require(o['versions']=={k:cfg[k+'_version'] for k in ('quspin','numpy','scipy')} and t['julia_version']==cfg['julia_version'],'Runtime version differs')
    times=[cfg['P08']['time_start']+(cfg['P08']['time_stop']-cfg['P08']['time_start'])*i/(cfg['P08']['time_count']-1) for i in range(cfg['P08']['time_count'])]
    spectral={(s,k) for s in cfg['P09']['spin'] for k in cfg['P09']['momentum_indices']}
    thermal={(b,k,f) for b in cfg['P10']['beta'] for k in cfg['P10']['momentum_indices'] for f in cfg['P10']['tau_fractions']}
    for data in (o,t):
        require([r['time'] for r in data['P08']['rows']]==times,'Partial/duplicate quench time grid')
        require(all(len(r['density'])==cfg['L'] for r in data['P08']['rows']),'Partial site density profile')
        rows=data['P09']['rows'];require(len(rows)==len(spectral) and {(r['spin'],r['momentum_index']) for r in rows}==spectral,'Partial/duplicate spectral grid')
        for r in rows:
            for branch in ('addition','removal'):
                require(sum(z['multiplicity'] for z in r[branch])==24,'IncompleteN±1 energy clusters')
        rows=data['P10']['rows'];require(len(rows)==len(thermal) and {(r['beta'],r['momentum_index'],r['tau_fraction']) for r in rows}==thermal,'Partial/duplicate thermal response grid')
        require(len(data['P10']['thermal'])==len(cfg['P10']['beta']) and {r['beta'] for r in data['P10']['thermal']}==set(cfg['P10']['beta']),'Incomplete thermal ensemble metadata')
    return {'P08_time_points':len(times),'P08_site_density_components':len(times)*cfg['L'],'P09_spin_momentum_channels':len(spectral),'P09_all_N_plus_minus_eigenstates_per_branch':24,'P10_beta_q_tau_points':len(thermal),'all_fixture_grids_complete':True}

def check_independent(audit,rawpaths,reported,replayed):
    require(audit['schema']=='spectral_receipt_independent_audit_v1' and audit['status']=='PASS','Independent scientific receipt audit notPASS')
    # The independent auditor must bind every raw result and both complete comparisons.
    require(audit['catalogue_ids']==list(IDS) and audit['source_archive_sha256']==SPECTRAL_ARCHIVE_SHA,'Independent audit scope differs')
    require(audit['auditor_sha256']==SPECTRAL_AUDITOR_SHA,'Independent auditor source differs')
    require('input_sha256' in audit,'Independent audit lacks input hash bindings')
    hashes=audit['input_sha256']
    for label,path in rawpaths.items():
        matches=[v for k,v in hashes.items() if Path(k).name==Path(path).name or k==label]
        require(sha(path) in matches,'Independent audit does not bind:'+label)
    replay_equal(reported,replayed)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base32',type=Path,default=BASE/'reporting/reports/full32_public_replay_v7')
    p.add_argument('--spectral-dir',type=Path,default=BASE/'results/spectral_dynamics_6458468')
    p.add_argument('--independent-audit',type=Path,required=True)
    p.add_argument('--independent-replay',type=Path,required=True)
    p.add_argument('--independent-auditor',type=Path,default=BASE/'reporting/audit_spectral_receipts.py')
    p.add_argument('--source-root',type=Path,default=PREP/'validated_package/HubbardRenyiDQMC')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();require(not a.output.exists(),'Choose a new output directory');a.output.mkdir(parents=True)
    try:
        base,manifest,closure=copy_base32(a.base32,a.output);ev=Evidence(a.output,manifest)
        ev.add('spectral_reporting_tool',__file__)
        source,freeze=spectral_release(ev);cfg=tomllib.loads((source/'manifest.toml').read_text())
        raw={name:ev.add('spectral_'+name.replace('.','_'),a.spectral_dir/name,digest) for name,digest in RAW_SHA.items()}
        o=read_json(raw['oracle.json']);t=tomllib.loads(raw['target.toml'].read_text());reported=read_json(raw['comparison.json'])
        manifest_sha=sha(source/'manifest.toml')
        require(o['manifest_sha256']==t['manifest_sha256']==reported['manifest_sha256']==manifest_sha,'Spectral manifest receipt differs')
        require(o['oracle_source_sha256']==sha(source/'oracle_quspin.py') and t['adapter_sha256']==sha(source/'target_ed.jl'),'Actual adapter source differs')
        for field,path in [('thermal_source_sha256','analysis/audit_u4_l4_thermal_ed.jl'),('ed_source_sha256','analysis/compute_s2_ed.jl')]:
            require(t[field]==cfg[field],'User source receipt differs');ev.add('spectral_user_source_'+Path(path).name,a.source_root/path,cfg[field])
        completeness=audit_raw(cfg,o,t)
        recorded_replay=read_json(ev.add('spectral_independent_replay',a.independent_replay,RAW_SHA['comparison.json']))
        independent=read_json(ev.add('spectral_independent_audit',a.independent_audit,SPECTRAL_AUDIT_SHA))
        ev.add('spectral_independent_auditor_source',a.independent_auditor,SPECTRAL_AUDITOR_SHA)
        check_independent(independent,{**raw,'independent_replay':a.independent_replay},reported,recorded_replay)
        for name in ('oracle.resources.txt','target.resources.txt','retrieval_sha256.txt'):
            ev.add('spectral_'+name.replace('.','_'),a.spectral_dir/name)
        ev.add('spectral_job_submission',BASE/'submissions/spectral_dynamics_6458468.json')
        ev.add('spectral_job_log',BASE/'logs/spectral_dynamics_6458468.log')
        replay=a.output/'analysis/spectral_comparator_replay.json'
        cmd=[sys.executable,str((source/'compare.py').resolve()),'--manifest',str((source/'manifest.toml').resolve()),'--oracle',str(raw['oracle.json'].resolve()),'--target',str(raw['target.toml'].resolve()),'--output',str(replay.resolve())]
        proc=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
        dump(a.output/'analysis/spectral_replay_process.json',{'command':cmd,'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr,'physics_executed':False})
        require(proc.returncode==0 and replay.is_file(),'Frozen spectral comparator replay failed')
        recomputed=read_json(replay);agreement=replay_equal(reported,recomputed);replay_equal(recorded_replay,recomputed)
        require(recomputed['schema']==cfg['schema'] and recomputed['manifest_sha256']==manifest_sha,'Comparator identity differs')
        rmap={r['id']:r for r in recomputed['rows']}
        require(len(recomputed['rows'])==3 and set(rmap)==set(IDS) and recomputed['distinct_ideas_passed']==3,'Incomplete spectral comparator')
        require(all(r['status']=='PASS' and r['eligibility_reasons']==[] for r in rmap.values()),'A spectral idea failed eligibility or frozen scientific gates; do not publish partial35ledger')
        rows=copy.deepcopy(base['rows']);before={r['id']:r for r in base['rows']};preserved={}
        for row in rows:
            if row['id'] not in IDS:preserved[row['id']]=row_hash(row);require(row==before[row['id']],'Original row changed');continue
            require(row['status']=='NOT_RUN','Spectral ID would replace prior scientific outcome')
            idea=row['id'];r=rmap[idea];gate=cfg['gates'];metric=[]
            for name,error in r['max_abs_errors'].items():
                tolerance=gate['pole_abs'] if name=='pole' else gate['residue_abs'] if name=='residue' else gate['cross_observable_abs']
                require(error<=tolerance,'Comparator metric exceeds declared bound')
                metric.append({'metric':name,'maximum_error':error,'absolute_tolerance':tolerance,'all_within_tolerance':True})
            units={'P08':(1,completeness['P08_time_points'],'complete real-time grid points'),'P09':(1,completeness['P09_spin_momentum_channels'],'complete spin/momentum spectra'),'P10':(len(cfg['P10']['beta']),completeness['P10_beta_q_tau_points'],'complete beta/momentum/tau points')}
            f,n,u=units[idea]
            row.update(status='PASS',reported_status='PASS',reason='',fixtures=f,comparisons=n,comparison_unit=u,metrics=metric,
                scope={'method':'Actual user ED/thermal APIs with new spectral/dynamics observable adapter versus independent QuSpin',
                       'existing_user_functions':{'P08':['build_space','apply_hamiltonian!'],'P09':['build_space','apply_hamiltonian!','ground_state_lanczos'],'P10':['full_fock_matrices','exact_thermal_density']}[idea],
                       'new_adapter':{'P08':'Fixed-N unitary CDW quench and local density/doublon measurement','P09':'Physical fermion operators, fullN±1 spin sectors and clustered Lehmann spectra','P10':'Connected Sz(q,tau) thermal Lehmann response'}[idea],
                       'production_DQMC_sampling_exercised':False,'production_DQMC_dynamics_exercised':False,'production_DMRG_exercised':False,
                       'fixture':'L4,t1,U4,PBC; exact frozen ensembles/normalizations in spectral manifest'},
                input_role='spectral_comparison_json',audit_issues=[],executed_fixture_scope=r['scope'])
        require(all(row==before[row['id']] for row in rows if row['id'] not in IDS),'Base row semantics changed')
        require(sum(before[r['id']]['status']=='PASS' and r==before[r['id']] for r in rows)==32,'Original32PASSrows not preserved')
        tally=counts(rows);require(tally=={'PASS':35,'FAIL':0,'UNRESOLVED':0,'NOT_RUN':65},'Unexpected35case counts')
        category={c:{s:sum(r['category']==c and r['status']==s for r in rows) for s in tally} for c in sorted({r['category'] for r in rows})}
        analysis={p.relative_to(a.output).as_posix():sha(p) for p in sorted((a.output/'analysis').rglob('*')) if p.is_file()}
        ledger={'generated_utc':datetime.now(timezone.utc).isoformat(),'tool_sha256':sha(__file__),'base32_ledger_sha256':BASE32_SHA,
                'counts':tally,'category_counts':category,'counting_rule':base['counting_rule'],
                'scope_limit':'New P08/P09/P10 additions exercise actual user ED/thermal APIs with explicit new observable adapters; no native DQMC dynamics or DMRG pass is inferred. Original32rows retain their exact scopes, including X06 sampling calibration.',
                'base32_preservation':{'original32PASSrows_unchanged':True,'all_other97row_objects_unchanged':True,'canonical_row_sha256':preserved,'all_base_files_copied':len(closure)},
                'spectral_audit':{'source_release':freeze['release'],'source_archive_sha256':SPECTRAL_ARCHIVE_SHA,'completeness':completeness,'comparator_replay':agreement,'independent_receipt_audit_status':independent['status'],'no_partial_idea_pass':True},
                'input_manifest':ev.manifest,'analysis_receipts':analysis,'rows':rows}
        dump(a.output/'input_manifest.json',ev.manifest);dump(a.output/'ledger_100.json',ledger)
        with (a.output/'ledger_100.csv').open('w',newline='')as f:
            fields=['id','category','title','status','reported_status','fixtures','comparisons','comparison_unit','reason','metrics'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
            for row in rows:w.writerow({k:json.dumps(row[k],sort_keys=True) if k=='metrics' else row.get(k,'') for k in fields})
        (a.output/'README.md').write_text('# 35 distinct externally anchored benchmarks passed\n\nExactly 100 catalogue IDs: 35 PASS, 0 FAIL, 0 UNRESOLVED, 65 NOT_RUN. Parameters and matrix, spectral or time-grid entries are not additional ideas.\n\nThe immutable 32v7 ledger and every base file are preserved; all 97 other row objects, including the 32 PASS rows, are unchanged. Only P08/P09/P10 are added after frozen source closure, complete-grid eligibility, actual-user source/version checks, the original and independent comparison receipts, and a fresh local replay of the exact frozen comparator. All authored source and paired raw results are copied or extracted in this report.\n\nThe additions concern a four-site Hubbard CDW quench, complete fermionic Lehmann spectra, and connected imaginary-time spin response. They use new observable adapters around actual user ED/thermal APIs and independently constructed QuSpin operators. They do not establish native DQMC dynamics or larger-system accuracy.\n\nSee [ledger_100.json](ledger_100.json), [input_manifest.json](input_manifest.json), [base-file closure](analysis/base32_file_closure.json) and [replayed comparator](analysis/spectral_comparator_replay.json). This tool performs no scientific simulation or cluster submission.\n')
        print(json.dumps({'counts':tally,'base32_preserved':True,'spectral_replay':'PASS','output':str(a.output)},indent=2))
    except Exception as exc:
        dump(a.output/'ROLLUP_FAILURE.json',{'status':'REJECTED_NO_NEW_LEDGER','reason':f'{type(exc).__name__}: {exc}','tool_sha256':sha(__file__),'physics_executed':False})
        raise
if __name__=='__main__':main()

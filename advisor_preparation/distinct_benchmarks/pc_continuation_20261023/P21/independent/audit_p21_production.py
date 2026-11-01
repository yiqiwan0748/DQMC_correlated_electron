#!/usr/bin/env python3
"""Independent local P21 production audit; no sampling or remote operations.

Uses the byte-pinned independent pilot paired-ratio math and the byte-pinned
independent X01 rank/ESS implementation. Native text parsing is the pilot's
independent grammar adapted explicitly to68bins and1024cycles, not the new
production parser. Production comparator replay is supplemental evidence only.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path,PurePosixPath
import re
import statistics
import subprocess
import sys
import tarfile
import tempfile

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
PILOT_TOOL=ROOT/'audit_p21_sampled_pilot.py'
PILOT_SHA='d7398ec452b437bede06985a0075bc584414c1d5ad1808d661284dc5f41faa93'
MIX_TOOL=ROOT/'audit_continuum_stage1.py'
MIX_SHA='6d68b99f06eeea2834f72e386e7b91fbe63bf8aa62f265c95947361eb761910c'
FREEZE_SHA='00c5eced240553bac20fbe3e41ec114a59f3a59916f470b2d57b98f7fb4e82ba'
ARCHIVE_SHA='ef35cce4e7670ef3c642946fac345055144c987fc6059887022bfc3be5da12a1'
REVIEW_SHA='0d0409e2806c7d8232479fecc1470b7169382147440991a2400d7e6e7d7a7afe'
JOB='6461010'
RELEASE='p21_sampled_production_v1'
CUTOFF=1792755300
CRITICAL=11.803638771292341
SCHEMA='P21_PRODUCTION_INDEPENDENT_AUDIT_V1'


def compact_audit(mixing):
    """Keep every failure, full input hashes and grid; aggregate repetitive PASSes."""
    prefixes=('independent paired ratio/blocking ','exact blocking decision ','exact component coordinate ',
              'independent rank/folded/Geyer diagnostics ','backend eligibility decision ',
              'independent Student7 pair interval ','exact real component eligibility ',
              'independent complex rectangle disk gate ','immutable raw user bin ',
              'structural source-bound diagnostic ','bound input unchanged at audit end ',
              'original transfer ')
    class Compact(mixing.Audit):
        def __init__(self):
            super().__init__();self.atomic_check_count=0;self.groups={}
        def check(self,name,ok,detail=None):
            self.atomic_check_count+=1
            prefix=next((p for p in prefixes if name.startswith(p)),None)
            if prefix:
                group=self.groups.setdefault(prefix,dict(checked=0,failures=0));group['checked']+=1;group['failures']+=int(not ok)
                if ok:return True
            return super().check(name,ok,detail)
        def finalize(self):
            for prefix,group in self.groups.items():
                super().check('aggregate complete grid: '+prefix,group['failures']==0,group)
    return Compact()


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(ok,reason):
    if not ok:raise ValueError(reason)


def source_module(path,want,name):
    require(sha(path)==want,'Approved independent math source changed: '+str(path))
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def eq(audit,pilot,label,left,right):
    try:pilot.equal(left,right,label)
    except ValueError as e:return audit.check(label,False,str(e))
    return audit.check(label,True)


def native_vectors(path,width):
    """Independent pilot vector grammar, with the declared68-bin budget."""
    rows=[];pattern=re.compile(r'\(\s*([^(),]+)\s*,\s*([^(),]+)\s*\)')
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        matches=list(pattern.finditer(line));require(len(matches)==width,'Exactly full native vector width')
        require(line[:matches[0].start()].strip()==str(width+1),'Native width marker')
        require(all(not line[x.end():y.start()].strip() for x,y in zip(matches,matches[1:])),'No unparsed native vector tokens')
        def number(text):
            x=float(text.strip().replace('D','E').replace('d','e'));require(math.isfinite(x),'Finite native scalar');return x
        values=[complex(number(m.group(1)),number(m.group(2))) for m in matches]
        rows.append((values,number(line[matches[-1].end():])))
    require(len(rows)==68,'Exactly68 original native bins including4burn')
    return rows


def protocol(audit,pilot,mixing,a):
    source=a.release_dir.resolve();freeze=audit.load(a.freeze)
    audit.check('approved production external freeze',sha(a.freeze)==FREEZE_SHA and freeze['archive_sha256']==ARCHIVE_SHA and freeze['schema']=='P21_SAMPLED_PRODUCTION_FREEZE_V1' and freeze['phase']=='PRODUCTION_SINGLE_LOOK' and len(freeze['files'])==52)
    audit.check('administrative freeze exact bytes',audit.file(source/'source_freeze.json').read_bytes()==Path(a.freeze).read_bytes())
    for name,want in freeze['files'].items():
        p=(source/name).resolve();require(p.is_relative_to(source),'Safe frozen source path')
        audit.check('frozen source '+name,sha(audit.file(p))==want)
    review=audit.load(a.source_review)
    audit.check('independent51member source review exact',sha(a.source_review)==REVIEW_SHA and freeze['peer_review_sha256']==REVIEW_SHA and review['status']=='APPROVED_FOR_SOURCE_FREEZE_ONLY' and dict(review['files'],SHA256SUMS=sha(source/'SHA256SUMS'))==freeze['files'])
    archive=ROOT/freeze['archive'];audit.check('exact52member archive hash',sha(audit.file(archive))==ARCHIVE_SHA)
    with tarfile.open(archive,'r:gz') as tar:
        seen=set()
        for m in tar.getmembers():
            rel=PurePosixPath(m.name);require(m.isfile() and not m.issym() and rel.parts[0]==RELEASE and '..' not in rel.parts,'Safe archive member')
            name='/'.join(rel.parts[1:]);require(name in freeze['files'] and name not in seen,'Unique expected archive member');seen.add(name)
            audit.check('archive source bytes '+name,hashlib.sha256(tar.extractfile(m).read()).hexdigest()==freeze['files'][name])
        audit.check('complete archive closure',seen==set(freeze['files']))
    cfg=audit.load(source/'manifest.toml');closure=audit.load(source/'closure.json');seeds=audit.load(source/'seeds.json')
    expected=dict(schema='P21_SAMPLED_PRODUCTION_V1',phase='PRODUCTION_SINGLE_LOOK',catalogue_id='P21',L=4,Ltau=10,t=1.,U=2.,mu_centered=0.,beta=1.,dtau=.1,n_stab=2,boundary='pbc',chains_per_backend=8,retained_bins=64,cycles_per_bin=1024,burn_bins=4,burn_cycles=4096,retained_cycles=65536,alf_NSweep=1024,alf_NBin=68,alf_discard_initial_bins=4,alf_Ltau_flag=1,alf_LOBS_ST=1,alf_LOBS_EN=10,user_update_calls_per_cycle=2,orders=['GT0','G0T'],component_parts=['real','imag'],real_components_total=1408,complex_components_per_order=352)
    audit.equal('fixed complete model/timegrid/budget',{k:cfg.get(k) for k in expected},expected)
    st=cfg['statistics'];audit.check('fixed one-look Bonferroni1408 Student7',st['looks']==1 and st['family_alpha']==.01 and st['real_intervals']==1408 and st['student_df']==7 and st['student_critical']==CRITICAL and st['complex_equivalence_margin']==.02 and st['complex_precision_half_radius']==.01)
    audit.check('independent Student7 elementary tail formula',abs(mixing.student7_cdf(CRITICAL)-(1-.01/(2*1408)))<16*math.ulp(1.))
    limits=dict(positive_sign_abs=1e-10,user_cumulative_deltaG_max=1e-8,user_endpoint_max_error=5e-10,native_green_max=1e-8,native_tau_max=1e-8,native_phase_max=1e-10,required_count_per_tau_per_bin=1024,structural_imaginary_raw_abs_max=0.,rhat_max=1.05,ess_min_per_chain=16.,ess_min_total=128.,blocking_error_relative_change_max=.5,real_constant_chain_allowed=False,requires_complete_1408_component_grid=True)
    audit.equal('all fixed numerical and statistical gates',{k:cfg['eligibility'].get(k) for k in limits},limits)
    audit.check('same proved observer and unchanged diagnostic math',sha(source/'UserMatrixObserver.jl')=='c6e3462ccaa54d2f19a0c58202f222da4f099aa431712247545608933e2fa807' and sha(source/'mixing_diagnostics.py')=='bc4496b44ea537cfb8310b7fae181a3dfd5b2a1093cf165d8c9c21bda04d8451')
    for name,want in closure['assets'].items():audit.check('preserved proof/planning/license '+name,sha(source/name)==want)
    for name,want in closure['user_base_pins'].items():audit.check('actual six native source pins '+name,sha(audit.file(a.native_source_dir/name))==want)
    rows=seeds['chains'];seedmap={(r['backend'],r['chain']):r['seed'] for r in rows}
    audit.check('all16 fixed fresh range-valid seeds',len(seedmap)==len(rows)==16 and set(seedmap)=={(b,i) for b in ('alf','user') for i in range(1,9)} and len(set(seedmap.values()))==16 and not(set(seedmap.values())&set(seeds['excluded_prior_seeds'])) and all(0<v<(2**31 if b=='alf' else 2**63) for (b,i),v in seedmap.items()))
    for p in source.iterdir():
        if p.suffix in ('.py','.jl','.sbatch','.toml','.in') or p.name=='SHA256SUMS':audit.check('raw LF executable/config '+p.name,b'\r' not in p.read_bytes())
    return cfg,closure,seedmap,freeze


def raw_file_list(result):
    paths=[result/name for name in ('execution.json','preflight.json','user_sampling.toml','user_sampling.resources.txt','user_sampling.stdout.log')]
    for i in range(1,9):
        paths += [result/f'alf_chain_{i}'/name for name in ('parameters','seeds.initial','seeds','info','P21_GT0_scal','P21_G0T_scal','P21_Counts_scal')]
        paths += [result/f'alf_chain_{i}{s}' for s in ('.resources.txt','.stdout.log')]
        paths.append(result/f'user_sampling.toml.chain{i}.toml')
        paths += [result/'user_sampling.toml.bins'/f'chain{i}_rawbin{j:03}.toml' for j in range(1,69)]
    require(len(paths)==len(set(paths))==629,'Exact629 explicit raw roles')
    return paths


def execution(audit,mixing,a,cfg,closure,freeze):
    directory=a.job_dir.resolve();result=directory/'results';record=audit.load(directory/'transfer.json')
    audit.check('original read-only job transfer',str(record['job_id'])==JOB and record.get('remote_writes') is False and record.get('binaries_copied') is False)
    seen=set()
    for row in record['files']:
        p=(directory/row['local']).resolve();require(p.is_relative_to(directory) and p not in seen,'Safe unique transfer path');seen.add(p)
        audit.check('original transfer '+row['local'],row['status']=='RETRIEVED' and audit.file(p).stat().st_size==row['bytes'] and sha(p)==row['sha256'])
    for name,want in freeze['files'].items():
        audit.check('retrieved executed frozen source '+name,sha(audit.file(directory/'source'/name))==want)
    audit.check('retrieved administrative freeze pin',sha(audit.file(directory/'source/source_freeze.json'))==FREEZE_SHA)
    scheduler=audit.file(directory/'scheduler.txt').read_text(encoding='utf-8-sig');audit.check('scheduler text exact transfer binding',scheduler==record['sacct'])
    rows=list(csv.DictReader(io.StringIO(scheduler),delimiter='|'));main=[r for r in rows if r.get('JobIDRaw')==JOB];batch=[r for r in rows if r.get('JobIDRaw')==JOB+'.batch']
    audit.check('exact scheduler identities completed0:0',len(main)==len(batch)==1 and main[0]['State']==batch[0]['State']=='COMPLETED' and main[0]['ExitCode']==batch[0]['ExitCode']=='0:0')
    require(main and batch,'Scheduler main/batch required');m=main[0]
    audit.check('fixed1CPU3500M120minute allocation',m['ReqCPUS']==m['AllocCPUS']=='1' and m['ReqMem']=='3500M' and mixing.duration_seconds(m['Timelimit'])==7200 and 0<mixing.duration_seconds(m['Elapsed'])<=7201)
    ex=audit.load(result/'execution.json');pre=audit.load(result/'preflight.json')
    audit.check('actual complete one-look identity',ex['schema']=='P21_PRODUCTION_EXECUTION_V1' and ex['status']=='SAMPLED_PRODUCTION_COMPLETE' and ex['job_id']==JOB and ex['one_prespecified_look'] is True and ex['scientific_completions']==0)
    expected=[f'alf_chain_{i}' for i in range(1,9)]+['user_sampling'];procs=ex['processes']
    audit.check('exact9 serial ordinary sampled process roles',[r['name'] for r in procs]==expected and all(r['exit_code']==0 and r['timed_out'] is False and r['proof_mode_present'] is False for r in procs))
    audit.check('absolute science cutoff and full120min actual-start admission',ex['scientific_cutoff_epoch']==CUTOFF and ex['hard_stop_epoch']==ex['actual_batch_start_epoch']+7200<=CUTOFF and ex['soft_stop_epoch']==ex['hard_stop_epoch']-5)
    audit.check('serial native/user endpoints before hard cutoff',all(math.isfinite(r['started_epoch']) and r['started_epoch']<r['ended_epoch']<=ex['hard_stop_epoch'] and r['soft_stop_epoch']==ex['soft_stop_epoch'] and r['timeout_seconds']>0 for r in procs) and all(x['ended_epoch']<=y['started_epoch'] for x,y in zip(procs,procs[1:])))
    resources={}
    for r in procs:
        cap=900 if r['name'].startswith('alf_') else 1800
        resources[r['name']]=mixing.resource_receipt(audit,result/(r['name']+'.resources.txt'),cap)
        audit.check('timed native process success '+r['name'],resources[r['name']]['exit_code']==0 and 0<r['elapsed_seconds']<=cap+5)
    audit.check('sum scientific process walls bounded by allocation',sum(r['wall_seconds'] for r in resources.values())<=mixing.duration_seconds(m['Elapsed'])+3)
    audit.check('preflight identity binary/provider/sourcefreeze',pre['status']=='SOURCE_RUNTIME_AND_PRIOR_PROOFS_VERIFIED' and pre['job_id']==JOB and pre['binary_sha256']==closure['binary_sha256']==freeze['native_binary_sha256'] and pre['provider_sha256']==closure['provider_sha256'] and pre['source_freeze_sha256']==FREEZE_SHA and pre['proof_mode_unset_in_sampled_children'] is True)
    pins={str(PurePosixPath(pre['proof_root'])/k):v for k,v in closure['remote_proof_pins'].items()};pins.update({str(PurePosixPath(pre['user_source_root'])/k):v for k,v in closure['user_base_pins'].items()})
    audit.equal('fresh native/proof/source byte-pin receipt',pre['actual_source_receipt_hashes'],pins)
    expected_files={k:v for k,v in freeze['files'].items() if k!='SHA256SUMS'}
    audit.equal('executed complete52member source closure',pre['release_files'],expected_files)
    audit.check('executed checksum/provider library closure',pre['release_sha256sums_sha256']==sha(a.release_dir/'SHA256SUMS') and pre['provider']==closure['provider'] and 'not found' not in pre['ldd'])
    provider_lines=[r for r in pre['ldd'].splitlines() if any(n in r for n in ('libopenblas','libblas','liblapack'))]
    audit.check('single pinned provider in actual linker receipt',bool(provider_lines) and all(str(PurePosixPath(closure['provider']).parent)+'/' in r for r in provider_lines))
    return ex,pre,resources,dict(zip(expected,procs))


def check_bin(row):
    require(row['counts']==[1024]*11 and all(type(x) is int for x in row['counts']),'Exact1024 actual callbacks at all11tau')
    require(len(row['sign'])==11 and all(math.isfinite(x) and abs(x/1024-1)<=1e-10 for x in row['sign']),'Full positive paired sign sums')
    for order in ('GT0','G0T'):
        for part in ('real','imag'):
            values=row[order+'.'+part];require(len(values)==352 and all(math.isfinite(x) for x in values),'Complete finite352vector')
            if part=='imag':require(all(x==0 for x in values),'Exact all-burn/retained raw structural imaginary zero')


def reconstruct_alf(audit,pilot,result,source,i,seed,process):
    folder=result/f'alf_chain_{i}';require((folder/'parameters').read_bytes()==(source/'parameters.in').read_bytes(),'Exact native parameters')
    require((folder/'seeds.initial').read_bytes()==(str(seed)+'\n').encode(),'Exact native initial seed')
    require(process['seed']==seed and process['initial_seed_sha256']==sha(folder/'seeds.initial') and process['parameters_sha256']==sha(folder/'parameters'),'Native process seed/parameter receipts')
    gt=native_vectors(folder/'P21_GT0_scal',352);gz=native_vectors(folder/'P21_G0T_scal',352);ct=native_vectors(folder/'P21_Counts_scal',11)
    bins=[];endpoint=0.
    for j,((g,sg),(z,sz),(counts,cs)) in enumerate(zip(gt,gz,ct),1):
        require(sg==sz and abs(sg-1)<=1e-10 and cs==1 and counts==[1024+0j]*11,'Native paired phase/count normalization')
        row={'counts':[1024]*11,'sign':[1024*sg]*11}
        for order,values in [('GT0',g),('G0T',z)]:
            row[order+'.real']=[1024*x.real for x in values];row[order+'.imag']=[1024*x.imag for x in values]
        check_bin(row)
        for spin in range(2):
            for r in range(4):
                for c in range(4):
                    a,b=pilot.index(0,spin,r,c),pilot.index(10,spin,r,c);identity=float(r==c)
                    endpoint=max(endpoint,abs(z[a]-g[a]+identity*sg),abs(g[b]+g[a]-identity*sg),abs(z[b]+g[a]))
        if j>4:bins.append(row)
    audit.check('native full68bin endpoints '+str(i),len(bins)==64 and endpoint<=5e-10)
    text=(folder/'info').read_text(encoding='utf-8')
    audit.check('native actual fresh seed '+str(i),re.findall(r'No initial configuration, Seed_in\s+(\d+)',text)==[str(seed)] and 'FIN TIME:' in text and 'f4c500cc' in text)
    for label,value in [('Sweeps',1024),('Bins',68)]:audit.check('native actual '+label+' '+str(i),re.findall(label+r'\s*:\s*(\d+)',text)==[str(value)])
    for name,pattern,limit in [('Green',r'Precision Green\s+Mean, Max\s*:\s*(\S+)\s+(\S+)',1e-8),('tau',r'Precision tau\s+Mean, Max\s*:\s*(\S+)\s+(\S+)',1e-8),('phase',r'Precision Phase, Max\s*:\s*(\S+)',1e-10)]:
        found=re.findall(pattern,text);require(len(found)==1,'One final native numerical diagnostic '+name)
        values=found[0] if isinstance(found[0],tuple) else (found[0],)
        values=[float(v.replace('D','E').replace('d','e')) for v in values]
        audit.check('native numerical '+str(i)+'/'+name,all(math.isfinite(v) and 0<=v<=limit for v in values))
    return bins,endpoint


def reconstruct_user(audit,pilot,result,source,seedmap,execution):
    path=result/'user_sampling.toml';target=audit.load(path);pilot.finite(target)
    audit.check('actual complete user production identity',target['schema']=='P21_SAMPLED_PRODUCTION_V1' and target['phase']=='PRODUCTION_SINGLE_LOOK' and target['status']=='SAMPLED_PRODUCTION_COMPLETE' and target['job_id']==JOB and target['scientific_completions']==0 and target['scientific_soft_stop_epoch']==execution['soft_stop_epoch'])
    for field,name in [('manifest_sha256','manifest.toml'),('driver_sha256','user_sampling.jl'),('observer_sha256','UserMatrixObserver.jl')]:audit.check('actual user source '+name,target[field]==sha(source/name))
    audit.check('actual user native package pins',target['original_kernel_sha256']=='b72bbbf158381b0f776c1fa259d532b069094c47e361d3490c19ea0f43ada4be' and target['observer_kernel_sha256']=='911a61642b23a4465639c3b3e4c63ff9d7cf05d51da288f3224f1775c0aba891' and target['julia_version']=='1.12.1' and target['jdqmcframework_version']=='1.2.6' and target['julia_threads']==target['blas_threads']==1)
    chains=target['chains'];require(len(chains)==8 and [r['chain'] for r in chains]==list(range(1,9)),'All8 fixed ordered user chains')
    output=[]
    for chain in chains:
        i=chain['chain'];require(chain['seed']==seedmap['user',i] and chain['status']=='SAMPLED_CHAIN_COMPLETE' and chain['burn_cycles']==4096 and chain['retained_cycles']==65536 and chain['update_calls']==2*69632,'Full new user chain/seed/cycle budget')
        eq(audit,pilot,'immutable whole-chain receipt '+str(i),audit.load(result/f'user_sampling.toml.chain{i}.toml'),chain)
        allrows=chain['burn_bins']+chain['bins'];require(len(chain['burn_bins'])==4 and len(chain['bins'])==64 and [r['raw_bin'] for r in allrows]==list(range(1,69)),'All68 raw user bins including4burn')
        deltas=[];rows=[];endpoint=0.
        for j,r in enumerate(allrows,1):
            require(r['chain']==i and r['seed']==chain['seed'] and r['role']==('BURN_EXCLUDED' if j<=4 else 'RETAINED') and r['bin']==(j if j<=4 else j-4),'User raw bin role/identity')
            require(r['primary_stage']=='initial(tau0) and pre_stabilization(tau1..M)' and r['measurement_origin']=='actual user maintained buffers; no dense recomputation' and r['ordering']=='tau-major,spin(up,down),row_i,column_j; all indices zero-based externally','Actual maintained moving-cut measurement')
            row={'counts':r['counts_by_tau'],'sign':r['sign_denominator_sums_by_tau']}
            for order in ('GT0','G0T'):
                for part in ('real','imag'):row[order+'.'+part]=r[order+'_numerator_'+part+'_sums']
            check_bin(row);deltas.append(r['cumulative_deltaG']);endpoint=max(endpoint,r['endpoint_contemporaneous_max_error'])
            eq(audit,pilot,f'immutable raw user bin chain{i}/{j}',audit.load(result/'user_sampling.toml.bins'/f'chain{i}_rawbin{j:03}.toml'),r)
            if j>4:rows.append(row)
        audit.check('user cumulative burn+retained error and endpoints '+str(i),all(math.isfinite(x) and 0<=x<=1e-8 for x in deltas) and all(x<=y for x,y in zip(deltas,deltas[1:])) and deltas[3]==chain['cumulative_deltaG_after_burn'] and deltas[-1]==chain['cumulative_deltaG'] and math.isfinite(endpoint) and 0<=endpoint<=5e-10 and 0<=chain['acceptance']<=1 and 0<chain['elapsed_seconds']<=1800)
        output.append(rows)
    return output


def interval(differences,structural):
    require(len(differences)==8 and all(math.isfinite(v) for v in differences),'Eight finite pair differences')
    mean=math.fsum(differences)/8;se=math.sqrt(math.fsum((v-mean)**2 for v in differences)/56);half=CRITICAL*se
    eligible=structural or se>0
    if structural:require(all(v==0 for v in differences) and mean==se==half==0,'Source-bound structural imaginary zero')
    return dict(pair_differences=differences,mean=mean,standard_error=se,half_width=half,lower=mean-half,upper=mean+half,structural_zero=structural,eligible=eligible,issues=[] if eligible else ['ZERO_BETWEEN_PAIR_VARIANCE_NONSTRUCTURAL_COMPONENT'])


def geometry(real,imag):
    half=math.hypot(real['half_width'],imag['half_width']);far=math.hypot(abs(real['mean'])+real['half_width'],abs(imag['mean'])+imag['half_width']);near=math.hypot(max(0.,abs(real['mean'])-real['half_width']),max(0.,abs(imag['mean'])-imag['half_width']))
    return dict(half_radius=half,farthest_rectangle_corner_norm=far,nearest_rectangle_point_norm=near,precision_met=half<=.01,rectangle_inside_equivalence_disk=far<=.02,rectangle_wholly_outside_equivalence_disk=near>.02)


def decision(eligible,grid):
    if not eligible or not all(r['precision_met'] for r in grid):return 'UNRESOLVED'
    if all(r['rectangle_inside_equivalence_disk'] for r in grid):return 'PASS'
    if any(r['rectangle_wholly_outside_equivalence_disk'] for r in grid):return 'FAIL'
    return 'UNRESOLVED'


def statistics_from_raw(audit,pilot,mixing,alf,user,saved):
    coordinates=[(o,p,t,s,i,j) for o in ('GT0','G0T') for p in ('real','imag') for t in range(11) for s in range(2) for i in range(4) for j in range(4)]
    require(len(saved['components'])==1408 and len(saved['complex_grid'])==704,'Full saved1408/704grid')
    components=[];indexed={};all_eligible=True
    for position,coordinate in enumerate(coordinates):
        order,part,tau,spin,i,j=coordinate;structural=part=='imag';ix=pilot.index(tau,spin,i,j)
        record=dict(order=order,part=part,tau_index=tau,tau=tau*.1,spin='up' if spin==0 else 'down',row=i,column=j)
        reported=saved['components'][position];eq(audit,pilot,'exact component coordinate '+str(position),{k:reported[k] for k in record},record)
        estimates={};eligible=True
        for backend,chains in [('alf',alf),('user',user)]:
            series=[];means=[];blocked=[];changes=[]
            for chain_id,bins in enumerate(chains,1):
                nums=[r[order+'.'+part][ix] for r in bins];dens=[r['sign'][tau] for r in bins]
                fine=pilot.paired_ratio_statistics(nums,dens);coarse=pilot.paired_ratio_statistics([x+y for x,y in zip(nums[::2],nums[1::2])],[x+y for x,y in zip(dens[::2],dens[1::2])])
                values=[x/y for x,y in zip(nums,dens)];series.append(values);means.append(fine['plug_in'])
                se,co=fine['standard_error'],coarse['standard_error'];change=abs(co/se-1) if se>0 else None;good=structural or (se>0 and co>0 and change<=.5);blocked.append(good);changes.append(change)
                expected=dict(full_ratio=fine['plug_in'],jackknife_bias_corrected=fine['bias_corrected'],standard_error=se,coarsened32_standard_error=co,observed_constant_bins=min(values)==max(values))
                reported_chain=reported[backend]['chains'][chain_id-1]
                eq(audit,pilot,f'independent paired ratio/blocking {position}/{backend}/{chain_id}',{k:reported_chain[k] for k in expected},expected)
                audit.check(f'exact blocking decision {position}/{backend}/{chain_id}',reported_chain['blocking']['eligible']==good)
            if structural:
                diagnostic=dict(eligible=True,status='SOURCE_BOUND_STRUCTURAL_ZERO',rhat=None,ess_total=None,ess_per_chain=[])
                audit.check(f'structural source-bound diagnostic {position}/{backend}',reported[backend]['mixing']['eligible'] is True and reported[backend]['mixing']['status']=='SOURCE_BOUND_STRUCTURAL_ZERO')
            else:
                diagnostic=mixing.own_mixing(series)
                eq(audit,pilot,f'independent rank/folded/Geyer diagnostics {position}/{backend}',{k:reported[backend]['mixing'][k] for k in diagnostic},diagnostic)
            good=bool(diagnostic['eligible'] and all(blocked));eligible &= good
            audit.check(f'backend eligibility decision {position}/{backend}',reported[backend]['eligible']==good)
            record[backend]=dict(chain_full_ratios=means,mean=math.fsum(means)/8,mixing=diagnostic,blocking_eligible=blocked,blocking_relative_changes=changes,eligible=good)
            estimates[backend]=means
        own=interval([r-l for l,r in zip(estimates['alf'],estimates['user'])],structural)
        eq(audit,pilot,'independent Student7 pair interval '+str(position),reported['difference_interval'],own)
        record['difference_interval']=own;record['eligible']=bool(eligible and own['eligible']);all_eligible &= record['eligible']
        audit.check('exact real component eligibility '+str(position),reported['eligible']==record['eligible']);components.append(record);indexed[coordinate]=record
    grid=[];position=0
    for order in ('GT0','G0T'):
        for tau in range(11):
            for spin in range(2):
                for i in range(4):
                    for j in range(4):
                        real,imag=[indexed[order,p,tau,spin,i,j]['difference_interval'] for p in ('real','imag')]
                        row=dict(order=order,tau_index=tau,tau=tau*.1,spin='up' if spin==0 else 'down',row=i,column=j,**geometry(real,imag))
                        eq(audit,pilot,'independent complex rectangle disk gate '+str(position),saved['complex_grid'][position],row);grid.append(row);position+=1
    status=decision(all_eligible,grid);candidate={'PASS':'PASS_CANDIDATE_REQUIRES_INDEPENDENT_AUDIT','FAIL':'FAIL_CANDIDATE_REQUIRES_INDEPENDENT_AUDIT','UNRESOLVED':'UNRESOLVED'}[status]
    audit.check('independent complete scientific decision',saved['status']==candidate and saved['mixing_blocking_and_pair_variance_eligible']==all_eligible and saved['all_complex_precision_met']==all(r['precision_met'] for r in grid) and saved['all_complex_rectangles_inside_margin']==all(r['rectangle_inside_equivalence_disk'] for r in grid) and saved['any_complex_rectangle_wholly_outside_margin']==any(r['rectangle_wholly_outside_equivalence_disk'] for r in grid))
    return dict(components=components,complex_grid=grid,globally_eligible=bool(all_eligible),all_precision_met=all(r['precision_met'] for r in grid),scientific_status=status,candidate_status=candidate,real_components=1408,complex_entries=704,independent_units=8,df=7,student_critical=CRITICAL,family_alpha=.01,simultaneous_real_intervals=1408)


def replay(audit,pilot,a,saved):
    require(all(r['passed'] for r in audit.checks),'Integrity rejection: do not execute any mutable production code')
    with tempfile.TemporaryDirectory(prefix='p21_production_analysis_replay_',dir=a.output.parent) as temp:
        output=Path(temp)/'comparison.json';run=subprocess.run([sys.executable,'-B',str(a.release_dir/'compare_production.py'),'--result-dir',str(a.job_dir/'results'),'--output',str(output)],capture_output=True,text=True,timeout=180)
        audit.check('frozen pure-analysis replay exit0',run.returncode==0,dict(exit_code=run.returncode,stderr=run.stderr[-2000:]))
        other=audit.load(output);agree=eq(audit,pilot,'full frozen production replay within64epsilon; exact labels',saved,other)
        # The temporary replay is diagnostic output, not a portable input receipt.
        audit.inputs.pop(str(output.resolve()),None)
        return dict(exit_code=run.returncode,saved_sha256=sha(a.job_dir/'results/comparison.json'),replay_sha256=sha(output),gates_and_labels_exact=agree,all_output_agrees_within_64_epsilon=agree,byte_identical=output.read_bytes()==(a.job_dir/'results/comparison.json').read_bytes())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--job-dir',type=Path,required=True)
    p.add_argument('--release-dir',type=Path,default=ROOT/RELEASE)
    p.add_argument('--freeze',type=Path,default=ROOT/'p21-sampled-production-v1-freeze.json')
    p.add_argument('--source-review',type=Path,default=ROOT/'next_stage_review/p21_production_v1_peer_review.json')
    p.add_argument('--native-source-dir',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.job_dir=a.job_dir.resolve();a.release_dir=a.release_dir.resolve();a.output=a.output.resolve();a.native_source_dir=(a.native_source_dir or a.job_dir/'native_source').resolve()
    require(not a.output.exists(),'Choose new audit output; retain failures');a.output.parent.mkdir(parents=True,exist_ok=True)
    pilot=source_module(PILOT_TOOL,PILOT_SHA,'p21_previous_independent_math');mixing=source_module(MIX_TOOL,MIX_SHA,'p21_independent_rank_math');audit=compact_audit(mixing)
    report=dict(schema=SCHEMA,job_id=JOB,catalogue_id='P21',created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),auditor_sha256=sha(__file__),approved_freeze_sha256=FREEZE_SHA,independent_math_sources={str(PILOT_TOOL):PILOT_SHA,str(MIX_TOOL):MIX_SHA},scientific_completions=0,scientific_sampling=False,remote_calls=False)
    try:
        audit.file(PILOT_TOOL);audit.file(MIX_TOOL)
        cfg,closure,seeds,freeze=protocol(audit,pilot,mixing,a)
        require(all(r['passed'] for r in audit.checks),'Source/protocol integrity rejected; production code replay is forbidden')
        ex,pre,resources,procs=execution(audit,mixing,a,cfg,closure,freeze);result=a.job_dir/'results'
        files=raw_file_list(result);hashes={q.relative_to(result).as_posix():sha(audit.file(q)) for q in files};saved=audit.load(result/'comparison.json');pilot.finite(saved)
        audit.check('full frozen comparison source/runtime identity',saved['schema']=='P21_SAMPLED_PRODUCTION_V1' and saved['phase']=='PRODUCTION_SINGLE_LOOK' and saved['job_id']==JOB and saved['scientific_completions']==0 and saved['manifest_sha256']==sha(a.release_dir/'manifest.toml') and saved['seed_sha256']==sha(a.release_dir/'seeds.json') and saved['comparison_tool_sha256']==sha(a.release_dir/'compare_production.py') and saved['receipt_io_sha256']==sha(a.release_dir/'receipt_io.py') and saved['mixing_tool_sha256']==sha(a.release_dir/'mixing_diagnostics.py'))
        audit.check('complete one-look whole-grid claimed scope',saved['catalogue_id']=='P21' and saved['independent_audit_required_before_catalogue_change'] is True and saved['one_prespecified_look'] is True and saved['real_components_retained']==1408 and saved['complex_entries_per_order']==352 and saved['chains_per_backend']==8 and saved['retained_bins_per_chain']==64 and saved['cycles_per_bin']==1024 and saved['burn_cycles']==4096 and saved['retained_cycles']==65536 and saved['student_df']==7 and saved['family_alpha']==.01 and saved['bonferroni_real_interval_count']==1408 and saved['student_critical']==CRITICAL and saved['complex_equivalence_margin']==.02 and saved['complex_precision_half_radius_max']==.01)
        audit.equal('exact629 raw whitelist/comparator hash closure',saved['input_sha256'],hashes)
        alf=[];errors=[]
        for i in range(1,9):
            audit.check('native executed proved binary '+str(i),procs[f'alf_chain_{i}']['argv'][-1]==pre['binary'])
            rows,endpoint=reconstruct_alf(audit,pilot,result,a.release_dir,i,seeds['alf',i],procs[f'alf_chain_{i}']);alf.append(rows);errors.append(endpoint)
        user=reconstruct_user(audit,pilot,result,a.release_dir,seeds,ex)
        eq(audit,pilot,'independent native endpoint maximum',saved['alf_endpoint_max_error'],max(errors))
        report['independent_reconstruction']=statistics_from_raw(audit,pilot,mixing,alf,user,saved)
        report['candidate_status']=saved['status'];report['resources']=resources
        if all(r['passed'] for r in audit.checks):
            report['frozen_comparator_replay']=replay(audit,pilot,a,saved)
        else:
            report['frozen_comparator_replay']=dict(skipped=True,reason='Independent integrity rejection; no production code executed')
        report['raw_input_count']=len(files)
    except Exception as error:audit.check('complete independent raw/statistical audit',False,type(error).__name__+': '+str(error))
    for name,want in list(audit.inputs.items()):audit.check('bound input unchanged at audit end '+name,Path(name).is_file() and sha(name)==want)
    audit.finalize()
    failures=[r for r in audit.checks if not r['passed']]
    report.update(audit_status='PASS' if not failures else 'UNRESOLVED',scientific_status=report.get('independent_reconstruction',{}).get('scientific_status','UNRESOLVED') if not failures else 'UNRESOLVED',check_count=len(audit.checks),atomic_check_count=audit.atomic_check_count,aggregate_grid_checks=audit.groups,failure_count=len(failures),failures=failures,checks=audit.checks,input_sha256=audit.inputs,pilot_math_scalar_checks=pilot.CHECKS,maximum_64epsilon_numeric_difference=pilot.MAX_NUMERIC_DIFFERENCE,limitations=['Native MonteCarlo/observer kernels are source-bound,not reimplemented.','Independent pilot parser grammar is explicitly budget-adapted to68rawbins/1024cycles;old pilot source and math remain unchanged.','Eight fixed independent chain-pair differences/7df;64bins do not create64Student replicates.','Nominal high-tail Student/CLT coverage is conditional on stationary approximately normal independent chain estimates;finiteRhat/ESS/blocking do not prove equilibrium,rare-tail coverage or precision.','Imaginary zero coverage rests on the fixed real-operator theorem/proof plus exact all-burn/all-retained rawzero tests;sampleconstancy alone is insufficient.','Native binary/provider remain onISAAC;their recorded live preflight hashes/linkage are byte-bound,not independently loaded or executed here.','This finite-dtau two-spin-block comparison includes analytic conserved cross-spinzero blocks;no general spin-mixing claim.','No ledger or checkpoint changed;separate root review/update is mandatory.'])
    with a.output.open('x',encoding='utf-8') as f:json.dump(mixing.plain(report),f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')
    print(json.dumps(dict(audit_status=report['audit_status'],scientific_status=report['scientific_status'],checks=len(audit.checks),failures=len(failures),output=str(a.output)),indent=2))
    if failures:raise SystemExit(2)


if __name__=='__main__':main()

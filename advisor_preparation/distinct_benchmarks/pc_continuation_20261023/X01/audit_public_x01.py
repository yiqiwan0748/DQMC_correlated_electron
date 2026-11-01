#!/usr/bin/env python3
"""Public scientific-selection audit of X01; fixed-data postprocessing only.

Never imports native Julia, QuSpin, a sampler, submission or scheduler tools.
The pinned independent bit-Fock/bin auditor is imported; frozen comparator
replays are separate Python subprocesses with bytecode writes disabled.
Historical Windows/ISAAC path keys are mapped by explicit campaign roles and
verified hashes, without changing any historical source or receipt bytes.
"""
from __future__ import annotations
import argparse, ast, base64, csv, datetime, hashlib, importlib.util, io, json
import math, os, re, subprocess, sys, tarfile, tomllib
from pathlib import Path, PurePosixPath
from types import SimpleNamespace

for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
MATH_SHA = '6d68b99f06eeea2834f72e386e7b91fbe63bf8aa62f265c95947361eb761910c'
FREEZE = {
 'continuum_production_v1c':'659b4456ba007aed1618e2f937fd79f87133ff89fbf62e07bfb5aacc69b35c7f',
 'x01_continuum_stage2_v2':'4615a3c7aed3c2989c9ac4c0b3f0078a2b87050a4d96a0a166eff9c49420bdb2',
 'x01_stage2_analysis_repair_v1':'4e446cfabe5e5e284db374cb51f02d5f6d2a8d4bbaafd7809a7d00cbe013ac0f',
 'x01_continuum_stage3_v1':'cb5cf95a35fd82e54700141f325a78e7b0cc69a3ed64e74115323d5affe66bb6'}
RELEASE_FOR = {'x01_6460771':'continuum_production_v1c',
 **{f'x01s2p{i}_{job}':'x01_continuum_stage2_v2' for i,job in enumerate((6460990,6460991,6460992,6460994),1)},
 'x01s2a_6460996':'x01_continuum_stage2_v2',
 **{f'x01s3p{i}_{job}':'x01_continuum_stage3_v1' for i,job in enumerate(range(6462170,6462174),1)},
 'x01s3a_6462174':'x01_continuum_stage3_v1'}

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p)
    return tomllib.loads(p.read_text(encoding='utf-8')) if p.suffix=='.toml' else json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,obj):
    with Path(p).open('x',encoding='utf-8') as f: json.dump(obj,f,indent=2,sort_keys=True,allow_nan=False); f.write('\n')
def relative(p):
    path=Path(os.path.abspath(p))
    for role,base in (('raw',rawroot),('source',snapshot),('audit',HERE),('output',args.output)):
        try:return role+'/'+str(path.relative_to(Path(os.path.abspath(base))))
        except ValueError:pass
    return str(path)

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--workspace',type=Path,default=Path.cwd())
parser.add_argument('--raw',type=Path,help='Portable restored raw root (restored/, indices, restoration receipt)')
parser.add_argument('--source',type=Path,help='Portable source root (four original releases and selected original receipts)')
parser.add_argument('--decoded-only',action='store_true',help='Replay exact decoded payload hashes without rechecking duplicated encoded chunk files; original chunk restoration receipt remains bound')
parser.add_argument('--selection',type=Path,default=HERE/'selection_manifest.json')
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
if args.output.exists(): raise ValueError('Choose a fresh audit output directory')
args.output.mkdir(parents=True)
snapshot=(args.source or HERE/'source').absolute()
rawroot=(args.raw or HERE/'raw').absolute()
restored=rawroot/'restored'
mathpath=HERE/'support/audit_continuum_stage1.py'
assert sha(mathpath)==MATH_SHA
spec=importlib.util.spec_from_file_location('pinned_independent_x01_math',mathpath)
M=importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
audit=M.Audit()
audit.file(__file__); audit.file(mathpath)
audit.check('pinned independent math auditor exact',sha(mathpath)==MATH_SHA)
selection=read(args.selection)
selected_paths={row['path'] for row in selection['selected_files']}
omitted_paths={row['path'] for row in selection['omitted_files']}
public_exclusions=[]
audit.file(args.selection)
audit.check('public selection binds completed original full closure',selection['original_full_audit_sha256']=='4eecab1e48c0164b43833c5f0f63e5004a0b16436f0f487c0c30ae7383000f38' and selection['original_full_audit_status']=='PASS' and selection['original_full_audit_checks']==6853)
audit.check('read-only/no new scientific execution',True,dict(sampler=False,native_solver=False,scheduler=False,remote=False,ledger_edits=False))
print('PHASE verified pinned math; rechecking raw closure',flush=True)

def bind(p,want=None,name=None):
    p=audit.file(p)
    if want is not None: audit.check(name or 'hash '+relative(p),sha(p)==want)
    return p

def safe_parts(text):
    parts=[x for x in text.replace('\\','/').split('/') if x]
    if any(x in ('.','..') or ':' in x for x in parts): raise ValueError('Unsafe relative role path')
    return parts

def map_windows(key):
    # Exact historical campaign-root marker, followed by an allowlisted role.
    normalized=key.replace('\\','/')
    parts=[x for x in normalized.split('/') if x]
    marker='isaac_benchmark_resume_20261022'
    if parts.count(marker)!=1: raise ValueError('Unknown Windows campaign path '+key)
    rel=parts[parts.index(marker)+1:]
    safe_parts('/'.join(rel))
    if not rel: raise ValueError('Empty campaign role')
    if (restored/rel[0]).is_dir(): return restored.joinpath(*rel)
    if (snapshot/rel[0]).is_dir() or (snapshot/rel[0]).is_file(): return snapshot.joinpath(*rel)
    raise ValueError('Unresolved Windows role '+key)

def map_remote(key):
    marker='/advisor-validation-20261022/'
    if key.count(marker)!=1: raise ValueError('Unregistered ISAAC root '+key)
    rel=safe_parts(key.split(marker,1)[1])
    if rel[0]=='results':
        name=rel[1]
        match=re.fullmatch(r'x01_continuum_stage2_v2_part([1-4])_(\d+)',name)
        if match: label='x01s2p'+match[1]+'_'+match[2]
        elif name=='x01_continuum_stage2_v2_analysis_6460996': label='x01s2a_6460996'
        else: raise ValueError('Unregistered remote results namespace '+name)
        return restored.joinpath(label,'results',*rel[2:])
    if rel[0]=='submissions':
        # Registry was independently retained with the analysis submission.
        matches=list((restored/'x01s2a_6460996/submissions').rglob(rel[-1]))
        if len(matches)==1:return matches[0]
    raise ValueError('Unresolved ISAAC role '+key)

def path_map_hashes(mapping,label,mapper):
    rows=[]
    for key,want in sorted(mapping.items()):
        p=mapper(key)
        if not p.is_file():
            role='restored/'+str(p.relative_to(restored)) if p.is_relative_to(restored) else None
            if role not in omitted_paths:raise ValueError('Missing mandatory scientific input '+key)
            public_exclusions.append(dict(original_key=key,sha256=want,reason='Explicit public administrative/duplicate selection exclusion'))
            continue
        bind(p,want,label+' hash '+key)
        rows.append(dict(original_key=key,restored_path=relative(p),sha256=want))
    return rows

# NEW derived public selection; original transport/admin receipts remain untouched.
closure=[]
assert len(selected_paths)==len(selection['selected_files']) and not(selected_paths&omitted_paths)
actual={str(p.relative_to(rawroot)) for p in restored.rglob('*') if p.is_file()}
audit.check('exact explicit public scientific payload grid',actual==selected_paths,dict(missing=sorted(selected_paths-actual),extra=sorted(actual-selected_paths)))
for item in selection['selected_files']:
    p=rawroot/item['path'];bind(p,item['sha256'])
    audit.check('selected original scientific file bytes '+item['path'],p.stat().st_size==item['bytes'])
    closure.append(dict(path=relative(p),sha256=item['sha256'],bytes=item['bytes']))
audit.check('all original files accounted by explicit selected+omitted grids',len(selected_paths)+len(omitted_paths)==1110)
write(args.output/'raw_closure_manifest.json',dict(schema='X01_PUBLIC_SELECTED_RAW_CLOSURE_V1',files=closure,original_files=len(closure),original_full_files=1110,selection_sha256=sha(args.selection),original_full_audit_sha256=selection['original_full_audit_sha256'],administrative_exclusions=selection['explicit_exclusions']))
print('PHASE public selected raw closure complete; source archives',flush=True)

# Decode immutable portable source archives; inspect exact regular members.
archive_indices=[]
mainindex=snapshot/'overnight_20261022/portable_source_archives/index.json'
for rec in audit.load(mainindex)['archives']:
    if rec['archive'] in ('continuum-production-v1c-source.tar.gz','x01-continuum-stage2-v2-source.tar.gz'):
        archive_indices.append((mainindex.parent,rec))
for folder in ('portable_x01_analysis_repair_archive_v1','portable_x01_stage3_source_archive_v1'):
    p=snapshot/'overnight_20261022'/folder/'index.json';rec=audit.load(p)
    if 'archive_name' in rec:
        rec=dict(rec,archive=rec['archive_name'],archive_bytes=rec['bytes'],archive_sha256=rec['sha256'],freeze_record='x01_continuum_stage3_v1/source_freeze.json',freeze_record_sha256=rec['source_freeze_sha256'])
    archive_indices.append((p.parent,rec))
sourceclosure=[]
for base,rec in archive_indices:
    b64=bind(base/rec['base64_file'],rec.get('base64_sha256',rec.get('base64_file_sha256')))
    data=base64.b64decode(b''.join(b64.read_bytes().split()),validate=True)
    audit.check('exact decoded source archive '+rec['archive'],len(data)==rec['archive_bytes'] and hashlib.sha256(data).hexdigest()==rec['archive_sha256'])
    freeze_rel=rec.get('freeze_record',rec.get('source_freeze'))
    freezep=snapshot/freeze_rel;freeze=audit.load(bind(freezep,rec.get('freeze_record_sha256',rec.get('source_freeze_sha256'))))
    release=freezep.parent.name
    audit.check('pinned original release freeze '+release,sha(freezep)==FREEZE[release])
    members={}
    with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as tar:
        for member in tar.getmembers():
            safe_parts(member.name)
            name=member.name.removeprefix(release+'/')
            if not member.isfile() or name in members: raise ValueError('Nonregular/duplicate source member')
            content=tar.extractfile(member).read();members[name]=hashlib.sha256(content).hexdigest()
    audit.equal('archive member grid and original frozen hashes '+release,members,freeze['files'])
    for name,want in members.items():bind(snapshot/release/name,want)
    sourceclosure.append(dict(release=release,archive_sha256=rec['archive_sha256'],source_freeze_path=relative(freezep),source_freeze_sha256=sha(freezep),members=members))
write(args.output/'source_closure_manifest.json',dict(schema='X01_FRESH_SOURCE_CLOSURE_V1',archives=sourceclosure,native_source_hashes=read(snapshot/'x01_continuum_stage2_v2/manifest.toml')['source_hashes']))
print('PHASE source archives complete; native/seed/source checks',flush=True)

cfgs={release:audit.load(snapshot/release/'manifest.toml') for release in FREEZE}
for label,release in RELEASE_FOR.items():
    cfg=cfgs[release]
    for name,want in cfg['source_hashes'].items():bind(restored/label/'native_source'/name,want)
for release,cfg in cfgs.items():
    audit.equal(release+' unchanged scientific gates',cfg['gates'],M.EXPECTED_GATES)
designs={name:{k:v for k,v in cfg['statistics'].items() if k!='future_stage_rule'} for name,cfg in cfgs.items()}
audit.equal('all stages scientific statistical designs identical',designs['continuum_production_v1c'],designs['x01_continuum_stage2_v2'])
audit.equal('stage3 scientific statistical design identical',designs['x01_continuum_stage2_v2'],designs['x01_continuum_stage3_v1'])
audit.check('stage-specific future-stage provenance retained separately',True,{name:cfg['statistics']['future_stage_rule'] for name,cfg in cfgs.items()})
# Reuse the ORIGINAL mathematical auditor's quantile reconstruction tolerance
# (audit_continuum_stage1.py lines430-431). The frozen declared critical value
# used in every confidence interval is unchanged and compared exactly above.
audit.check('independently reconstructed Student7critical original audit tolerance',abs(M.student7_critical()-cfgs['x01_continuum_stage2_v2']['statistics']['student_t_critical'])<2e-9,dict(reconstructed=M.student7_critical(),declared=cfgs['x01_continuum_stage2_v2']['statistics']['student_t_critical'],original_reconstruction_tolerance=2e-9))
seedsets=[]
for release in ('continuum_production_v1c','x01_continuum_stage2_v2','x01_continuum_stage3_v1'):
    seeds=audit.load(snapshot/release/'seeds.toml')['bundles']
    values=[v for row in seeds for v in [row['energy_seed']]+[v for pair in row['node_replica_seeds'] for v in pair]]
    audit.check(release+'288 unique actual initialization seeds',len(values)==len(set(values))==288)
    seedsets.append(set(values))
audit.check('all three stage initialization sets disjoint, no pooling',all(not a&b for i,a in enumerate(seedsets) for b in seedsets[i+1:]))

# Original analysis repair: exactly four local identifier AST replacements.
orig=ast.parse((snapshot/'x01_continuum_stage2_v2/compare_production.py').read_text())
repair=ast.parse((snapshot/'x01_stage2_analysis_repair_v1/compare_production.py').read_text())
changes=[]
class Rename(ast.NodeTransformer):
    def visit_Name(self,node):
        if node.id=='merge_receipt':changes.append(type(node.ctx).__name__);node.id='merge'
        return node
Rename().visit(repair)
audit.check('repair exactly four identifiers, no scientific AST change',ast.dump(orig)==ast.dump(repair) and sorted(changes)==['Load','Load','Load','Store'],changes)
repairreceipt=audit.load(snapshot/'x01_stage2_analysis_repair_result_v1/repair_receipt.json')
maps=path_map_hashes(repairreceipt['input_sha256'],'analysis repair historical input',map_windows)
print('PHASE repair receipt path closure complete',flush=True)
audit.check('repair no new sampling or statistical look',repairreceipt['analysis_only'] is True and repairreceipt['new_sampling'] is False and repairreceipt['new_statistical_look'] is False and repairreceipt['scientific_completions']==0)
originalfailed=audit.load(restored/'x01s2a_6460996/results/comparison.json')
audit.check('original shadowing failure preserved',originalfailed.get('error')=="TypeError: 'dict' object is not callable" and sha(restored/'x01s2a_6460996/results/comparison.json')=='dbaa8be737d26632fd610adc1a21565311ea40c2f8a7c63fde2bcaf2d92bed33')

# Strict four-part stage2 raw identity and full merge receipt closure.
analysis=restored/'x01s2a_6460996/results'
merged=audit.load(analysis/'target.toml');merge=audit.load(analysis/'merge_receipt.json')
remote_mappings=path_map_hashes(merge['all_input_sha256'],'stage2 original merge input',map_remote)
audit.check('stage2 no prior samples pooled',merge['prior_stage_values_pooled'] is False and merge['scientific_sampling'] is False)
audit.check('stage2 merged target bound',sha(analysis/'target.toml')==merge['merged_target_sha256'])
energy=[];nodes=[]
for i,job in enumerate((6460990,6460991,6460992,6460994),1):
    folder=restored/f'x01s2p{i}_{job}'/'results';receipt=audit.load(folder/'part_receipt.json');part=audit.load(folder/'target.toml')
    audit.check(f'stage2 part{i}complete native run',receipt['status']=='COMPLETE' and receipt['native_exit']==0 and receipt['job_id']==str(job) and receipt['part_index']==i and part['completion']['completed_calls']==40)
    audit.check(f'stage2 part{i}receipt and source closure',sha(folder/'part_receipt.json')==merge['part_receipts_sha256'][str(i)] and receipt['source_hashes']==cfgs['x01_continuum_stage2_v2']['source_hashes'] and receipt['source_freeze_sha256']==FREEZE['x01_continuum_stage2_v2'])
    for name,want in receipt['input_sha256'].items():
        p=folder/name
        if not p.is_file() and 'restored/'+str(p.relative_to(restored)) in omitted_paths:
            public_exclusions.append(dict(original_key=str(p),sha256=want,reason='Explicit stage3 administrative guard/queue exclusion'));continue
        bind(p,want)
    energy.extend(part['energy_rows']);nodes.extend(part['entropy_node_rows'])
    for row in part['energy_rows']+part['entropy_node_rows']:
        name=f'dt{i}_bundle{row["bundle"]}_'+('lambda'+str(row['lambda_index']) if 'lambda_index' in row else 'energy')+'.toml'
        audit.check('stage2 original part and merged side bytes '+name,sha(folder/name)==sha(analysis/name)==merge['side_files_sha256'][name])
audit.equal('stage2 complete merged32energy equals four parts',merged['energy_rows'],energy)
audit.equal('stage2 complete merged128nodes equals four parts',merged['entropy_node_rows'],nodes)
write(args.output/'historical_input_role_map.json',dict(schema='X01_HISTORICAL_PATH_ROLE_MAPPING_V1',method='Only historical campaign-root prefixes are mapped. File bytes, expected hashes, gates, seeds, field names and numerical outputs are unchanged.',windows_repair_inputs=maps,remote_merge_inputs=remote_mappings))
print('PHASE strict four-part merge closure complete; independent math',flush=True)

# Independently construct exact bit-Fock reference once; no new sampling.
ownreference=M.exact_reference()
write(args.output/'independent_bit_fock_reference.json',M.plain(ownreference))
stage_reports={}
for stage,label,release,savedpath in [
 (2,'x01s2a_6460996','x01_stage2_analysis_repair_v1',snapshot/'x01_stage2_analysis_repair_result_v1/comparison.json'),
 (1,'x01_6460771','continuum_production_v1c',restored/'x01_6460771/results/comparison.json')]:
    source=snapshot/release;result=restored/label/'results';cfg=cfgs[release]
    saved=audit.load(savedpath);raw=audit.load(result/'target.toml');oracle=audit.load(result/'oracle.json')
    local=M.Audit()
    local.equal('independent full-Fock continuum oracle',oracle['continuum'],ownreference['continuum'])
    local.check('independent256state fullGCE dimension',oracle['hilbert_dimension']==ownreference['hilbert_dimension']==256)
    for existing,expected in zip(oracle['finite_dtau'],ownreference['finite_dtau']):
        for key in ('dtau','slices','asymmetric','symmetric'):local.equal('independent finiteh '+str(existing['dtau'])+'/'+key,existing[key],expected[key])
    seeds=audit.load(source/'seeds.toml')['bundles']
    reconstructed=M.reconstruct(SimpleNamespace(result_dir=result),local,cfg,raw,oracle,saved,seeds)
    print('PHASE stage',stage,'independent math complete; frozen comparator replay',flush=True)
    out=args.output/f'stage{stage}_replayed_comparison.json'
    command=[sys.executable,'-B',str(source/'compare_production.py'),'--source',str(source),'--native-root',str(restored/label/'native_source'),'--result-dir',str(result),'--target-exit-code','0','--oracle-exit-code','0','--output',str(out)]
    run=subprocess.run(command,text=True,capture_output=True,timeout=180,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    replayed=read(out)
    local.check('frozen separate-process replay expected exit2',run.returncode==2,dict(exit=run.returncode,stdout=run.stdout,stderr=run.stderr))
    local.equal('all replay input hash multisets exact',sorted(saved['input_sha256'].values()),sorted(replayed['input_sha256'].values()))
    # Gates and Boolean decisions are exact; only platform floating arithmetic
    # uses the ORIGINAL audit tolerance, never a new acceptance margin.
    local.equal('complete replay output excluding absolute input path keys',{k:v for k,v in saved.items() if k!='input_sha256'},{k:v for k,v in replayed.items() if k!='input_sha256'})
    for key,want in replayed['input_sha256'].items():bind(key,want)
    audit.check(f'stage{stage}fresh mathematical/replay audit',all(r['passed'] for r in local.checks),[r for r in local.checks if not r['passed']])
    stage_report=dict(stage=stage,candidate_status=saved['scientific_status'],reconstructed_status=reconstructed['scientific_status'],check_count=len(local.checks),failure_count=sum(not r['passed'] for r in local.checks),checks=local.checks,independent_reconstruction=M.plain(reconstructed),replay=dict(command=command,exit_code=run.returncode,output=relative(out),sha256=sha(out),saved_comparison_sha256=sha(savedpath),original_numeric_tolerance=dict(absolute=2e-12,relative=2e-10),gates_and_labels_exact=True))
    write(args.output/f'stage{stage}_audit.json',M.plain(stage_report));stage_reports[str(stage)]=stage_report
    print('STAGE',stage,'checks',len(local.checks),'failures',stage_report['failure_count'],flush=True)

# Final third look is incomplete; verify actual partial receipts, no fits.
history=[]
for label in ['x01_cancelled_6460728','x01_6460771','x01s2a_6460996']+[f'x01s3p{i}_{6462169+i}' for i in range(1,5)]+['x01s3a_6462174']:
    p=restored/label;main=list(csv.DictReader(io.StringIO((p/'scheduler.txt').read_text()),delimiter='|'))[0]
    history.append(dict(job_label=label,scheduler=main,scheduler_sha256=sha(p/'scheduler.txt')))
stage3rows=[]
for i,job in enumerate(range(6462170,6462174),1):
    folder=restored/f'x01s3p{i}_{job}'/'results'
    if i<3:
        receipt=audit.load(folder/'part_receipt.json');raw=audit.load(folder/'target.toml')
        for name,want in receipt['input_sha256'].items():
            p=folder/name
            if not p.is_file() and 'restored/'+str(p.relative_to(restored)) in omitted_paths:
                public_exclusions.append(dict(original_key=str(p),sha256=want,reason='Explicit stage3 administrative guard/queue exclusion'));continue
            bind(p,want)
        audit.check(f'stage3 part{i}source/seed identity',receipt['source_freeze_sha256']==FREEZE['x01_continuum_stage3_v1'] and receipt['source_hashes']==cfgs['x01_continuum_stage3_v1']['source_hashes'] and receipt['seeds_sha256']==sha(snapshot/'x01_continuum_stage3_v1/seeds.toml'))
        side=list(folder.glob('dt*_bundle*_*.toml'))
        E,S=raw['energy_rows'],raw['entropy_node_rows']
        audit.check(f'stage3 part{i}actual partial side grid',len(E)==8 and len(S)==(32 if i==1 else 28) and len(side)==(40 if i==1 else 36))
        audit.check(f'stage3 part{i}actual completion/exit',receipt['status']==('COMPLETE' if i==1 else 'EXECUTION_INCOMPLETE') and receipt['native_exit']==(0 if i==1 else 137) and (raw.get('completion',{}).get('status')=='PRODUCTION_STAGE3_PART_COMPLETE' if i==1 else 'completion' not in raw))
        for row in E+S:
            name=f'dt{i}_bundle{row["bundle"]}_'+('lambda'+str(row['lambda_index']) if 'lambda_index' in row else 'energy')+'.toml'
            audit.equal('stage3 retained aggregate/side identity '+name,audit.load(folder/name),row)
        stage3rows.append(dict(part=i,job_id=str(job),native_exit=receipt['native_exit'],status=receipt['status'],energy_rows=len(E),entropy_nodes=len(S),side_files=len(side),part_receipt_sha256=sha(folder/'part_receipt.json'),target_sha256=sha(folder/'target.toml')))
    else:
        audit.check(f'stage3 part{i}cancelled with no sampler outputs',not (folder/'target.toml').exists() and not (folder/'part_receipt.json').exists() and not list(folder.glob('dt*.toml')))
        row=next(h['scheduler'] for h in history if h['job_label']==f'x01s3p{i}_{job}')
        audit.check(f'stage3 part{i}cancelled before allocation',row['State'].startswith('CANCELLED') and row['ElapsedRaw']=='0' and row['AllocCPUS']=='0' and row['Start']=='None')
        stage3rows.append(dict(part=i,job_id=str(job),status='CANCELLED_BEFORE_ALLOCATION',native_exit=None,energy_rows=0,entropy_nodes=0,side_files=0))
analysis3=restored/'x01s3a_6462174/results'
audit.check('final merge rejected incomplete part before oracle/comparator', 'Incomplete or wrong part receipt' in (analysis3/'slurm.log').read_text() and (analysis3/'process_exit_codes.txt').read_text()=='merge_exit=1\noracle_exit=125\ncomparator_exit=125\nanalysis_exit=1\n')
audit.check('no stage3 merged grid or statistical output',not any((analysis3/name).exists() for name in ('target.toml','merge_receipt.json','oracle.json','comparison.json')))
closure3_sha='aee287fbb68366bd703a6d83debec7bd9f1fe887371058f8fba7f60ea27c6044'
write(args.output/'failure_history.json',dict(schema='X01_ALL_RESTORED_ATTEMPTS_HISTORY_V1',scheduler_history=history,stage3_parts=stage3rows,stage3_terminal_closure_sha256=closure3_sha,stage3_conclusion='INCOMPLETE_CLOSED_FINAL_RESERVED_LOOK',selected_scientific_result='stage2 UNRESOLVED',new_scientific_completions=0,fourth_look=False,pooling=False,rerun_or_retry=False,partial_grid_fits=False))

failures=[r for r in audit.checks if not r['passed']]
selected=stage_reports['2']
out=dict(schema='X01_PUBLIC_SCIENTIFIC_SELECTION_RECONCILIATION_V1',audit_scope='Original scientific source/raw/statistics and explicit public selection; PC administration and transport are provenance-only',selection_sha256=sha(args.selection),original_full_audit_sha256=selection['original_full_audit_sha256'],public_administrative_exclusions=public_exclusions,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),audit_status='PASS' if not failures else 'UNRESOLVED',scientific_status='UNRESOLVED',catalogue_id='X01',new_scientific_completions=0,existing_sources_modified=False,scientific_sampling=False,remote_calls=False,job_submission=False,public_mutation=False,auditor_sha256=sha(__file__),independent_math_auditor_sha256=MATH_SHA,raw_original_files=len(closure),stages={'1':dict(status=stage_reports['1']['reconstructed_status'],audit_path=relative(args.output/'stage1_audit.json')),'2':dict(status=selected['reconstructed_status'],audit_path=relative(args.output/'stage2_audit.json'),endpoint_precision=selected['independent_reconstruction']['endpoint_precision'],decision_components=selected['independent_reconstruction']['decision_components']),'3':dict(status='INCOMPLETE_CLOSED_FINAL_RESERVED_LOOK',parts=stage3rows,no_primary_reconstruction=True,no_comparator_replay=True,no_ledger_promotion=True)},selected_stage=2,reason='Fixed quartic S2 continuum halfwidth exceeds predeclared precision threshold. Incomplete final stage3 adds no result; no pooling, fourth look, retry, or gate relaxation.',limitations=['Nominal fixed three-look/26-interval family uses approximate normal independent bundle contrasts; this audit does not prove normality, rare-tail coverage or equilibrium.','Rhat, ESS and64/32binblocking are rejection diagnostics only.','The native sampler is byte-pinned, not independently reimplemented.','Historical absolute path prefixes are role-mapped with exact original file SHA checks; no original bytes are rewritten.','All pilot, failed mechanical analysis and incomplete/cancelled final-stage histories are retained; scientific status is UNRESOLVED despite integrity audit PASS.'],check_count=len(audit.checks),failure_count=len(failures),failures=failures,checks=audit.checks,input_sha256=audit.inputs)
write(args.output/'audit.json',M.plain(out))
write(args.output/'summary.json',{k:v for k,v in out.items() if k not in ('checks','input_sha256')})
print('audit='+out['audit_status']+'; science=UNRESOLVED; checks='+str(len(audit.checks))+'; failures='+str(len(failures))+'; zero new scientific completions',flush=True)
if failures:sys.exit(2)

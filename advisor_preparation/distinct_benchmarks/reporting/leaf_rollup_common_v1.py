#!/usr/bin/env python3
"""Read-only leaf-report helpers. No scientific computation or cluster calls."""
import copy,csv,hashlib,json,math,shutil,subprocess,sys,tarfile
from datetime import datetime,timezone
from pathlib import Path,PurePosixPath

ROOT=Path(__file__).resolve().parents[1]
PREP=ROOT.parent
IDS={f'{p}{i:02d}' for p in ('A','J','P','X') for i in range(1,26)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read(p):
    return json.loads(Path(p).read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON:'+x)))
def dump(p,data):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def relative(s):
    p=PurePosixPath(s);require(not p.is_absolute() and '..' not in p.parts and str(p) not in ('','.'),'Unsafe relative path')
    return Path(*p.parts)
def canonical(r):return hashlib.sha256(json.dumps(r,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def tally(rows):
    require(len(rows)==100 and {r['id'] for r in rows}==IDS,'Not exactly100uniqueIDs')
    c={s:sum(r['status']==s for r in rows) for s in ('PASS','FAIL','UNRESOLVED','NOT_RUN')};require(sum(c.values())==100,'Invalid status');return c
def finite(x):
    if isinstance(x,dict):return all(finite(v) for v in x.values())
    if isinstance(x,list):return all(finite(v) for v in x)
    return not isinstance(x,float) or math.isfinite(x)
def equal(a,b):
    changes=[]
    def visit(x,y,path):
        require(type(x) is type(y),'Replay type:'+path)
        if isinstance(x,dict):
            require(x.keys()==y.keys(),'Replay keys:'+path)
            for k in x:visit(x[k],y[k],path+'/'+k)
        elif isinstance(x,list):
            require(len(x)==len(y),'Replay length:'+path)
            for i,(u,v) in enumerate(zip(x,y)):visit(u,v,path+'/'+str(i))
        elif isinstance(x,float):
            d=abs(x-y);require(math.isfinite(x) and math.isfinite(y) and d<=64*sys.float_info.epsilon*max(1.,abs(x),abs(y)),'Replay float:'+path)
            if d:changes.append(dict(path=path,absolute_difference=d,reported=x,replayed=y))
        else:require(x==y,'Replay value:'+path)
    visit(a,b,'');return dict(nonfloats_exact=True,float_allowance='64*binary64_epsilon*max(1,abs(a),abs(b))',roundoff_differences=changes)

class Evidence:
    def __init__(self,out,manifest):self.out=out;self.manifest=manifest
    def add(self,role,path,expected=None):
        path=Path(path);require(path.is_file() and not path.is_symlink(),'Missing/nonregular evidence:'+str(path));path=path.resolve()
        digest=sha(path);require(expected is None or digest==expected,'SHA differs:'+role)
        require(role not in self.manifest,'Duplicate role:'+role)
        rel=Path('inputs')/role/path.name;dst=self.out/rel;require(not dst.exists(),'Duplicate destination')
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dst)
        try:origin=path.relative_to(PREP.resolve()).as_posix()
        except ValueError:origin='external/'+role+'/'+path.name
        self.manifest[role]=dict(original_path=origin,evidence_copy=rel.as_posix(),sha256=digest,bytes=dst.stat().st_size)
        return dst

def copy_base(path,out,ledger_sha,manifest_sha,passed,label):
    path=Path(path);require(sha(path/'ledger_100.json')==ledger_sha and sha(path/'input_manifest.json')==manifest_sha,'Immutable base identity differs')
    base=read(path/'ledger_100.json');manifest=read(path/'input_manifest.json')
    require(base['input_manifest']==manifest and tally(base['rows'])==base['counts']=={'PASS':passed,'FAIL':0,'UNRESOLVED':0,'NOT_RUN':100-passed},'Base manifest/counts inconsistent')
    for role,r in manifest.items():
        p=path/relative(r['evidence_copy']);require(p.is_file() and not p.is_symlink() and sha(p)==r['sha256'] and p.stat().st_size==r['bytes'],'Base input differs:'+role)
    for rel,h in base['analysis_receipts'].items():require(sha(path/relative(rel))==h,'Base analysis receipt differs:'+rel)
    closure=[]
    for p in sorted(path.rglob('*')):
        require(not p.is_symlink(),'Base contains symlink')
        if not p.is_file():continue
        rel=p.relative_to(path);dest=rel if len(rel.parts)>1 else Path('inputs')/(label+'_metadata')/rel
        q=out/dest;require(not q.exists(),'Base copy collision');q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        closure.append(dict(base_relative_path=rel.as_posix(),evidence_copy=dest.as_posix(),sha256=sha(p),bytes=p.stat().st_size))
    for p in sorted(path.iterdir()):
        if p.is_file():
            rel=Path('inputs')/(label+'_metadata')/p.name
            manifest[label+'_metadata_'+p.name.replace('.','_')]=dict(original_path=p.resolve().relative_to(PREP.resolve()).as_posix(),evidence_copy=rel.as_posix(),sha256=sha(p),bytes=p.stat().st_size)
    dump(out/'analysis'/f'{label}_file_closure.json',closure)
    return base,manifest,closure

def frozen_release(ev,role,freeze_path,freeze_sha,archive_path,archive_sha,root_name,ids):
    f=read(ev.add(role+'_freeze',freeze_path,freeze_sha));archive=ev.add(role+'_archive',archive_path,archive_sha)
    require(f['archive_sha256']==archive_sha and f['distinct_ids']==list(ids),'Release identity/IDs')
    seen={};root=ev.out/'frozen_sources'/role
    with tarfile.open(archive) as tar:
        for m in tar:
            if m.isdir():continue
            require(m.isfile(),'Nonregular archive member');r=relative(m.name);k=r.as_posix()
            require(r.parts[0]==root_name and k in f['files'] and k not in seen,'Unexpected/duplicate frozen member')
            b=tar.extractfile(m).read();h=hashlib.sha256(b).hexdigest();require(h==f['files'][k],'Frozen member hash:'+k)
            p=root/r;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b);seen[k]=h
    require(seen==f['files'],'Incomplete release closure')
    return root/root_name,f

def independent(ev,role,audit_path,auditor_path,replay_path,ids,raw,archive_sha,schema):
    audit=read(ev.add(role+'_independent_audit',audit_path));auditor=ev.add(role+'_independent_auditor',auditor_path)
    replay=ev.add(role+'_independent_replay',replay_path)
    require(audit['schema']==schema and audit['status']=='PASS' and audit['catalogue_ids']==list(ids),'Independent audit scope/status')
    require(audit['auditor_sha256']==sha(auditor) and audit['source_archive_sha256']==archive_sha,'Independent auditor/release binding')
    hashes=audit['input_sha256']
    for name,path in {**raw,'independent_replay':replay}.items():
        require(sha(path) in [h for k,h in hashes.items() if Path(k).name==Path(path).name or k==name],'Independent audit missing binding:'+name)
    return audit,read(replay)

def replay_comparator(ev,role,command,reported,independent_replay):
    output=ev.out/'analysis'/f'{role}_comparator_replay.json';cmd=command(output)
    p=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
    dump(ev.out/'analysis'/f'{role}_replay_process.json',dict(command=cmd,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,physics_executed=False))
    require(p.returncode==0 and output.exists(),'Frozen comparator replay failed:'+role)
    data=read(output);agreement=equal(reported,data);equal(independent_replay,data)
    return data,agreement

def finish(out,base,ev,rows,changed,expected_passed,closure,tool,scope,audits,base_label,base_hash):
    before={r['id']:r for r in base['rows']};preserved={}
    for r in rows:
        if r['id'] not in changed:
            require(r==before[r['id']],'Base row changed:'+r['id']);preserved[r['id']]=canonical(r)
        else:require(before[r['id']]['status']=='NOT_RUN' and r['status']=='PASS','Not a new complete passing idea')
    for r in base['rows']:
        if r['status']=='PASS':require(r==next(x for x in rows if x['id']==r['id']),'Prior PASS changed')
    c=tally(rows);require(c=={'PASS':expected_passed,'FAIL':0,'UNRESOLVED':0,'NOT_RUN':100-expected_passed},'Unexpected completion count')
    receipt={p.relative_to(out).as_posix():sha(p) for p in sorted((out/'analysis').rglob('*')) if p.is_file()}
    data=dict(generated_utc=datetime.now(timezone.utc).isoformat(),tool_sha256=sha(tool),counts=c,
        category_counts={category:{s:sum(r['category']==category and r['status']==s for r in rows) for s in c} for category in sorted({r['category'] for r in rows})},
        counting_rule=base['counting_rule'],scope_limit=scope,base_ledger_sha256=base_hash,
        base_preservation=dict(label=base_label,all_prior_PASS_rows_unchanged=True,all_other_row_objects_unchanged=True,canonical_preserved_row_sha256=preserved,all_base_files_copied=len(closure)),
        tranche_audits=audits,input_manifest=ev.manifest,analysis_receipts=receipt,rows=rows)
    dump(out/'input_manifest.json',ev.manifest);dump(out/'ledger_100.json',data)
    with (out/'ledger_100.csv').open('w',newline='') as f:
        fields=['id','category','title','status','reported_status','fixtures','comparisons','comparison_unit','reason','metrics'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(r[k],sort_keys=True) if k=='metrics' else r.get(k,'') for k in fields})
    (out/'README.md').write_text(f'# {expected_passed} distinct externally anchored benchmarks passed\n\nExactly100catalogueIDs: {expected_passed}PASS,0FAIL,0UNRESOLVED,{100-expected_passed}NOT_RUN. Parameters, matrix entries and gates are not additional ideas.\n\nEvery base file and all other row objects are preserved, including all historical failures. New completions require frozen source/archive closure, full receipt grids, actual target source/runtime checks, independent audit and a fresh replay of each exact frozen comparator. Original/replayed comparisons and all authored source are retained. Gate metrics keep their actual inequality directions and units; positivity and lower-bound checks are not absolute errors.\n\n'+scope+'\n\nSee [ledger](ledger_100.json), [input hashes](input_manifest.json) and analysis/. No physics or cluster job is run by this reporter.\n')
    print(json.dumps(dict(counts=c,output=str(out),physics_executed=False),indent=2));return data

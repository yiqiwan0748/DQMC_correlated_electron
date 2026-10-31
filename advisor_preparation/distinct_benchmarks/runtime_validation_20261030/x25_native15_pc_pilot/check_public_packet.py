"""Check a scalar pilot packet; no model, eigensolve or reference execution."""
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import stat

NAMES={'README.md','AUDIT_SCOPE.md','results.json','native_criteria.json','sources.json',
       'scientific_cases.csv','component_timings.csv','check_public_packet.py','MANIFEST.json'}
FILE_LIMIT=1024*1024
TOTAL_LIMIT=4*1024*1024


def require(value,message):
    if not value:raise ValueError(message)


def read_bounded(path):
    before=path.lstat()
    require(stat.S_ISREG(before.st_mode)and not path.is_symlink()and before.st_size<=FILE_LIMIT,'Bounded regular packet file required')
    fd=os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK)
    try:
        start=os.fstat(fd)
        require(stat.S_ISREG(start.st_mode)and start.st_size<=FILE_LIMIT and (start.st_dev,start.st_ino)==(before.st_dev,before.st_ino),'Packet file identity differs')
        data=bytearray()
        while True:
            block=os.read(fd,min(65536,FILE_LIMIT+1-len(data)))
            if not block:break
            data.extend(block);require(len(data)<=FILE_LIMIT,'Packet file exceeds bound')
        end=os.fstat(fd)
        require((start.st_dev,start.st_ino,start.st_size,start.st_mtime_ns,start.st_ctime_ns)==(end.st_dev,end.st_ino,end.st_size,end.st_mtime_ns,end.st_ctime_ns)and len(data)==start.st_size,'Packet file changed while read')
        return bytes(data)
    finally:os.close(fd)


def load(raw):
    def pairs(items):
        result={}
        for k,v in items:
            require(k not in result,'Duplicate JSON key')
            result[k]=v
        return result
    def bad(_):raise ValueError('Nonfinite JSON')
    value=json.loads(raw,object_pairs_hook=pairs,parse_constant=bad)
    def walk(x):
        if type(x)is float:require(math.isfinite(x),'Nonfinite scalar')
        elif type(x)is list:
            for y in x:walk(y)
        elif type(x)is dict:
            for y in x.values():walk(y)
    walk(value);return value


def main():
    here=Path(__file__).parent
    require({p.name for p in here.iterdir()}==NAMES,'Exact public file closure differs')
    stats=[p.lstat()for p in here.iterdir()]
    require(all(stat.S_ISREG(s.st_mode)and s.st_size<=FILE_LIMIT for s in stats)and sum(s.st_size for s in stats)<=TOTAL_LIMIT,'Packet size/type bound differs')
    bodies={};total=0
    for name in sorted(NAMES):
        body=read_bounded(here/name);total+=len(body);require(total<=TOTAL_LIMIT,'Actual packet total exceeds bound');bodies[name]=body
    manifest=load(bodies['MANIFEST.json'])
    require(set(manifest['files'])==NAMES-{'MANIFEST.json'},'Manifest coverage differs')
    for name,item in manifest['files'].items():
        body=bodies[name]
        require(item=={'bytes':len(body),'sha256':hashlib.sha256(body).hexdigest()},'Packet hash differs: '+name)
    d=load(bodies['results.json'])
    require(set(d)=={'schema','status','scientific_credit','benchmark_ideas_added','catalogue_status','full_X25_status','independent_QuSpin_reference','external_ground_state_or_response_certificates','campaign_ledger','model','native_solver','case_count','cases','timing','resources','saved_audit_scope'},'Public scientific field closure differs')
    require(d['schema']=='X25_NATIVE15_PUBLIC_PILOT_V1'and d['status']=='PILOT_SAVED_AUDITS_PASS_ZERO_CREDIT','Pilot scope differs')
    require(type(d['scientific_credit'])is int and d['scientific_credit']==0 and d['benchmark_ideas_added']==0,'Pilot credit differs')
    require(all(d[k]=='NOT_RUN'for k in ['catalogue_status','full_X25_status','independent_QuSpin_reference','external_ground_state_or_response_certificates']),'Unperformed work relabeled')
    require(d['campaign_ledger']==dict(PASS=50,FAIL=0,UNRESOLVED=2,NOT_RUN=48),'Ledger differs')
    expected=[(f'{i+1:02d}',h,n,4-n,(70,448,784,448,70)[n])for i,(h,n)in enumerate((h,n)for h in(0.,-.02,.02)for n in range(5))]
    cases=d['cases'];require(len(cases)==d['case_count']==15,'All15 cases required')
    for row,case in zip(cases,expected):
        require(tuple(row[k]for k in ['case_id','h','N_up','N_down','dimension'])==case,'Case changed/reordered')
        require(row['native_method']=='lanczos'and row['native_converged']is True,'Native method/convergence differs')
    science=list(csv.DictReader(io.StringIO(bodies['scientific_cases.csv'].decode('utf-8'),newline='')))
    require(len(science)==15,'Scientific CSV incomplete')
    for a,b in zip(science,cases):
        for k,v in a.items():
            observed=v if k=='case_id'else int(v)if k in ['N_up','N_down','dimension','iterations']else float(v)
            require(observed==b[k],'Scientific scalar differs: '+k)
    timings=list(csv.DictReader(io.StringIO(bodies['component_timings.csv'].decode('utf-8'),newline='')))
    require(len(timings)==len(d['timing']['per_case'])==15,'Timing CSV incomplete')
    for a,b,c in zip(timings,d['timing']['per_case'],cases):
        require(a['case_id']==b['case_id']==c['case_id'],'Timing case differs')
        for k,v in a.items():require((v if k=='case_id'else float(v))==b[k],'Timing scalar differs: '+k)
    require(d['timing']['native_ground_API_is_Lanczos_only']is False and d['timing']['compilation_separately_measured']is False,'Timing scope changed')
    require(d['resources']['per_Julia_or_per_case_RSS']=='NOT_RECORDED','Memory scope changed')
    criteria=load(bodies['native_criteria.json'])
    require(criteria['scientific_credit']==0 and criteria['external_comparison']=='NOT_RUN'and all(r['failed']==0 for r in criteria['criteria']),'Criterion scope differs')
    require({r['kind']for r in criteria['criteria']}=={'saved_array_consistency_diagnostic','frozen_native_gate'},'Criterion kinds differ')
    source=load(bodies['sources.json'])
    require(all(set(r)=={'logical_snapshot_filename','role','bytes','sha256'}and not r['logical_snapshot_filename'].startswith('/')for r in source['references']),'Source reference fields differ')
    print(json.dumps({'status':'PASS_PACKET_HASHES_AND_SCALAR_CONSISTENCY_ONLY','case_count':15,'scientific_credit':0,'full_X25_status':'NOT_RUN','external_reference':'NOT_RUN'},sort_keys=True))


if __name__=='__main__':main()

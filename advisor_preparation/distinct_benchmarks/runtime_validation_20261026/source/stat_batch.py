#!/usr/bin/env python3
"""Independent full-width Python/C metadata, exact decimal strings; no numerics."""
import argparse, hashlib, json, os, subprocess, time, tomllib
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def statrow(p):
    x=os.stat(p)
    return {k:str(v) for k,v in dict(inode=x.st_ino,device=x.st_dev,major=os.major(x.st_dev),minor=os.minor(x.st_dev),bytes=x.st_size,mtime_ns=x.st_mtime_ns,ctime_ns=x.st_ctime_ns).items()}

def main():
    a=argparse.ArgumentParser();a.add_argument('--paths',type=Path,required=True);a.add_argument('--helper',type=Path,required=True);a.add_argument('--helper-sha256',required=True);a.add_argument('--output',type=Path,required=True);x=a.parse_args()
    if x.output.exists():raise RuntimeError('Fresh metadata output required')
    report=dict(schema='P17_INDEPENDENT_FULLWIDTH_PYTHON_C_STAT_V1',status='RUNNING',scientific_completions=0,rows=[],helper_sha256=x.helper_sha256)
    def save():x.output.write_text(json.dumps(report,indent=2)+'\n')
    def cstat(p):
        if sha(x.helper)!=x.helper_sha256:raise RuntimeError('Compiled helper bytes changed')
        q=subprocess.run([str(x.helper),str(p)],capture_output=True,text=True,timeout=5)
        return dict(exit_code=q.returncode,result=json.loads(q.stdout),stderr=q.stderr)
    tools=tomllib.loads((Path(__file__).resolve().parent/'metadata/tools.toml').read_text())['tools']
    paths=tomllib.loads(x.paths.read_text())['paths']
    def gnu_batch():
        tool=tools['stat'];q=Path(tool['path'])
        if str(q.resolve())!=tool['realpath'] or q.stat().st_size!=tool['bytes'] or sha(q)!=tool['sha256']:raise RuntimeError('Fixed GNU stat changed')
        args=[tool['path'],'--dereference','--printf=%d|%i|%s\n','--']+[str(Path(v).resolve())for v in paths]
        got=subprocess.run(args,capture_output=True,text=True,timeout=15)
        lines=got.stdout.splitlines()
        parsed=[line.split('|')for line in lines]
        valid=got.returncode==0 and len(parsed)==len(paths)and all(len(v)==3 and all(z.isdecimal()and 0<=int(z)<2**64 for z in v)for v in parsed)
        after=sha(q)
        if after!=tool['sha256']or str(q.resolve())!=tool['realpath']or q.stat().st_size!=tool['bytes']:raise RuntimeError('Fixed GNU stat changed after capture')
        return dict(utility_after_sha256=after,exit_code=got.returncode,raw=got.stdout,stderr=got.stderr,rows=parsed,complete_exact_order_and_decimal_strings=valid,utility_sha256=tool['sha256'])
    report['GNU_before']=gnu_batch()
    save()
    for path in paths:
        row={'path':path};report['rows'].append(row)
        try:
            q=Path(path).resolve();row['realpath']=str(q)
            row['python_before']=statrow(q);row['C_before']=cstat(q)
            row['sha256']=sha(q)
            row['python_after']=statrow(q);row['C_after']=cstat(q)
            row['python_before_after_equal']=row['python_before']==row['python_after']
        except Exception as e:row['error']=type(e).__name__+': '+str(e)
        save()
    report['GNU_after']=gnu_batch()
    report['status']='OBSERVATIONS_COMPLETE';report['helper_after_sha256']=sha(x.helper);save()

if __name__=='__main__':main()

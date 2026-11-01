#!/usr/bin/env python3
"""Replay fixed failed-pilot/repair science only; no producers/jobs/network."""
import argparse,hashlib,importlib.util,json,os,shutil,subprocess,sys,tarfile,textwrap,tomllib
from pathlib import Path,PurePosixPath
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from evidence_io import extract_verified
HERE=Path(__file__).resolve().parent
AUDITOR_SHA='e2e69b19c29d7d0aace5fb16780bc35c9f91b410563a49493be7dce4c58dd893'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return tomllib.loads(Path(p).read_text()) if Path(p).suffix=='.toml' else json.loads(Path(p).read_text())
def source_releases(out):
    index=read(HERE/'sources/index.json');releases={}
    for row in index['releases']:
        freeze=read(HERE/row['freeze']);archive=HERE/row['archive']
        assert sha(HERE/row['freeze'])==row['freeze_sha256'] and sha(archive)==row['archive_sha256'] and archive.stat().st_size==row['archive_bytes']
        root=out/row['release'];root.mkdir();seen=set();prefix=freeze.get('archive_prefix',row['release']+'/')
        with tarfile.open(archive,'r:gz') as tar:
            for member in tar:
                name=member.name.removeprefix(prefix);parts=PurePosixPath(name)
                assert member.isfile() and not parts.is_absolute() and '..' not in parts.parts and '\\' not in name and ':' not in name and name not in seen and name in freeze['files']
                data=tar.extractfile(member).read();assert hashlib.sha256(data).hexdigest()==freeze['files'][name];data.decode('utf-8');assert b'\0' not in data
                target=root/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data);seen.add(name)
        assert seen==set(freeze['files']);shutil.copyfile(HERE/row['freeze'],root/'source_freeze.json');releases[row['release']]=root
    return releases
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True);a=parser.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    assert np.__version__=='2.4.6'
    (out/'source').mkdir()
    releases=source_releases(out/'source')
    original=extract_verified(HERE/'failed_pilot/data.tar.gz',read(HERE/'failed_pilot/data_manifest.json'),out/'failed')
    repair=extract_verified(HERE/'repaired_pilot/data.tar.gz',read(HERE/'repaired_pilot/data_manifest.json'),out/'repaired')
    source=releases['pyscf_two_orbital_v1'];analysis=releases['pyscf_pilot_analysis_v2'];fresh=out/'fresh_replay';fresh.mkdir()
    cmd=[sys.executable,'-B',str(analysis/'strict_numpy_bool_json.py'),'--original-source',str(source),'--pilot-dir',str(original),'--output',str(fresh/'comparison.json'),'--receipt',str(fresh/'repair_receipt.json')]
    run=subprocess.run(cmd,capture_output=True,text=True,timeout=120,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    (fresh/'stdout.txt').write_text(run.stdout);(fresh/'stderr.txt').write_text(run.stderr)
    assert run.returncode==0,'Original frozen comparator/strict serializer rejects retained pilot'
    # Preserve the independent auditor itself unchanged. Administrative
    # transfer/remote-inode checks are not rerun publicly. Its exact original
    # scientific/resource block and original need/close functions run below.
    apath=HERE/'independent/audit_recovered_pilot.py';assert sha(apath)==AUDITOR_SHA;text=apath.read_text()
    defs=text[text.index(' def need('):text.index(' before=load(')]
    maths=text[text.index(' original=load('):text.index(' inputs={')]
    function='def original_scientific_fragment():\n checks=[]\n'+defs+maths+' return dict(checks=checks,summaries=summaries,project=project,projectpeak=projectpeak,replay_max=replay_max)\n'
    logical=out/'scientific_roles';logical.mkdir()
    # These are byte-exact original scientific files, including the actual
    # failed workflow receipt. No original execution status is synthesized.
    for dest,src in [('original_results',original),('results',repair),('independent_replay',fresh),('original_source',source)]:
        shutil.copytree(src,logical/dest)
    spec=importlib.util.spec_from_file_location('original_pilot_math',apath);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    scope=dict(np=np,math=module.math,re=module.re,load=module.load,coo=module.coo,p=logical)
    exec(compile(function,str(apath)+'#unchanged_scientific_resource_fragment','exec'),scope)
    result=scope['original_scientific_fragment']();assert all(row['passed'] for row in result['checks'])
    record=read(HERE/'repaired_pilot/independent_audit_summary.json');old=read(original/'execution.json')
    assert old['status']=='UNRESOLVED' and [row['exit_code'] for row in old['processes']]==[0,0,1]
    summary=dict(schema='PYSCF_PUBLIC_PILOT_RECOVERY_REPLAY_V1',audit_status='PASS',catalogue_status={'P23':'NOT_RUN','P24':'NOT_RUN'},scientific_completions=0,producer_reruns=False,scientific_solver_executed=False,original_execution_status='UNRESOLVED',original_process_exits=[0,0,1],original_independent_audit_sha256=record['original_complete_audit_sha256'],original_mathematical_auditor_sha256=AUDITOR_SHA,unchanged_scientific_fragment_sha256=hashlib.sha256(maths.encode()).hexdigest(),check_count=len(result['checks']),checks=result['checks'],projected_total_seconds=sum(result['project'].values()),maximum_projected_peak_mib=max(result['projectpeak'].values()),original_resource_policy_unchanged=True,public_transport_or_root_admission_checks_replayed=False,production_admitted=False)
    (out/'replay_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n');print(json.dumps({k:summary[k] for k in ('audit_status','catalogue_status','scientific_completions','check_count','projected_total_seconds','maximum_projected_peak_mib')},indent=2))
if __name__=='__main__':main()

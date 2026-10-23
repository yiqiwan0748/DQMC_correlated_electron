#!/usr/bin/env python3
"""Replay retained full-production operators/eigensystems; no producers or jobs."""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from evidence_io import extract_verified
from replay_pilot import source_releases
HERE=Path(__file__).resolve().parent
AUDITOR_SHA='4b33bd92bbef0a719f47a11d65a0db0a297b0374096efde97ccb62f43d54b0b7'
ORIGINAL_AUDIT_SHA='0ca893c030f775dd6d912f4a1b89e8a4e3cc3339e7180b23ccb0dba73e068ab6'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
 assert np.__version__=='2.4.6'
 assert sha(HERE/'independent/audit_pyscf_production.py')==AUDITOR_SHA and sha(HERE/'production/original_independent_audit.json')==ORIGINAL_AUDIT_SHA
 (out/'releases').mkdir();releases=source_releases(out/'releases')
 raw=extract_verified(HERE/'production/data.tar.gz',load(HERE/'production/data_manifest.json'),out/'raw')
 # All original scientific execution bytes are retained directly. The source
 # archive is independently verified before any retained analysis is executed.
 source=releases['pyscf_two_orbital_production_v2'];source.rename(raw/'source')
 replay=raw/'independent_replay';replay.mkdir()
 cmd=[sys.executable,'-B',str(raw/'source/compare.py'),'--mode','production','--target',str(raw/'results/target.toml'),'--reference',str(raw/'results/reference.json'),'--output',str(replay/'comparison.json')]
 result=subprocess.run(cmd,capture_output=True,text=True,timeout=300,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
 (out/'comparator.stdout.txt').write_text(result.stdout);(out/'comparator.stderr.txt').write_text(result.stderr)
 assert result.returncode==0,'Frozen full comparator did not pass; no catalogue promotion'
 files={str(q.relative_to(raw)):{'sha256':sha(q),'bytes':q.stat().st_size} for role in ['source','results','independent_replay'] for q in sorted((raw/role).rglob('*')) if q.is_file()}
 manifest=out/'public_scientific_input_manifest.json';write(manifest,dict(schema='PUBLIC_PYSCF_SCIENTIFIC_INPUT_MANIFEST_V1',files=files))
 cmd=[sys.executable,'-B',str(HERE/'independent/audit_pyscf_production.py'),'--receipt-root',str(raw),'--output',str(out/'public_independent_audit.json'),'--public-input-manifest',str(manifest),'--expected-public-manifest-sha256',sha(manifest)]
 run=subprocess.run(cmd,capture_output=True,text=True,timeout=300,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
 (out/'auditor.stdout.txt').write_text(run.stdout);(out/'auditor.stderr.txt').write_text(run.stderr)
 assert run.returncode==0,'Independent full scientific audit did not pass; no catalogue promotion'
 audit=load(out/'public_independent_audit.json');recorded=load(HERE/'production/original_independent_audit.json');comp=load(raw/'results/comparison.json')
 assert audit['audit_status']==audit['integrity']=='PASS' and audit['catalogue_status']==recorded['catalogue_status']=={'P23':'PASS','P24':'PASS'} and audit['complete_case_count']==7 and audit['complete_manifold_count']==56 and audit['pilot_completion_credit']==0
 # Public and original private audits differ in transport/hash-check counts;
 # their complete independently reconstructed scientific case metrics agree.
 assert audit['case_metrics']==recorded['case_metrics'] and audit['frozen_comparator_checks']==recorded['frozen_comparator_checks']
 summary=dict(schema='PYSCF_FULL_PUBLIC_REPLAY_V1',audit_status='PASS',integrity='PASS',catalogue_status=audit['catalogue_status'],scientific_completions=2,pilot_completion_credit=0,complete_case_count=7,complete_manifold_count=56,public_independent_check_count=audit['check_count'],original_independent_check_count=recorded['check_count'],frozen_comparator_checks=len(comp['checks']),all_frozen_gates_pass=all(c['passed'] for c in comp['checks']),all_case_metrics_exact= audit['case_metrics']==recorded['case_metrics'],case_metrics=audit['case_metrics'],original_audit_sha256=ORIGINAL_AUDIT_SHA,auditor_sha256=AUDITOR_SHA,source_archive_sha256=audit['source_archive_sha256'],source_freeze_sha256=audit['source_freeze_sha256'],execution_original_sha256=audit['execution_original_sha256'],retained_comparison_sha256=sha(raw/'results/comparison.json'),retained_target_sha256=sha(raw/'results/target.toml'),retained_reference_sha256=sha(raw/'results/reference.json'),data_archive_sha256=sha(HERE/'production/data.tar.gz'),fresh_comparator_sha256=sha(replay/'comparison.json'),fresh_public_independent_audit_sha256=sha(out/'public_independent_audit.json'),public_manifest_sha256=sha(manifest),maximum_comparator_replay_float_difference=audit['local_vs_ISAAC_replay_max_float_difference'],scientific_producer_executed=False,remote_transport_checks_replayed=False,original_scientific_execution_retained=True)
 write(out/'replay_summary.json',summary);print(json.dumps({k:summary[k] for k in ['audit_status','catalogue_status','scientific_completions','complete_case_count','complete_manifold_count','public_independent_check_count','frozen_comparator_checks']},indent=2))
if __name__=='__main__':main()

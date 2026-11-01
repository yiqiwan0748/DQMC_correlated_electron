#!/usr/bin/env python3
"""Replay published P21 raw-bin analysis; never run ALF, Julia or Monte Carlo."""
import argparse,gzip,hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
from evidence_io import extract_verified,sha

HERE=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((HERE/'production_data_manifest.json').read_text())
    raw=extract_verified(HERE/manifest['archive'],manifest,out/'raw')
    report=json.loads((HERE/'audit_summary.json').read_text());stored=HERE/'independent_statistics.json.gz'
    if sha(stored)!=report['independent_statistics_gzip_sha256']:raise ValueError('Reference statistics changed')
    expected=json.loads(gzip.decompress(stored.read_bytes()))
    fresh=out/'independent_audit.json'
    argv=[sys.executable,'-B',str(HERE/'independent/audit_p21_production.py'),'--job-dir',str(raw),'--release-dir',str(raw/'source'),'--freeze',str(HERE/'production_freeze.json'),'--source-review',str(HERE/'independent_source_review.json'),'--native-source-dir',str(raw/'native_source'),'--output',str(fresh)]
    run=subprocess.run(argv,capture_output=True,text=True,timeout=420)
    (out/'audit.stdout.txt').write_text(run.stdout);(out/'audit.stderr.txt').write_text(run.stderr)
    if run.returncode:raise ValueError('Independent arithmetic audit rejected; preserved output in '+str(out))
    result=json.loads(fresh.read_text())
    if result['audit_status']!='PASS' or result['scientific_status']!='UNRESOLVED':raise ValueError('Scientific result or integrity differs')
    spec=importlib.util.spec_from_file_location('p21_independent_comparison',HERE/'independent/audit_p21_sampled_pilot.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    m.equal(result['independent_reconstruction'],expected,'Complete independent real/complex statistics')
    summary={'schema':'published_p21_scientific_replay_v1','audit_status':'PASS','scientific_status':'UNRESOLVED','scientific_completions':0,'all_1408_real_and_704_complex_statistics_agree':True,'scalar_checks':m.CHECKS,'maximum_numeric_difference':m.MAX_NUMERIC_DIFFERENCE,'gates_and_labels_exact':True,'interval_arithmetic_rule':'Existing64*binary64epsilon*max(1,abs(a),abs(b)); scientific thresholds unchanged','scientific_selection':'Three controller submission/queue metadata files excluded; scientific source/raw bytes and629estimator roles unchanged. Derived transfer manifest and root mappings are explicit.','physics_executed':False}
    (out/'replay_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()

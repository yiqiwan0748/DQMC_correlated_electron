#!/usr/bin/env python3
"""Replay finite-fixture X09 certification from existing source/receipts only."""
import argparse,json,subprocess,sys
from pathlib import Path
from evidence_io import extract_verified
HERE=Path(__file__).resolve().parent
DB=HERE.parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
 manifest=json.loads((HERE/'data_manifest.json').read_text());raw=extract_verified(HERE/manifest['archive'],manifest,out/'raw')
 tranche=DB/'ground_certification_precision_v2';old=DB/'ground_certification';base=DB/'reporting/reports/ground_initial_v10/ledger_100.json'
 fresh=out/'comparison.json'
 subprocess.run([sys.executable,'-B',str(tranche/'compare.py'),'--target',str(raw/'results/target.toml'),'--oracle',str(raw/'results/oracle.json'),'--manifest',str(tranche/'manifest.toml'),'--source-root',str(raw/'native_source'),'--output',str(fresh)],check=True)
 if json.loads(fresh.read_text())!=json.loads((raw/'results/comparison.json').read_text()):raise ValueError('Frozen comparator changed')
 audit=out/'audit.json'
 subprocess.run([sys.executable,'-B',str(HERE/'audit_ground_precision_recovery.py'),'--raw',str(raw),'--tranche',str(tranche),'--original-tranche',str(old),'--base-ledger',str(base),'--output',str(audit)],check=True)
 report=json.loads(audit.read_text())
 if report['status']!='PASS' or report['scientific_status']!={'X09':'PASS','X12':'PASS'} or report['new_distinct_passes']!=1:raise ValueError('Certification rejected')
 print(json.dumps({'audit_status':'PASS','scientific_status':'PASS','catalogue_id':'X09','new_distinct_passes':1,'X12_not_counted_again':True,'checks':report['check_count'],'solver_executed':False},indent=2))
if __name__=='__main__':main()

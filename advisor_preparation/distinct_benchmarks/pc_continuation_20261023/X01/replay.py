#!/usr/bin/env python3
"""Replay X01's fixed scientific text data; no sampling, jobs or network calls."""
import argparse,json,subprocess,sys
from pathlib import Path
from evidence_io import extract_verified
HERE=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
 out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
 manifest=json.loads((HERE/'data_manifest.json').read_text())
 data=extract_verified(HERE/manifest['archive'],manifest,out/'data')
 subprocess.run([sys.executable,'-B',str(HERE/'audit_public_x01.py'),'--raw',str(data/'raw'),'--source',str(data/'source'),'--decoded-only','--selection',str(HERE/'selection_manifest.json'),'--output',str(out/'audit')],check=True)
 d=json.loads((out/'audit/audit.json').read_text())
 if d['audit_status']!='PASS' or d['scientific_status']!='UNRESOLVED' or d['new_scientific_completions']!=0 or d['selected_stage']!=2:raise ValueError('Unexpected scientific replay')
 receipt={k:d[k] for k in ['schema','audit_status','scientific_status','catalogue_id','new_scientific_completions','selected_stage','reason','check_count','failure_count','auditor_sha256','independent_math_auditor_sha256','original_full_audit_sha256']}
 receipt['physics_executed']=False;receipt['source_raw_bytes_changed']=False
 (out/'replay_summary.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()

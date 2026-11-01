#!/usr/bin/env python3
"""Layout adapter for the frozen post-failure Fraction analysis; zero benchmark credit."""
import argparse,hashlib,io,json,subprocess,sys,tarfile,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(b):return hashlib.sha256(b).hexdigest()
def member(archive,name):
 with tarfile.open(fileobj=io.BytesIO(archive.read_bytes()),mode='r:gz') as source:
  selected=[m for m in source if m.name==name]
  if len(selected)!=1 or not selected[0].isfile() or selected[0].issym() or selected[0].islnk():raise ValueError('Missing/unsafe unique saved input')
  return source.extractfile(selected[0]).read()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--expected-manifest-sha256',required=True);parser.add_argument('--output-directory',type=Path,required=True);args=parser.parse_args()
 if sha((ROOT/'PUBLIC_INPUT_MANIFEST.json').read_bytes())!=args.expected_manifest_sha256:raise SystemExit('Reviewed public manifest differs')
 out=args.output_directory.resolve();out.mkdir(parents=True,exist_ok=False)
 # The metadata-only checker binds all selected members before posthoc math.
 subprocess.run([sys.executable,'-B',str(ROOT/'tools/replay_integrity.py'),'--expected-manifest-sha256',args.expected_manifest_sha256,'--output',str(out/'metadata_integrity.json')],check=True,timeout=60)
 release=ROOT/'conditioning/analysis_release_v1';freeze=json.loads((release/'source_freeze.json').read_text())
 for name,digest in freeze['files'].items():
  if sha((release/name).read_bytes())!=digest:raise ValueError('Changed separately frozen analysis source')
 with tempfile.TemporaryDirectory(prefix='x14-posthoc-input-map-') as temporary:
  staging=Path(temporary);original=staging/'original_source';original.mkdir()
  archive=ROOT/'archives/x14-original-scientific-source.tar.gz'
  # Exact original input bytes are placed under the paths required by the
  # NEW saved-context analyzer. This does not run the old renderer/guard.
  for name in ['render_quspin_reference.py','manifest.toml']:(original/name).write_bytes(member(archive,'spectral_peak_shift_v1/'+name))
  (original/'source_freeze.json').write_bytes((ROOT/'freezes/original_source_freeze.json').read_bytes())
  trace=staging/'trace.json';trace.write_bytes(member(ROOT/'archives/x14-trace-diagnostic-scientific-data.tar.gz','trace_diagnostic/results/diagnostic/diagnostic.json'))
  command=[sys.executable,'-B',str(release/'analyze_saved_context.py'),'--trace',str(trace),'--renderer',str(original/'render_quspin_reference.py'),'--output',str(out/'conditional_conditioning_replay.json')]
  with (out/'stdout.txt').open('x') as so,(out/'stderr.txt').open('x') as se:result=subprocess.run(command,stdout=so,stderr=se,check=False,timeout=60)
 if result.returncode!=0:raise ValueError('Frozen posthoc saved-context analysis failed')
 expected=(ROOT/'conditioning/conditional_conditioning_v1.json').read_bytes();actual=(out/'conditional_conditioning_replay.json').read_bytes()
 if actual!=expected:raise ValueError('Frozen posthoc result differs byte-for-byte')
 receipt={'schema':'PUBLIC_X14_POSTHOC_CONDITIONING_REPLAY_V1','status':'PASS_CONDITIONAL_DIAGNOSTIC_ONLY','public_manifest_sha256':args.expected_manifest_sha256,'result_sha256':sha(actual),'analysis_source_sha256':sha((release/'analyze_saved_context.py').read_bytes()),'exact_frozen_report_bytes_equal':True,'input_mapping':'Exact public original renderer/manifest/freeze and saved trace staged temporarily for new posthoc analyzer; no original guard or renderer executed.','original_reference_runtime_replayed':False,'root_census_recomputed':False,'original_probe_or_gate_or_padding_changed':False,'new_fixture_eigensolver_sampling':False,'catalogue_status':'NOT_RUN','scientific_completions':0,'scope':'Only the retained failed bracket under the frozen normal-binary64 conditional arithmetic model; not hardware interval certification, a full root census or a solver repair.'}
 with (out/'replay_receipt.json').open('x') as stream:json.dump(receipt,stream,indent=2,allow_nan=False);stream.write('\n')
 print(receipt['status'])
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Public metadata/member integrity only; never renders curves or reconstructs roots."""
import argparse,hashlib,io,json,math,tarfile,tomllib
from pathlib import Path,PurePosixPath
ROOT=Path(__file__).resolve().parents[1]
def sha(b):return hashlib.sha256(b).hexdigest()
def load(name):return json.loads((ROOT/name).read_text())
def safe(name):
 p=PurePosixPath(name);return not p.is_absolute() and '..' not in p.parts and '\\' not in name and ':' not in name
def unpack(body):
 rows={}
 with tarfile.open(fileobj=io.BytesIO(body),mode='r:gz') as t:
  for m in t:
   if not m.isfile() or m.issym() or m.islnk() or not safe(m.name) or m.name in rows:raise ValueError('Unsafe/duplicate public member')
   rows[m.name]=t.extractfile(m).read()
 return rows
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--expected-manifest-sha256',required=True);parser.add_argument('--output',type=Path);args=parser.parse_args();manifest_bytes=(ROOT/'PUBLIC_INPUT_MANIFEST.json').read_bytes();assert sha(manifest_bytes)==args.expected_manifest_sha256;manifest=json.loads(manifest_bytes);checks=[]
 def need(ok,name):
  checks.append({'name':name,'passed':bool(ok)})
  if not ok:raise ValueError(name)
 need(manifest['schema']=='PUBLIC_X14_FAILURE_PACKET_INPUT_MANIFEST_V1','Named public manifest schema')
 for name,row in manifest['files'].items():
  need(safe(name),'Safe public file '+name);p=ROOT/name;need(p.is_file() and not p.is_symlink() and p.stat().st_size==row['bytes'] and sha(p.read_bytes())==row['sha256'],'Exact public byte '+name)
 index=load('archive_members.json');archives={};nested_members=0
 def member_tree(rows,entries,label):
  nonlocal nested_members
  need(set(rows)=={r['member'] for r in entries} and len(entries)==len(rows),'Exact complete selected member set '+label)
  for entry in entries:
   body=rows[entry['member']];need(sha(body)==entry['sha256'] and len(body)==entry['bytes'],'Exact selected member '+label+'::'+entry['member'])
   if entry.get('nested_scientific_source_archive'):
    nested_members+=len(entry['nested_members']);member_tree(unpack(body),entry['nested_members'],label+'::'+entry['member'])
   else:body.decode('utf8')
 for name,row in index['archives'].items():
  body=(ROOT/name).read_bytes();need(sha(body)==row['sha256'] and len(body)==row['bytes'],'Exact public archive '+name);archives[name]=unpack(body);member_tree(archives[name],row['members'],name)
 selection=load('selection.json');need(selection['catalogue_id']=='X14' and selection['catalogue_status']=='NOT_RUN' and selection['scientific_completions']==0 and not selection['full_original_diagnostic_guard_replay_supported'] and not selection['original_full_diagnostic_archive_published'] and selection['conditioning_diagnostic_included'] is True,'Honest excluded failure/selection scope')
 original=archives['archives/x14-original-scientific-source.tar.gz'];freeze=load('freezes/original_source_freeze.json');prefix=freeze['prefix'].rstrip('/')+'/';need(sha((ROOT/'archives/x14-original-scientific-source.tar.gz').read_bytes())==selection['original_source_archive_sha256'] and len(original)==len(freeze['files'])==39,'Exact original39-member scientific source')
 for name,digest in freeze['files'].items():need(sha(original[prefix+name])==digest,'Exact original scientific freeze member '+name)
 selected=archives['archives/x14-diagnostic-source-scientific-selection.tar.gz'];diagfreeze=load('freezes/original_full_diagnostic_freeze.json');dp=diagfreeze['prefix'].rstrip('/')+'/';expected=set(diagfreeze['files'])-{'prior/original_reference.stderr.txt'};need(len(expected)==len(selected)==12 and set(selected)=={dp+name for name in expected},'Exact selected12of13 diagnostic source roles')
 for name in expected:need(sha(selected[dp+name])==diagfreeze['files'][name],'Exact selected diagnostic source member '+name)
 projected=load('derived/exception_projections.json');need(len(projected['projections'])==4 and projected['traceback_bodies_published'] is False and all(q['derived_not_original'] is True and len(q['original_sha256'])==64 for q in projected['projections']),'Four clearly derived traceback projections')
 failure=archives['archives/x14-original-failure-scientific-data.tar.gz'];trace=archives['archives/x14-trace-diagnostic-scientific-data.tar.gz'];get=lambda rows,name:json.loads(rows[name])
 need('original_failure/results/reference.json' not in failure and 'original_failure/results/comparison.json' not in failure,'Original absent reference/comparison remain absent')
 e=get(failure,'original_failure/results/execution.json');need(e['schema']=='X14_READONLY_ANALYSIS_EXECUTION_V1' and e['job_id']=='6466312' and e['status']=='UNRESOLVED' and [(p['name'],p['exit_code']) for p in e['processes']]==[('native',0),('reference',1)] and e['scientific_completions']==0,'Original native success/reference failure/no comparator faithfully retained')
 n=get(failure,'original_failure/results/native.json');cfg=tomllib.loads(original[prefix+'manifest.toml'].decode());old=tomllib.loads(original[prefix+'prior/target.toml'].decode());selected_old=next(row for row in old['P09']['rows'] if (row['spin'],row['momentum_index'])==(cfg['pilot_spin'],cfg['pilot_momentum']));specs=[]
 for branch in cfg['branches']:
  groups=selected_old['addition']+selected_old['removal'] if branch=='combined' else selected_old[branch];specs.append((f"{cfg['pilot_spin']}_k{cfg['pilot_momentum']}_{branch}",'physical',branch,groups))
 for name in ['single','equal','unequal','rescaled']:
  poles=cfg['controls'][name+'_pole' if name=='single' else name+'_poles'];weights=cfg['controls'][name+'_weight' if name=='single' else name+'_weights'];specs.append(('control_'+name,'control','control',[dict(omega=p,weight=w,multiplicity=1,span=0.) for p,w in zip(poles,weights)]))
 expected_rows=[(ident,role,branch,groups,eta) for ident,role,branch,groups in specs for eta in cfg['etas']]
 need(n['schema']=='X14_RENDERED_RAW_V1' and n['backend']=='actual_native_prior_ED_new_numpy_renderer' and n['mode']=='pilot' and n['scientific_completions']==0 and n['producer_or_sampler_executed'] is False and len(n['rows'])==35 and len(expected_rows)==35,'Retained native35row pilot metadata only')
 for row,(ident,role,branch,groups,eta) in zip(n['rows'],expected_rows):
  need((row['id'],row['role'],row['branch'],row['eta'])==(ident,role,branch,eta) and row['groups']==groups and row['poles']==[g['omega'] for g in groups] and row['weights']==[g['weight'] for g in groups],'Exact retained prior pole/weight metadata '+ident+'/'+str(eta));need(all(math.isfinite(v) and v>=0 for v in row['weights']) and math.isfinite(row['full_weight']) and row['full_weight']>0,'Finite nonnegative retained weights metadata '+ident+'/'+str(eta))
 d=get(trace,'trace_diagnostic/results/diagnostic/diagnostic.json');de=get(trace,'trace_diagnostic/results/execution.json');pins=json.loads(original[prefix+'reference_library_pins.json']);need(de['job_id']=='6466728' and de['status']=='CAPTURED_ORIGINAL_FAILURE' and de['process']['exit_code']==1 and de['scientific_completions']==0 and d['original_failed_job']=='6466312' and d['diagnostic_status']=='CAPTURED_ORIGINAL_FAILURE' and d['original_source_unchanged'] is True and d['history_truncated'] is False and d['reference_library_sha256_before']==d['reference_library_sha256_after']==pins['files'],'Original unchanged trace capture metadata and30runtime hashes')
 need(all(d[k] is False for k in ['extra_scientific_probes','native_renderer_executed','new_eigensolver_or_sampler_executed','original_functions_replaced','original_padding_or_gates_changed','automatic_remedy_or_catalogue_promotion']) and d['catalogue_status']=='NOT_RUN' and d['scientific_completions']==0 and len(d['completed_original_reference_rows'])==0 and len(d['captured_failures'])==1,'Trace-only instrumentation/no missing output invented/no credit')
 conditioning=load('conditioning/conditional_conditioning_v1.json');peer=load('conditioning/independent_replay_review.json');cf=load('conditioning/analysis_release_v1/source_freeze.json')
 for name,digest in cf['files'].items():need(sha((ROOT/'conditioning/analysis_release_v1'/name).read_bytes())==digest,'Frozen separate conditioning source '+name)
 need(sha((ROOT/'conditioning/conditional_conditioning_v1.json').read_bytes())=='49852efe9f022232f067ee9abe36eab85801fdb6e8d69b232de97a151ad16a8a' and conditioning['trace_sha256']==sha(trace['trace_diagnostic/results/diagnostic/diagnostic.json']) and conditioning['analysis_source_sha256']==cf['files']['analyze_saved_context.py'] and peer['status']=='PASS_CONDITIONAL_DIAGNOSTIC_ONLY' and peer['original_diagnostic_sha256']==peer['fresh_replay_sha256']==sha((ROOT/'conditioning/conditional_conditioning_v1.json').read_bytes()) and conditioning['scientific_completions']==0 and conditioning['catalogue_status']=='NOT_RUN','Separate approved conditional diagnostic identity only; no math performed by integrity replay')
 result={'schema':'PUBLIC_X14_FAILURE_METADATA_REPLAY_V1','status':'PASS','check_count':len(checks),'checks':checks,'public_manifest_sha256':sha(manifest_bytes),'selected_archive_count':len(archives),'nested_scientific_source_members_checked':nested_members,'native_metadata_rows':35,'retained_physical_rows':15,'retained_control_rows':20,'reference_output_reconstructed':False,'comparator_replayed':False,'root_census_recomputed_or_verified':False,'curve_or_peak_numerical_agreement_tested':False,'conditioning_conclusion_tested':False,'catalogue_status':'NOT_RUN','scientific_completions':0,'scope':'Only immutable selected bytes/source/member/metadata integrity. This does not turn the failed/incomplete X14 pilot into a scientific result.'}
 if args.output:
  with args.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
 print(json.dumps({k:result[k] for k in ['status','check_count','catalogue_status','scientific_completions']}))
if __name__=='__main__':main()

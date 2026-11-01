"""Public saved-byte/field consistency check; no scientific calculation or network."""
import argparse,csv,hashlib,json,math,re
from pathlib import Path,PurePosixPath
P=Path(__file__).resolve().parent
STAGE='FINAL_SCIENTIFIC_RECORD'
RECEIPT='Public commit/readback evidence and the active ledger are maintained separately; this artifact does not itself claim upload completion.'
SCIENCE='PASS_ROOT_RECONCILED_FINITE_SCOPE'
def need(v,message):
 if not v:raise ValueError(message)
def pin(p):
 b=p.read_bytes();return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def load(p):
 def pairs(rows):
  d={}
  for k,v in rows:need(k not in d,'Duplicate JSON key');d[k]=v
  return d
 return json.loads(p.read_bytes(),object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in()).throw(ValueError('Nonfinite JSON token '+x)))
def same(a,b):
 if type(a)is not type(b):return False
 if type(a)is float:return a.hex()==b.hex()
 if type(a)is dict:return set(a)==set(b)and all(same(a[k],b[k])for k in a)
 if type(a)is list:return len(a)==len(b)and all(same(x,y)for x,y in zip(a,b))
 return a==b

def finite(x):
 if type(x)is float:need(math.isfinite(x),'Nonfinite saved float')
 if type(x)is dict:
  for y in x.values():finite(y)
 if type(x)is list:
  for y in x:finite(y)
def key(x):return x['backend'],x['arm'],x['time_numerator']
def canonical(name):
 p=PurePosixPath(name);need(not p.is_absolute()and str(p)==name and '..'not in p.parts and '.'not in p.parts and '\\'not in name,'Safe relative member required');return p
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--repository-root');args=parser.parse_args()
m=load(P/'MANIFEST.json');need(m['schema']=='P13_FINAL_PUBLIC_SCIENTIFIC_PACKET_MANIFEST_V3'and m['record_creation_stage']==STAGE and m['publication_receipt_scope']==RECEIPT,'Final scientific record manifest identity differs');files=m['files'];actual={p.relative_to(P).as_posix()for p in P.rglob('*')if p.is_file()};need(actual==set(files)|{'MANIFEST.json'},'Closed final file set differs')
for name,want in files.items():
 canonical(name);q=P/name;need(q.is_file()and not q.is_symlink()and q.resolve().is_relative_to(P)and pin(q)==want,'Public file bytes differ: '+name)
for p in P.rglob('*'):need(not p.is_symlink(),'Public packet contains symlink')
r=load(P/'comparison_scientific_projection.json');w=load(P/'WORKFLOW_STATUS.json');scope=load(P/'AUDIT_SCOPE.json');refs=load(P/'SOURCE_REFERENCES.json');audit=load(P/'independent_audit_scientific_summary.json');raw=load(P/'independent_raw_state_checks.json');src=load(P/'AUDIT_SOURCE_REFERENCES.json');final=load(P/'FINAL_RECONCILIATION.json')
for x in (r,w,scope,refs,audit,raw,src,final):finite(x)
need(all(x['record_creation_stage']==STAGE and x['publication_receipt_scope']==RECEIPT and x['scientific_status']==SCIENCE and x['catalogue_id']=='P13'and x['distinct_scientific_questions']==1 and x['proposed_publication_credit']==1 and x['active_ledger_updated_by_packet']is False for x in(w,scope,refs,audit,raw,src,final)), 'Final scientific/publication/one-question scope differs')
need(w['root_final_reconciliation_complete']is True and w['final_scientific_PASS_decision']is True and final['root_final_reconciliation_complete']is True and final['final_recorded_association_checks']=={'passed':1890,'failed':0},'Final reconciliation decision differs')
need('ledger_counts'not in w and 'catalogue_status'not in w and 'scientific_credit'not in w,'Final artifact must keep active ledger separate')
need(r['catalogue_status']=='NOT_RUN'and r['scientific_credit']==0 and w['original_workflow']['catalogue_status']=='NOT_RUN'and w['original_workflow']['scientific_credit']==0,'Original historical zero-credit/pending fields must remain unchanged')
need(w['independent_full_132_case_artifact_and_raw_purity_audit']['root_final_reconciliation_complete']is True and w['independent_full_132_case_artifact_and_raw_purity_audit']['physical_case_outputs_in_this_final_packet']==132 and w['independent_full_132_case_artifact_and_raw_purity_audit']['all133_original_durable_files_byte_verified']is True,'Complete full audit association differs')
need(final['actual_full_receipt_byte_fingerprint']=={'bytes':1285752,'sha256':'5dee207cb1a047092506dd4aa69fa6a8a08f052ec5939d3d73d7dafab52e3b4b'}and final['earlier_failures_preserved']is True and final['no_new_scientific_calculation_by_this_record']is True,'Accepted original evidence/scope association differs')
need(r['status']=='PASS_FULL_PENDING_INDEPENDENT_AUDIT'and r['comparison_run_completed']is True and r['pilot_only']is False,'Original report identity differs');need('runtime_before'not in r and 'runtime_after'not in r,'Private runtime dictionary present');need(len(r['gates'])==11692 and all(g['status']=='PASS'for g in r['gates'])and r['failed_gate_names']==[],'Original ordered gate record differs');need(len(r['state_certificates'])==132 and len(r['model_certificates'])==4 and len(r['initial_state_certificates'])==2 and len(r['growth_rates'])==4,'Original certificate coverage differs');need(r['conditional_floating_certificate_not_interval_arithmetic']is True and r['arithmetic']['conditional_not_interval']is True,'Original conditional arithmetic scope differs')
coverage={(b,a,k)for b in ('native','tenpy')for a in ('coarse','primary')for k in range(33)};certs={key(x):x for x in r['state_certificates']};need(len(certs)==132 and set(certs)==coverage,'Original state keys differ')
with (P/'time_series.csv').open(newline='')as f:rows=list(csv.DictReader(f))
need(len(rows)==132,'CSV coverage differs')
for row,c in zip(rows,r['state_certificates']):
 need((row['backend'],row['arm'],int(row['time_numerator']),int(row['time_denominator']))==(c['backend'],c['arm'],c['time_numerator'],c['time_denominator']),'CSV key differs')
 for col,value in [('entropy_center',c['physical_entropy_centered_on_actual_library_value']['center']),('entropy_radius',c['physical_entropy_centered_on_actual_library_value']['radius']),('library_S2_original',c['library_S2_original'])]:need(same(float(row[col]),value),'CSV saved binary64 field differs')
 need(same(json.loads(row['bond_ranks_original']),c['bond_ranks_original']),'CSV original bond ranks differ')
for fit in r['growth_rates']:need(fit['time_numerators']==list(range(4,13))and fit['weight_denominator']==15 and fit['weight_numerators']==[-8,-6,-4,-2,0,2,4,6,8],'Fixed growth fit/window differs')
need(audit['original_full_audit_status']=='PASS_P13_A03_INDEPENDENT_ARTIFACT_AND_RAW_S2_AUDIT_PENDING_ROOT_RECONCILIATION'and audit['full_attempt_error_count']==0 and audit['physical_case_count']==132 and audit['physical_time_numerators']==list(range(33))and audit['own_runtime_before_and_after_checks_passed']is True,'Actual full audit identity/coverage differs');need(audit['original_full_receipt_byte_fingerprint']=={'bytes':1285752,'sha256':'5dee207cb1a047092506dd4aa69fa6a8a08f052ec5939d3d73d7dafab52e3b4b'},'Original full receipt fingerprint differs')
for pass_record in audit['whole_byte_passes'].values():need(pass_record=={'members':1262,'uncompressed_bytes':692120347,'all_member_hashes_match':True},'Original full member census differs')
need(audit['original_census_members']==1000 and audit['original_census_bytes']==673011348 and audit['original_source_inventory_members']==197 and audit['source_context_admission_join_completed']is True,'Original source/census association differs');need(audit['original_callback_body_associations']=={'source_callback_body_equalities':414,'native_parts':70,'reference_parts':70,'comparison_checkpoints':274,'every_intermediate_checkpoint_prefix_and_phase_checked':True},'Original chunk coverage differs')
scalar=audit['independent_exact_scalar_reconciliation'];need(scalar['states']==132 and scalar['ordered_gate_occurrences']==11692 and len(scalar['growth_fits'])==4 and scalar['original_conditional_longdouble_certificates_recomputed']is False and audit['original_conditional_certificates_recomputed']is False,'Scalar audit scope differs')
need(raw['case_count']==len(raw['cases'])==len(raw['original_durable_case_byte_receipts'])==132 and raw['time_numerators']==list(range(33))and raw['time_denominator']==8,'Actual raw case coverage differs');need(set(map(key,raw['cases']))==coverage and len(set(map(key,raw['cases'])))==132,'Actual raw case keys differ');need({key(x)for x in raw['original_durable_case_byte_receipts']}==coverage,'Original raw-case fingerprint keys differ')
for x in raw['original_durable_case_byte_receipts']:need(type(x['bytes'])is int and x['bytes']>0 and re.fullmatch('[0-9a-f]{64}',x['sha256'])is not None,'Actual raw-case fingerprint schema differs')
case_byte_receipts={key(x):x for x in raw['original_durable_case_byte_receipts']}
for c in raw['cases']:
 original_case_bytes=(json.dumps(c,sort_keys=True,indent=2,allow_nan=False)+'\n').encode();need({'bytes':len(original_case_bytes),'sha256':hashlib.sha256(original_case_bytes).hexdigest()}=={k:case_byte_receipts[key(c)][k]for k in('bytes','sha256')},'Original durable numerical-case byte fingerprint differs')
 need(set(c)=={'backend','arm','time_numerator','time_denominator','independent_raw_export_metric','source_bound_original_tensor_join'}and c['time_denominator']==8,'Raw case field schema differs');q=c['independent_raw_export_metric'];j=c['source_bound_original_tensor_join'];need(q['status']=='PASS_ADDITIONAL_RAW_EXPORT_PHYSICAL_CHECK_ONLY_PENDING_WHOLE_AUDIT'and q['method']=='EXACT_INTEGER_QUANTIZED_FULL4096_GRAM_PURITY_AND_DIRECTED_FIXED_POINT_LOG_V1'and q['quantization_bits']==96 and q['log_bits']==128 and q['log_series_terms']==64,'Recorded raw metric method differs');need(q['original_MPS_tensor_reconstruction_and_evolution_bounds_recomputed']is False and q['normalized_raw_export_only_not_MPS_or_evolution_certificate']is True and q['no_charge_projection_or_observed_phase_alignment']is True and q['no_libm_log_sqrt_or_floating_Gram_arithmetic']is True,'Raw metric scope differs');need(j['status']=='PASS_CONDITIONAL_ORIGINAL_TENSOR_IDENTITY_JOIN_ONLY_PENDING_WHOLE_AUDIT'and j['state_key_from_original_certificate']==list(key(c))and j['original_MPS_tensors_or_eigensystem_evolution_recomputed']is False and j['original_ball_radii_independently_recomputed']is False and j['source_bound_original_tensor_export_identity_error_recomputed']is False,'Conditional original tensor join differs');need(q['scientific_credit']==j['scientific_credit']==0,'Raw subroutine credit differs')
need(len(refs['files'])==8 and len(src['files'])==4,'Scientific source count differs')
for x in src['files']:need(pin(P/x['path'])=={k:x[k]for k in('bytes','sha256')},'Pure source bytes differ')
if args.repository_root:
 root=Path(args.repository_root).resolve()
 for x in refs['files']:
  q=canonical(x['repository_path']);need(str(q).startswith('advisor_preparation/distinct_benchmarks/runtime_validation_20261028/p13_quench_pilot_rayleigh_recomparison/source/'),'Referenced source prefix differs');need(pin(root/str(q))=={k:x[k]for k in('bytes','sha256')},'Existing referenced source bytes differ')
print(json.dumps({'status':'PASS_PUBLIC_SAVED_BYTE_FIELD_AND_FINAL_FINITE_SCOPE_CONSISTENCY','files':len(files)+1,'original_ordered_gate_occurrences':11692,'original_state_records':132,'actual_additional_raw_cases':132,'scientific_status':SCIENCE,'proposed_publication_credit':1,'active_ledger_updated_by_verifier':False,'physical_computations':0,'referenced_source_files_checked':8 if args.repository_root else 0},sort_keys=True))

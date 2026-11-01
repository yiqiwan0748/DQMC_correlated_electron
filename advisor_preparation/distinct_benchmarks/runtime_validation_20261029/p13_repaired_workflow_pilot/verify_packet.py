# Public byte/status/CSV/source consistency only; no physics or network.
import argparse,csv,gzip,hashlib,io,json,math
from pathlib import Path,PurePosixPath
P=Path(__file__).resolve().parent
RAW_PIN={'bytes':1338924,'sha256':'e9c955a994947d6ccf267f120987ed6066c6ad747245a7b6d59ba7995d432008'}
GZIP_PIN={'bytes':20982,'sha256':'c388cd0da7231be58a37121024c86147d91d3685ada2844d41bcd5ee2a7b61bb'}
def pin(b):return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
def decode(b):
 def pairs(rows):
  d={}
  for k,v in rows:
   if k in d:raise ValueError('Duplicate JSON key')
   d[k]=v
  return d
 return json.loads(b,object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in()).throw(ValueError(x)))
def finite(x):
 if type(x)is float and not math.isfinite(x):raise ValueError('Nonfinite number')
 if type(x)is dict:
  for v in x.values():finite(v)
 if type(x)is list:
  for v in x:finite(v)
def same(a,b):
 if type(a)is not type(b):return False
 if type(a)is float:return a.hex()==b.hex()
 if type(a)is list:return len(a)==len(b)and all(same(x,y)for x,y in zip(a,b))
 if type(a)is dict:return set(a)==set(b)and all(same(a[k],b[k])for k in a)
 return a==b
def main():
 a=argparse.ArgumentParser();a.add_argument('--repository-root');args=a.parse_args();m=decode((P/'MANIFEST.json').read_bytes());expected=m['files'];assert set(q.name for q in P.iterdir())==set(expected)|{'MANIFEST.json'}
 assert len(expected)==7 and all((P/n).is_file()and not(P/n).is_symlink()and pin((P/n).read_bytes())==v for n,v in expected.items())
 compressed=(P/'comparison_scientific_projection.json.gz').read_bytes();assert pin(compressed)==GZIP_PIN
 with gzip.GzipFile(fileobj=io.BytesIO(compressed))as stream:raw=stream.read(RAW_PIN['bytes']+1);assert stream.read(1)==b''
 assert pin(raw)==RAW_PIN==m['projection_uncompressed'];r=decode(raw);w=decode((P/'WORKFLOW_STATUS.json').read_bytes());refs=decode((P/'SOURCE_REFERENCES.json').read_bytes());finite(r);finite(w)
 assert r['status']==w['original_comparison_report_status']=='PASS_PILOT_ZERO_CREDIT_PENDING_INDEPENDENT_AUDIT'and r['comparison_run_completed']is True and r['pilot_only']is True
 assert w['original_workflow_status']==r['status']and w['original_scheduler_state']=='COMPLETED'and w['original_scheduler_exit']=='0:0'and w['pilot_workflow_completed_successfully']is True and w['original_output_census_completed']is True and w['workflow_errors']==[]
 assert r['catalogue_status']==w['catalogue_status']=='NOT_RUN'and r['scientific_credit']==w['scientific_credit']==w['distinct_ideas_added']==0 and w['full_curve_status']=='NOT_RUN'
 assert w['ledger_counts']=={'PASS':50,'FAIL':0,'UNRESOLVED':2,'NOT_RUN':48}
 assert len(r['gates'])==8889 and all(g['status']=='PASS'for g in r['gates'])and r['failed_gate_names']==[]and len(r['state_certificates'])==20 and len(r['model_certificates'])==4 and len(r['initial_state_certificates'])==2 and r['growth_rates']==[]
 assert r['growth_rate_status']=='NOT_RUN_PILOT_FIVE_TIMES_DO_NOT_COVER_FIXED_NINE_POINT_WINDOW'and r['nontrivial_full_curve_status']=='NOT_RUN_PILOT_ONLY'
 assert 'runtime_before'not in r and 'runtime_after'not in r and r['conditional_floating_certificate_not_interval_arithmetic']is True and r['arithmetic']['conditional_not_interval']is True
 assert w['saved_association_audit']['checks']==9321 and w['saved_association_audit']['failed_checks']==0 and w['independent_saved_receipt_review']['focused_checks']==152 and w['independent_saved_receipt_review']['failed_checks']==0
 with(P/'time_series.csv').open(newline='')as f:rows=list(csv.DictReader(f))
 assert len(rows)==20
 for row,s in zip(rows,r['state_certificates']):
  b=s['physical_entropy_centered_on_actual_library_value'];assert row['backend']==s['backend']and row['arm']==s['arm']and int(row['time_numerator'])==s['time_numerator']and int(row['time_denominator'])==s['time_denominator']
  assert same(float(row['entropy_center']),b['center'])and same(float(row['entropy_radius']),b['radius'])and same(float(row['library_S2_original']),s['library_S2_original'])and same(json.loads(row['bond_ranks_original']),s['bond_ranks_original'])
 assert len(refs['files'])==8
 for x in refs['files']:
  assert x['repository_path']=='advisor_preparation/distinct_benchmarks/runtime_validation_20261028/p13_quench_pilot_rayleigh_recomparison/source/'+x['name']and x['relative_path']=='../../runtime_validation_20261028/p13_quench_pilot_rayleigh_recomparison/source/'+x['name']and type(x['bytes'])is int and x['bytes']>0 and len(x['sha256'])==64
 if args.repository_root:
  root=Path(args.repository_root).resolve()
  for x in refs['files']:
   q=PurePosixPath(x['repository_path']);assert not q.is_absolute()and '..'not in q.parts and str(q).startswith('advisor_preparation/distinct_benchmarks/runtime_validation_20261028/p13_quench_pilot_rayleigh_recomparison/source/');p=root/str(q);assert p.resolve().is_relative_to(root)and p.is_file()and not p.is_symlink()and pin(p.read_bytes())=={k:x[k]for k in('bytes','sha256')}
 print(json.dumps({'status':'PASS_PUBLIC_BYTE_STATUS_AND_SAVED_FIELD_CONSISTENCY','public_files':8,'reported_gates':8889,'state_records':20,'full_P13':'NOT_RUN','scientific_credit':0,'physical_computations':0,'referenced_source_files_checked':8 if args.repository_root else 0}))
if __name__=='__main__':main()

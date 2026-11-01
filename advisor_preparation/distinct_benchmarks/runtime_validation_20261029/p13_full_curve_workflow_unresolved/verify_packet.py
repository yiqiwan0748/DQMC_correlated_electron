# Public byte/field consistency checker; no scientific computation or network.
import argparse,csv,hashlib,json,math
from pathlib import Path,PurePosixPath
P=Path(__file__).resolve().parent

def pin(p):
 raw=p.read_bytes();return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
def load(p):
 def pairs(rows):
  d={}
  for k,v in rows:
   if k in d:raise ValueError('Duplicate JSON key')
   d[k]=v
  return d
 return json.loads(p.read_bytes(),object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in()).throw(ValueError(x)))
def same(a,b):
 if type(a)is not type(b):return False
 if type(a)is float:return a.hex()==b.hex()
 if type(a)is dict:return set(a)==set(b)and all(same(a[k],b[k])for k in a)
 if type(a)is list:return len(a)==len(b)and all(same(x,y)for x,y in zip(a,b))
 return a==b

def finite(x):
 if type(x)is float and not math.isfinite(x):raise ValueError('Nonfinite saved number')
 if type(x)is dict:
  for v in x.values():finite(v)
 if type(x)is list:
  for v in x:finite(v)
a=argparse.ArgumentParser();a.add_argument('--repository-root');args=a.parse_args()
m=load(P/'MANIFEST.json');expected=m['files'];assert set(x.name for x in P.iterdir())==set(expected)|{'MANIFEST.json'}
assert all((P/n).is_file()and not(P/n).is_symlink()and pin(P/n)==v for n,v in expected.items())
r=load(P/'comparison_scientific_projection.json');w=load(P/'WORKFLOW_STATUS.json');refs=load(P/'SOURCE_REFERENCES.json');finite(r);finite(w)
assert r['status']==w['original_comparison_report_status']=='PASS_FULL_PENDING_INDEPENDENT_AUDIT'and r['comparison_run_completed']is True and r['pilot_only']is False
assert w['publication_status']=='NUMERICAL_COMPARISON_WITH_WORKFLOW_FAILURE'and w['original_workflow_status']=='UNRESOLVED'and w['original_scheduler_state']=='FAILED'and w['original_scheduler_exit']=='2:0'and w['full_workflow_PASS']is False
assert w['catalogue_status']==r['catalogue_status']=='NOT_RUN'and w['scientific_credit']==r['scientific_credit']==w['distinct_ideas_added']==0
assert len(w['original_full_workflow_predicates_unmet'])==3 and all(x['pass']is False for x in w['original_full_workflow_predicates_unmet'])
assert len(r['gates'])==11692 and all(g['status']=='PASS'for g in r['gates'])and r['failed_gate_names']==[]
assert 'runtime_before'not in r and 'runtime_after'not in r and len(r['state_certificates'])==132 and len(r['growth_rates'])==4 and len(r['model_certificates'])==4 and len(r['initial_state_certificates'])==2
assert r['conditional_floating_certificate_not_interval_arithmetic']is True and r['arithmetic']['conditional_not_interval']is True
with(P/'time_series.csv').open(newline='')as f:rows=list(csv.DictReader(f))
assert len(rows)==132
for row,s in zip(rows,r['state_certificates']):
 b=s['physical_entropy_centered_on_actual_library_value']
 assert row['backend']==s['backend']and row['arm']==s['arm']and int(row['time_numerator'])==s['time_numerator']and int(row['time_denominator'])==s['time_denominator']
 assert same(float(row['entropy_center']),b['center'])and same(float(row['entropy_radius']),b['radius'])and same(float(row['library_S2_original']),s['library_S2_original'])and same(json.loads(row['bond_ranks_original']),s['bond_ranks_original'])
assert len(refs['files'])==8
if args.repository_root:
 root=Path(args.repository_root).resolve()
 for x in refs['files']:
  q=PurePosixPath(x['repository_path']);assert not q.is_absolute()and '..'not in q.parts and str(q).startswith('advisor_preparation/distinct_benchmarks/runtime_validation_20261028/p13_quench_pilot_rayleigh_recomparison/source/')
  assert pin(root/str(q))=={k:x[k]for k in ('bytes','sha256')}
print(json.dumps({'status':'PASS_PUBLIC_BYTE_STATUS_AND_SAVED_FIELD_CONSISTENCY','public_files':len(expected)+1,'recorded_gates':11692,'workflow_predicates_unmet':3,'full_P13':'NOT_RUN','scientific_credit':0,'physical_computations':0,'referenced_source_files_checked':8 if args.repository_root else 0}))

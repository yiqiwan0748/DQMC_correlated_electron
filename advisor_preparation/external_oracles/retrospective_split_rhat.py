#!/usr/bin/env python3
"""Recompute split-Rhat from immutable saved bin receipts; starts no simulations."""
import csv,hashlib,json,math,statistics as st,tomllib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
RESULTS=ROOT.parent/'isaac/results'
MANIFEST=ROOT.parent/'validation_design/benchmark_manifest.toml'
CASES={c['id']:c for c in tomllib.loads(MANIFEST.read_text())['cases']}

def rhat(chains):
 n=len(chains[0]);assert n%2==0 and all(len(c)==n for c in chains)
 h=n//2;parts=[c[:h] for c in chains]+[c[h:] for c in chains]
 W=st.mean(st.variance(c) for c in parts);B=h*st.variance(st.mean(c) for c in parts)
 tiny=1e-24*max(1.,max(v*v for c in parts for v in c))
 if W<=tiny:return 1. if B<=tiny else math.inf
 return math.sqrt(max(1.,((h-1)/h*W+B/h)/W))

def scalar_column(c):
 ob=c['observable']
 return {'energy_uncentered':0,'energy_centered':1,'density_total':2,'double_occupancy_total':3,'spin_z_correlation':4 if c.get('sites')==[0,1] else 5,'onebody_up':6 if c.get('sites')==[0,1] else 7}[ob]

def pseudo(chain,c):
 if c['observable']=='S2':
  nodes=chain['nodes'];n=len(nodes[0]['numerator_bins'])
  def estimate(drop=None):
   v=0.
   for node in nodes:
    a=node['numerator_bins'];b=node['denominator_bins']
    ratio=(sum(a)-(a[drop] if drop is not None else 0))/(sum(b)-(b[drop] if drop is not None else 0))
    if ratio<=0:raise ValueError('Nonpositive entropy ratio')
    v-=math.log(ratio)
   return v
  center=estimate();return [n*center-(n-1)*estimate(b) for b in range(n)]
 k=scalar_column(c);a=[v[k] for v in chain['numerator_bins']];b=chain['denominator_bins'];n=len(a);sa=sum(a);sb=sum(b);center=sa/sb
 return [n*center-(n-1)*(sa-a[i])/(sb-b[i]) for i in range(n)]

rows=[];inputs=[]
for file in [RESULTS/name for name in ['external_atomic_mc_4096.toml','external_thermal_mc_4096.toml','mc_runtime_pilot_512.toml']]:
 try:d=tomllib.loads(file.read_text())
 except (ValueError,UnicodeError):continue
 found=False
 for c in d.get('cases',[]):
  for arm in c.get('arms',[]):
   if arm.get('method')!='DQMC' or 'raw_bins' not in arm:continue
   found=True;cfg=CASES.get(c['id'],c);reason=''
   try:
    chains=[pseudo(chain,cfg) for chain in arm['raw_bins']];value=rhat(chains)
    status='PASS_DIAGNOSTIC_ONLY' if value<=1.05 else 'REJECT_MIXING'
   except (ValueError,KeyError,ZeroDivisionError,AssertionError) as ex:
    value=None;status='UNRESOLVED_DIAGNOSTIC';reason=str(ex)
   affected=c['id'] in ['ATOM-05','ATOM-06','ATOM-07','ATOM-08'] and 'lazy-local-proposals' not in d.get('mc_protocol_revision','')
   rows.append({'receipt':str(file.relative_to(ROOT.parent)),'slurm_job_id':d.get('slurm_job_id'),'case_id':c['id'],'original_arm_status':arm.get('status'),'split_rhat':value if value is None or math.isfinite(value) else str(value),'limit':1.05,'retrospective_mixing_status':status,'known_atomic_kernel_fault_affects_fixture':affected,'eligible_for_corrected_sampler_summary':False,'reason':reason})
 if found:inputs.append({'path':str(file.relative_to(ROOT.parent)),'sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'protocol':d.get('mc_protocol_revision','ORIGINAL_NO_LAZY_LOCAL_PROPOSALS')})
report={'status':'RETROSPECTIVE_DIAGNOSTIC_ONLY_ORIGINAL_RECEIPTS_UNCHANGED','simulations_started':0,'diagnostic':'Classic split-Rhat on paired leave-one-bin-out jackknife pseudovalues, matching run_external_validation.jl; not rank-normalized/folded, not a proof of equilibration','limit':1.05,'source_inputs':inputs,'rows':rows,'counts':{s:sum(r['retrospective_mixing_status']==s for r in rows) for s in ['PASS_DIAGNOSTIC_ONLY','REJECT_MIXING','UNRESOLVED_DIAGNOSTIC']},'known_fault_rows':sum(r['known_atomic_kernel_fault_affects_fixture'] for r in rows),'reporting_rule':'No original arm, including an original PASS or retrospective diagnostic pass, is merged into a corrected-kernel validation summary.'}
(ROOT/'retrospective_split_rhat.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
with (ROOT/'retrospective_split_rhat.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps({'counts':report['counts'],'known_fault_rows':report['known_fault_rows'],'inputs':len(inputs),'rows':len(rows),'rejected':[(r['case_id'],r['split_rhat']) for r in rows if r['retrospective_mixing_status']=='REJECT_MIXING']},indent=2))

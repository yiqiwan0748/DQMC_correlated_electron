#!/usr/bin/env python3
"""Independently audit saved final external campaign bins/counts; starts no jobs."""
import hashlib,json,math,statistics as st,tarfile,tomllib
from pathlib import Path
R=Path(__file__).resolve().parent.parent
FILES=['external_atomic_corrected_4096.toml','external_thermal_corrected_4096.toml','external_atomic_corrected_16384.toml','external_thermal_corrected_16384.toml']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def effective(x):
 n=len(x);mean=st.mean(x);y=[v-mean for v in x];g0=sum(v*v for v in y)/n
 if g0<=2.220446049250313e-16*max(1.,max(v*v for v in x)):return float(n)
 tau=1.
 for lag in range(1,min(n//4,16)+1):
  corr=sum(y[i]*y[i+lag] for i in range(n-lag))/((n-lag)*g0)
  if corr<=0:break
  tau+=2*corr
 return n/tau

def rhat(chains):
 h=len(chains[0])//2;parts=[c[:h] for c in chains]+[c[h:] for c in chains];W=st.mean(st.variance(c) for c in parts);B=h*st.variance(st.mean(c) for c in parts);tiny=1e-24*max(1.,max(v*v for c in parts for v in c))
 if W<=tiny:return 1. if B<=tiny else math.inf
 return math.sqrt(max(1.,((h-1)/h*W+B/h)/W))

def chain_values(raw,c):
 entropy=c['observable']=='S2'
 if entropy:nodes=[(n['numerator_bins'],n['denominator_bins']) for n in raw['nodes']]
 else:
  k={'energy_uncentered':0,'energy_centered':1,'density_total':2,'double_occupancy_total':3,'spin_z_correlation':4 if c.get('sites')==[0,1] else 5,'onebody_up':6 if c.get('sites')==[0,1] else 7}[c['observable']]
  nodes=[([v[k] for v in raw['numerator_bins']],raw['denominator_bins'])]
 n=len(nodes[0][0]);assert n==64
 def estimate(drop=None,pair=False):
  out=0.
  for aa,bb in nodes:
   a=[aa[i]+aa[i+1] for i in range(0,n,2)] if pair else aa;b=[bb[i]+bb[i+1] for i in range(0,n,2)] if pair else bb
   ratio=(sum(a)-(a[drop] if drop is not None else 0))/(sum(b)-(b[drop] if drop is not None else 0))
   if entropy:assert ratio>0;out-=math.log(ratio)
   else:out=ratio
  return out
 center=estimate();loo=[estimate(i) for i in range(n)];pseudo=[n*center-(n-1)*v for v in loo];value=st.mean(pseudo);se=math.sqrt((n-1)/n*sum((v-st.mean(loo))**2 for v in loo))
 loo2=[estimate(i,pair=True) for i in range(n//2)];se2=math.sqrt((n//2-1)/(n//2)*sum((v-st.mean(loo2))**2 for v in loo2))
 es=[]
 for a,b in nodes:
  ratio=sum(a)/sum(b) if entropy else value;es.append(effective([x-ratio*y for x,y in zip(a,b)]))
 return value,se,se2,pseudo,min(es)

def t7critical(alpha):
 def tail(q):
  t=math.atan(q/math.sqrt(7));return .5-(10*t+7.5*math.sin(2*t)+1.5*math.sin(4*t)+math.sin(6*t)/6)/(10*math.pi)
 lo=0.;hi=100.
 for _ in range(100):
  m=(lo+hi)/2
  if tail(m)>alpha/2:lo=m
  else:hi=m
 return (lo+hi)/2

def main():
 manifest=R/'validation_design/benchmark_manifest.toml';configs={c['id']:c for c in tomllib.loads(manifest.read_text())['cases']}
 with tarfile.open(R/'validation_source/corrected-source.tar.gz') as t:hashes={m.name:hashlib.sha256(t.extractfile(m).read()).hexdigest() for m in t.getmembers() if m.isfile()}
 final={};inputs=[];history=[];critical=t7critical(.01/400)
 for name in FILES:
  path=R/'isaac/results'/name;d=tomllib.loads(path.read_text());assert 'finished_utc' in d,name
  assert d['mc_protocol_revision']=='lazy-local-proposals-v1; split-Rhat<=1.05; alpha=0.01/400'
  assert d['manifest_sha256']==sha(manifest) and d['runner_sha256']==hashes['validation_design_corrected/run_external_validation.jl']
  assert d['mc_measurements_per_seed'] in [4096,16384]
  for n,h in d['source_hashes'].items():assert h==hashes['renyi_corrected/'+n],n
  for n,h in d['all_src_sha256'].items():assert h==hashes['renyi_corrected/'+n],n
  inputs.append({'path':'isaac/results/'+name,'sha256':sha(path),'job_id':d['slurm_job_id'],'measurements_per_chain':d['mc_measurements_per_seed']})
  for case in d['cases']:
   c=configs[case['id']];a=next(a for a in case['arms'] if a['method']=='DQMC');raw=a['raw_bins'];assert [x['seed'] for x in raw]==[1009,2011,3011,4001,5003,6007,7001,8009]
   chain=[chain_values(x,c) for x in raw];vals=[x[0] for x in chain];mean=st.mean(vals);se=st.stdev(vals)/math.sqrt(8);half=critical*se;rh=rhat([x[3] for x in chain]);eff=sum(x[4] for x in chain);blocking=max(abs(x[2]-x[1])/max(x[1],1e-12) for x in chain);err=abs(mean-case['reference'])
   assert max(abs(x-y) for x,y in zip(vals,a['seed_values']))<2e-10 and abs(mean-a['value'])<2e-10 and abs(se-a['standard_error'])<2e-10
   assert abs(half-a['simultaneous_half_width'])<2e-9 and abs(rh-a['split_rhat'])<2e-8
   # Tiny-variance channels have floating-point-sensitive relative blocking ratios;
   # retain the registered diagnostic gate as well as the independently recomputed one.
   good=all(math.isfinite(x) for x in [mean,half,rh]) and rh<=1.05 and eff>=128 and blocking<=.5 and a['blocking_max_relative_change']<=.5 and a['estimated_effective_bins']>=128
   status='PASS' if good and err+half<=a['equivalence_margin'] else 'FAIL' if good and err-half>a['equivalence_margin'] else 'UNRESOLVED'
   assert status==a['status'],(case['id'],status,a['status'],blocking,a['blocking_max_relative_change'],eff)
   row={'id':case['id'],'receipt':name,'measurements_per_chain':d['mc_measurements_per_seed'],'status_recomputed':status,'value_recomputed':mean,'reference':case['reference'],'standard_error_recomputed':se,'simultaneous_half_width_recomputed':half,'error_plus_half_width':err+half,'equivalence_margin':a['equivalence_margin'],'split_rhat_recomputed':rh,'effective_bins_recomputed':eff,'blocking_relative_change_recomputed':blocking,'registered_blocking_relative_change':a['blocking_max_relative_change'],'maximum_seed_mean_reconstruction_error':max(abs(x-y) for x,y in zip(vals,a['seed_values']))}
   final[case['id']]=row;history.append(row)
 assert set(final)=={k for k in configs if k.startswith(('QT-','ATOM-'))} and len(final)==24
 assert all(r['status_recomputed']=='PASS' for r in final.values())
 deterministic=R/'isaac/results/external_deterministic_attempt01.toml';dd=tomllib.loads(deterministic.read_text());assert dd['manifest_sha256']==sha(manifest)
 with tarfile.open(R/'validation_source/baseline-source.tar.gz') as t:
  for n,h in dd['source_hashes'].items():assert hashlib.sha256(t.extractfile('renyi/'+n).read()).hexdigest()==h,n
 oraclepath=R/'isaac/results/quspin_attempt02/oracle_values.toml';oracles=tomllib.loads(oraclepath.read_text())['values'];assert dd['oracle_sha256']==sha(oraclepath)
 armmap={};fixed_count=0;max_scaled_error=0.
 for c in dd['cases']:
  for a in c['arms']:
   if a['method']=='DQMC' and c['id'] in final:continue
   assert a['status']=='PASS' and math.isfinite(a['value']) and math.isfinite(a['reference'])
   error=abs(a['value']-a['reference']);assert error<=a['tolerance'] and abs(error-a['absolute_error'])<1e-12
   if c['id'] in oracles:assert abs(a['reference']-oracles[c['id']])<1e-12
   if 'residual' in a:assert a['residual']<=1e-9
   if 'variance' in a:assert -1e-9<=a['variance']<=1e-7
   if 'norm' in a:assert abs(a['norm']-1)<=1e-9
   fixed_count+=1;max_scaled_error=max(max_scaled_error,error/a['tolerance'])
   armmap[(c['id'],a['method'])]=a
 for cid,row in final.items():armmap[(cid,'DQMC')]={'status':row['status_recomputed']}
 counts={s:0 for s in ['PASS','FAIL','UNRESOLVED']};arms={s:0 for s in counts}
 for cid,c in configs.items():
  statuses=[armmap[(cid,m)]['status'] for m in c['target_methods']]
  for s in statuses:arms[s]+=1
  counts['PASS' if all(s=='PASS' for s in statuses) else 'FAIL' if 'FAIL' in statuses else 'UNRESOLVED']+=1
 summarypath=R/'isaac/results/summary/summary.json';summary=json.loads(summarypath.read_text())
 assert counts==summary['case_counts']=={'PASS':100,'FAIL':0,'UNRESOLVED':0}
 assert all(summary['arm_counts'][s]==n for s,n in arms.items()) and arms=={'PASS':188,'FAIL':0,'UNRESOLVED':0}
 summarycases={c['id']:c for c in summary['cases']}
 for cid,row in final.items():
  arm=next(a for a in summarycases[cid]['arms'] if a['method']=='DQMC');assert arm['receipt']==row['receipt'] and arm['status']=='PASS'
 assert fixed_count==164
 report={'status':'PASS_FINAL_CORRECTED_MC_RAW_BINS_PROVENANCE_AND_CAMPAIGN_COUNTS','scientific_or_timing_jobs_started':0,'method':'Independent Python standard-library paired jackknife, split-Rhat, blocking/influence-bin ESS,Student-t7 arithmetic; source hashes checked against frozen corrected archive. No original stochastic PASS reused.','manifest_sha256':sha(manifest),'input_receipts':inputs,'corrected_stochastic_cases_audited':len(final),'stochastic_inspections_audited':len(history),'critical_value_recomputed':critical,'latest_corrected_rows':list(final.values()),'all_corrected_inspection_history':history,'deterministic_receipt_sha256':sha(deterministic),'deterministic_source_hashes_match_frozen_baseline':True,'fixed_method_arms_arithmetic_and_available_residual_variance_norm_gates_checked':fixed_count,'maximum_fixed_error_divided_by_tolerance':max_scaled_error,'QuSpin_reference_receipt_sha256':sha(oraclepath),'combined_case_counts':counts,'combined_method_arm_counts':arms,'summary_sha256':sha(summarypath),'limits':['Counts cover the initial100 finite-system comparisons, not100distinct completed research ideas.','Small matched external targets do not certify all temperatures,sizes,determinant-zero/sign-problem regimes or thermodynamic limits.','Student-t intervals retain the registered nominal independent-chain approximation; diagnostics do not prove global equilibrium.'],'analysis_script_sha256':sha(Path(__file__))}
 (R/'external_oracles/final_campaign_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k not in ['latest_corrected_rows','all_corrected_inspection_history','input_receipts']},indent=2))
if __name__=='__main__':main()

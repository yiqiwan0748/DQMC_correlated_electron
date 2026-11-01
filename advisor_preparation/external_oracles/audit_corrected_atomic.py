#!/usr/bin/env python3
"""Audit the first corrected atomic receipts; launches no sampler or job."""
import hashlib,json,math,statistics as st,tarfile,tomllib
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
INPUT=ROOT/'isaac/results/external_atomic_corrected_4096.toml'
DIAG=ROOT/'isaac/results/atomic_corrected_kernel_diagnostic_v2.toml'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def split_rhat(chains):
 h=len(chains[0])//2;parts=[c[:h] for c in chains]+[c[h:] for c in chains]
 W=st.mean(st.variance(c) for c in parts);B=h*st.variance(st.mean(c) for c in parts);tiny=1e-24*max(1.,max(v*v for c in parts for v in c))
 if W<=tiny:return 1. if B<=tiny else math.inf
 return math.sqrt(max(1.,((h-1)/h*W+B/h)/W))
def pseudo(raw,observable):
 if observable=='S2':
  nodes=raw['nodes'];n=len(nodes[0]['numerator_bins'])
  def value(drop=None):
   out=0.
   for node in nodes:
    a=node['numerator_bins'];b=node['denominator_bins'];ratio=(sum(a)-(a[drop] if drop is not None else 0))/(sum(b)-(b[drop] if drop is not None else 0));assert ratio>0
    out-=math.log(ratio)
   return out
  center=value();return [n*center-(n-1)*value(k) for k in range(n)]
 k={'energy_centered':1,'density_total':2,'double_occupancy_total':3}[observable]
 a=[v[k] for v in raw['numerator_bins']];b=raw['denominator_bins'];n=len(a);sa=sum(a);sb=sum(b);center=sa/sb
 return [n*center-(n-1)*(sa-a[i])/(sb-b[i]) for i in range(n)]
def main():
 d=tomllib.loads(INPUT.read_text());dg=tomllib.loads(DIAG.read_text());assert d['mc_protocol_revision']=='lazy-local-proposals-v1; split-Rhat<=1.05; alpha=0.01/400'
 assert dg['runner_sha256']==sha(ROOT/'isaac/diagnostics/diagnose_atomic_odd_sweep_v2.jl')=='402f3614afd71ddb87f0977fcc778d9a314b40ff4bd7828443f28597d7f04dad'
 with tarfile.open(ROOT/'validation_source/corrected-source.tar.gz') as t:
  hashes={m.name:hashlib.sha256(t.extractfile(m).read()).hexdigest() for m in t.getmembers() if m.isfile()}
 for name,digest in d['source_hashes'].items():assert digest==hashes['renyi_corrected/'+name],name
 for name,digest in d['all_src_sha256'].items():assert digest==hashes['renyi_corrected/'+name],name
 assert d['runner_sha256']==hashes['validation_design_corrected/run_external_validation.jl']
 rows=[]
 for c in d['cases']:
  a=next(a for a in c['arms'] if a['method']=='DQMC');vals=a['seed_values'];mean=st.mean(vals);se=st.stdev(vals)/math.sqrt(len(vals));half=a['critical_value']*se;error=abs(mean-c['reference']);rh=split_rhat([pseudo(r,c['observable']) for r in a['raw_bins']])
  assert abs(mean-a['value'])<1e-12 and abs(se-a['standard_error'])<1e-12 and abs(half-a['simultaneous_half_width'])<1e-11
  assert abs(rh-a['split_rhat'])<1e-10
  good=math.isfinite(mean) and math.isfinite(half) and a['blocking_max_relative_change']<=.5 and a['estimated_effective_bins']>=128 and rh<=1.05
  status='PASS' if good and error+half<=a['equivalence_margin'] else 'FAIL' if good and error-half>a['equivalence_margin'] else 'UNRESOLVED'
  assert status==a['status'],c['id']
  rows.append({'id':c['id'],'observable':c['observable'],'status_recomputed':status,'value':mean,'reference':c['reference'],'absolute_error':error,'simultaneous_half_width':half,'error_plus_half_width':error+half,'equivalence_margin':a['equivalence_margin'],'split_rhat_recomputed':rh})
 kernel=dg['kernel_checks'];keys=['maximum_G00_vs_exact_field_error','maximum_Gcurrent_vs_exact_field_error','maximum_G00_vs_fresh_stabilized_error','maximum_propagator_vs_field_log_error','maximum_reported_deltaG']
 report={'status':'PASS_FILE_PROVENANCE_AND_CORRECTED_ATOMIC_RESULT_ARITHMETIC','scientific_jobs_started_by_audit':0,'receipt_sha256':sha(INPUT),'diagnostic_receipt_sha256':sha(DIAG),'corrected_source_hashes_match_frozen_archive':True,'all_src_files_verified':len(d['all_src_sha256']),'consumer_hash_matches_frozen_archive':True,'job_id':d['slurm_job_id'],'counts':{s:sum(r['status_recomputed']==s for r in rows) for s in ['PASS','FAIL','UNRESOLVED']},'rows':rows,'actual_kernel_diagnostic_rows':len(kernel),'actual_kernel_maximum_errors':{k:max(r[k] for r in kernel) for k in keys},'all_nstab2_trajectories_match_nstab1':all(r['matches_nstab1_trajectory'] for r in kernel),'every_site_revisited_alternating_class':all(all(v>0 for v in r['measurement_alternating_visits_per_site']) for r in kernel),'minimum_alternating_visits':min(v for r in kernel for v in r['measurement_alternating_visits_per_site']),'scope':'Audits saved new ISAAC receipts and their arithmetic/provenance. Seven atomic comparisons pass at4096; energy ATOM05 remains unresolved through precision. This does not certify untested regimes or rigorous finite-MCMC confidence coverage.','analysis_script_sha256':sha(Path(__file__))}
 (ROOT/'external_oracles/corrected_atomic_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
if __name__=='__main__':main()

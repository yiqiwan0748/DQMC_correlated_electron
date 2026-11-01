#!/usr/bin/env python3
"""Recompute every replicate from raw bins, then apply the predeclared X06 gates."""
import argparse,hashlib,importlib.metadata,json,math,statistics as st,tomllib
from pathlib import Path
from scipy.stats import beta,binomtest,t

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def jk(nums,dens):
 n=len(nums);sn,sd=sum(nums),sum(dens);center=sn/sd;loo=[(sn-x)/(sd-y) for x,y in zip(nums,dens)];lm=st.mean(loo)
 return n*center-(n-1)*lm,math.sqrt((n-1)/n*sum((x-lm)**2 for x in loo)),[n*center-(n-1)*x for x in loo]
def effective(x):
 n=len(x);y=[a-st.mean(x) for a in x];g=sum(a*a for a in y)/n
 if g<=2.220446049250313e-16*max(1.,max(a*a for a in x)):return float(n)
 tau=1.
 for k in range(1,min(n//4,16)+1):
  r=sum(y[i]*y[i+k] for i in range(n-k))/((n-k)*g)
  if r<=0:break
  tau+=2*r
 return n/tau
def rhat(chains):
 h=len(chains[0])//2;split=[x[:h] for x in chains]+[x[h:] for x in chains];W=st.mean(st.variance(x) for x in split);B=h*st.variance(st.mean(x) for x in split);tiny=1e-24*max(1.,max(x*x for a in split for x in a))
 return (1. if B<=tiny else math.inf) if W<=tiny else math.sqrt(max(1.,((h-1)/h*W+B/h)/W))
def close(x,y):
 x,y=float(x),float(y)
 if not math.isfinite(x) or not math.isfinite(y):return x==y or (math.isnan(x) and math.isnan(y))
 return abs(x-y)<=2e-10*max(1.,abs(x),abs(y))
def audit_rep(row,cfg,seeds,mode,ref):
 rid=row['id'];assert row['mode']==mode;raw=row['raw_bins'];assert len(raw)==cfg['chains_per_replicate']
 means=[];ses=[];ses2=[];eff=[];pseudo=[]
 for j,x in enumerate(raw):
  assert x['chain']==j+1 and x['seed']==seeds[mode][str(rid)][j]
  nums,dens=x['numerator_bins'],x['denominator_bins'];assert len(nums)==len(dens)==cfg['bins']
  assert all(math.isfinite(v) for v in nums+dens)
  v,se,p=jk(nums,dens);v2,se2,_=jk([sum(nums[k:k+2]) for k in range(0,len(nums),2)],[sum(dens[k:k+2]) for k in range(0,len(dens),2)])
  means.append(v);ses.append(se);ses2.append(se2);eff.append(effective([a-v*b for a,b in zip(nums,dens)]));pseudo.append(p)
  for actual,recorded in [(v,row['chain_means'][j]),(se,row['chain_block_se'][j]),(se2,row['chain_double_block_se'][j]),(eff[-1],row['chain_effective_bins'][j])]:assert close(actual,recorded)
 value=st.mean(means);se=st.stdev(means)/math.sqrt(len(means));half=cfg['student_t7_critical']*se;rh=rhat(pseudo);change=max(abs(a-b)/max(a,1e-12) for a,b in zip(ses,ses2));total=sum(eff)
 good=all(math.isfinite(x) for x in [value,se,half,rh,change,total]) and rh<=cfg['split_rhat_max'] and change<=cfg['blocking_relative_change_max'] and total>=cfg['effective_bins_min'];covered=abs(value-ref)<=half
 for actual,name in [(value,'value'),(se,'standard_error'),(half,'half_width_95'),(rh,'split_rhat'),(change,'blocking_max_relative_change')]:assert close(actual,row[name]),f'Replicate{rid}:{name}'
 assert good==row['quality_pass'] and covered==row['covered']
 return {'id':rid,'value':value,'standard_error':se,'half_width_95':half,'covered':covered,'quality_pass':good,'split_rhat':rh if math.isfinite(rh) else str(rh),'blocking_max_relative_change':change if math.isfinite(change) else str(change),'effective_bins':total if math.isfinite(total) else str(total),'elapsed_seconds':row['elapsed_seconds']}
def main():
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['pilot','final'],required=True);p.add_argument('--inputs',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();own=Path(__file__).parent;cfgp=own/'manifest.toml';cfg=tomllib.loads(cfgp.read_text());seeds=tomllib.loads((own/'seeds.toml').read_text());sources=tomllib.loads((own/'source_hashes.toml').read_text())['files'];ref=cfg['oracle_reference'];n=cfg['pilot_replicates' if a.mode=='pilot' else 'final_replicates'];frozen=json.loads((own/'frozen_inputs.json').read_text())
 for name,digest in frozen['files'].items():assert sha(own/name)==digest,('Frozen input changed',name)
 assert abs(t.ppf(.975,7)-cfg['student_t7_critical'])<1e-12
 out={'id':'X06','schema':cfg['schema'],'mode':a.mode,'status':'UNRESOLVED','manifest_sha256':sha(cfgp),'inputs':{str(x):sha(x) for x in a.inputs},'analysis_source_sha256':sha(__file__),'scipy_version':importlib.metadata.version('scipy'),'reference':ref,'scope':'Single interacting hopping fixture and actual eight-chain95% interval; no universalcoverage or Bonferroni-tail calibration claim.'}
 try:
  rows=[];seen=set();metas=[]
  for inp in a.inputs:
   data=tomllib.loads(inp.read_text());meta=data['metadata'];metas.append(meta)
   assert meta['mode']==a.mode and meta['manifest_sha256']==sha(cfgp) and meta['seeds_sha256']==sha(own/'seeds.toml') and meta['driver_sha256']==sha(own/'run_coverage.jl')
   assert meta['oracle_sha256']==cfg['oracle_sha256'] and meta['source_hashes']==sources and meta['reference']==ref and meta['blas_threads']==1 and meta['julia_version']=='1.12.1'
   assert data['completion']['replicate_count']==meta['last_replicate']-meta['first_replicate']+1
   assert len(data['replicates'])==data['completion']['replicate_count']
   assert [r['id'] for r in data['replicates']]==list(range(meta['first_replicate'],meta['last_replicate']+1))
   for row in data['replicates']:
    assert row['id'] not in seen;seen.add(row['id']);rows.append(audit_rep(row,cfg,seeds,a.mode,ref))
  assert seen==set(range(1,n+1)),f'Need exactly1..{n}; missing{sorted(set(range(1,n+1))-seen)}'
  rows.sort(key=lambda x:x['id']);out['replicates']=rows;out['metadata']=metas;out['raw_audit']='PASS';out['quality_failures']=[r['id'] for r in rows if not r['quality_pass']];out['mean_sampling_seconds_per_replicate']=st.mean(r['elapsed_seconds'] for r in rows)
  if a.mode=='pilot':
   out['status']='PILOT_ONLY';out['planning_full400_sampling_seconds']=400*out['mean_sampling_seconds_per_replicate'];out['planning_full400_seconds_at3x_pilot']=3*out['planning_full400_sampling_seconds'];out['planning_wall_cap_seconds_per_shard']=60*cfg['resources']['full_wall_cap_minutes']
  else:
   k=sum(r['covered'] for r in rows);test=binomtest(k,n,cfg['coverage_null'],alternative='two-sided');alpha=1-cfg['coverage_cp_confidence'];cp=[0. if k==0 else beta.ppf(alpha/2,k,n-k+1),1. if k==n else beta.ppf(1-alpha/2,k+1,n-k)]
   errors=[r['value']-ref for r in rows];bias=st.mean(errors);bias_se=st.stdev(errors)/math.sqrt(n);bh=t.ppf((1+cfg['bias_interval_confidence'])/2,n-1)*bias_se;ci=[bias-bh,bias+bh];margin=cfg['bias_equivalence_margin'];coverage_ok=test.pvalue>=cfg['binomial_significance'] and cp[0]>=cfg['coverage_equivalence_lower'] and cp[1]<=cfg['coverage_equivalence_upper'];bias_ok=ci[0]>=-margin and ci[1]<=margin
   out.update(covered_replicates=k,total_replicates=n,observed_coverage=k/n,binomial_pvalue=float(test.pvalue),coverage_clopper_pearson99=[float(x) for x in cp],coverage_gate_pass=bool(coverage_ok),mean_bias=bias,bias_standard_error=bias_se,mean_bias_interval99=[float(x) for x in ci],bias_equivalence_gate_pass=bool(bias_ok))
   if out['quality_failures']:out['reason']='Quality gate failed; all400 intervals retained without filtering or replacement'
   elif test.pvalue<cfg['binomial_significance'] or ci[0]>margin or ci[1]<-margin:out['status']='FAIL';out['reason']='Resolved coverage or bias discrepancy at the predeclared gate'
   elif coverage_ok and bias_ok:out['status']='PASS'
   else:out['reason']='Coverage/bias equivalence precision unresolved at fixed400replicates'
 except Exception as ex:out['reason']=f'{type(ex).__name__}: {ex}'
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with a.output.open('x') as f:json.dump(out,f,indent=2,allow_nan=False)
 print(json.dumps({k:v for k,v in out.items() if k not in ['replicates','metadata','inputs']},indent=2));raise SystemExit(1 if out['status']=='FAIL' else 0)
if __name__=='__main__':main()

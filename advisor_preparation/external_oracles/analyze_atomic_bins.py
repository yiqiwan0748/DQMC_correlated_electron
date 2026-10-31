#!/usr/bin/env python3
"""Read-only diagnosis of the saved atomic MC bins; launches no sampler."""
import hashlib,json,math,statistics as st,tomllib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
INPUT=ROOT.parent/'isaac/results/external_atomic_mc_4096.toml'
d=tomllib.loads(INPUT.read_text())
result={'input_sha256':hashlib.sha256(INPUT.read_bytes()).hexdigest(),'scientific_runs_started':0,'diagnostic_scope':'64 saved bin means per seed, each averaging64 sweeps. Sweep-level autocorrelation/ESS cannot be reconstructed.','fixtures':{}}
for cid,U,beta,mu in [('ATOM-01',2.,1.,0.),('ATOM-05',4.,.5,.3)]:
 c=next(c for c in d['cases'] if c['id']==cid)
 arm=next(a for a in c['arms'] if a['method']=='DQMC')
 chains=[];series=[];max_energy_identity=0.;max_density_identity=0.;max_offdiag=0.
 for r in arm['raw_bins']:
  x=[v[1]/den for v,den in zip(r['numerator_bins'],r['denominator_bins'])]
  series.append(x);mean=st.mean(x);var=st.variance(x)
  acfs={str(lag):sum((x[i]-mean)*(x[i+lag]-mean) for i in range(len(x)-lag))/((len(x)-lag)*var) for lag in (1,2,3,4)}
  blockse={}
  for size in (1,2,4,8):
   blocked=[st.mean(x[i:i+size]) for i in range(0,len(x),size)]
   blockse[str(size)]=st.stdev(blocked)/math.sqrt(len(blocked))
  chains.append({'seed':r['seed'],'energy_mean':mean,'energy_within_seed_SE':math.sqrt(var/len(x)),'early_half_mean':st.mean(x[:32]),'late_half_mean':st.mean(x[32:]),'min_bin_mean':min(x),'max_bin_mean':max(x),'acf_of_saved_bins':acfs,'contiguous_block_SE':blockse})
  for raw,den in zip(r['numerator_bins'],r['denominator_bins']):
   v=[q/den for q in raw];Eunc,Ecenter,N,D=v[:4]
   max_energy_identity=max(max_energy_identity,abs(Eunc-U*D),abs(Ecenter-U*(D-N/2+1)))
   max_density_identity=max(max_density_identity,abs(N-4-(1-math.exp(-2*beta*mu))*D))
   max_offdiag=max(max_offdiag,abs(v[6]),abs(v[7]))
 n=len(series[0]);m=len(series);means=[st.mean(x) for x in series];W=st.mean(st.variance(x) for x in series);B=n*st.variance(means)
 Rhat=math.sqrt(((n-1)/n*W+B/n)/W)
 weights=[1/ch['energy_within_seed_SE']**2 for ch in chains];weighted=sum(w*x for w,x in zip(weights,means))/sum(weights)
 Q=sum(w*(x-weighted)**2 for w,x in zip(weights,means))
 result['fixtures'][cid]={'chains':chains,'classic_Rhat_on_saved_bins':Rhat,'between_to_within_seed_mean_variance_ratio':st.variance(means)/(W/n),'inverse_variance_heterogeneity_Q':Q,'heterogeneity_df':m-1,'reported_effective_bins':arm.get('estimated_effective_bins'),'reported_between_seed_SE':arm['standard_error'],'naive_within_bin_pooled_SE':math.sqrt(W/(n*m)),'max_atomic_energy_identity_error':max_energy_identity,'max_atomic_density_D_identity_error':max_density_identity,'max_spurious_offdiagonal_C':max_offdiag}
result['interpretation']=[
 'Second fixture has stable between-seed plateaus. Six seeds have energy near-1.90; seeds4001/5003 near-1.51. Blocking within a chain does not reveal the missing between-mode motion.',
 'The much larger between-seed error is a useful warning; increasing only samples inside the same trapped trajectories is not justified before inspecting transition coverage/initialization/stabilization.',
 'Atomic energy and density-versus-double-occupancy identities hold, and offdiagonal propagators vanish numerically. This supports internal measurement consistency but does not prove the sampled field distribution is correct.',
 'Saved bins average64 sweeps, so a high-frequency periodic defect or unsampled component could be invisible. No sweep-level ESS or exact cause is inferred from these bins.',
 'First fixture is a control with no analogous plateau; differences in time-slice count, stabilization and field dynamics need a fresh targeted diagnostic before declaring a cause.'
]
(ROOT/'atomic_bin_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
for key,value in result['fixtures'].items():
 print(key,'Rhat',value['classic_Rhat_on_saved_bins'],'Q',value['inverse_variance_heterogeneity_Q'],'between/within',value['between_to_within_seed_mean_variance_ratio'],'energy identity max',value['max_atomic_energy_identity_error'])

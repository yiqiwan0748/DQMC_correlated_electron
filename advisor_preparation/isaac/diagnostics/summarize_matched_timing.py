#!/usr/bin/env python3
"""Turn a matched timing receipt into an auditable public summary. No benchmark runs."""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics,tomllib
from pathlib import Path
from datetime import datetime,timezone
EXPECTED_SIZES=[8,16,32,64,128]
FROZEN_BENCHMARK_SHA256='56ce5e3c395df4aeb1fe0b7654e315a81e280c74b2f65bb0b0ab65a010d06c49'
LARGE_CUT_BENCHMARK_SHA256='aefc720f6626f020b455b0eb80d398b96a910b627b17cb1968984ac3bd4c4e8d'
CORRECTED_SOURCE_HASHES={'src/HubbardRenyiDQMC.jl':'336f4393cd69b04ffc1b8c42f97207d719a4ab88cbbd1d343a8088234bffc721','src/dqmc_kernel.jl':'b72bbbf158381b0f776c1fa259d532b069094c47e361d3490c19ea0f43ada4be','src/grover_rank1_update.jl':'79064a0aa868ff7b5c665876bc5dd6eeeee17da2684b7a04807891c81d0d358e'}
DIRECT='without_Grover_cache_direct_fused';OPT='with_rank1_Grover_cache'
BOOL_CHECKS=['sweep_boundary_HS_fields_exact','sweep_acceptance_rates_exact','determinant_signs_exact','next_RNG_draws_exact']
NUM_CHECKS=['maximum_Green_absolute_error','maximum_logdet_absolute_error','maximum_measurement_absolute_error','maximum_measurement_relative_error','maximum_reported_deltaG']
def finite(x):return isinstance(x,(int,float)) and math.isfinite(x)
def median(xs):return statistics.median(xs) if xs else None
def fmt(x,d=6):return f'{x:.{d}g}' if finite(x) else '—'
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,list):return [clean(v) for v in x]
 if isinstance(x,float) and not math.isfinite(x):return None
 return x
def process_row(row,index,duplicates,exploration=False):
 raw=row.get('raw_timings',[]);checks=row.get('correctness_checks',[])
 dr=[r for r in raw if r.get('variant')==DIRECT];op=[r for r in raw if r.get('variant')==OPT]
 dt=[r['node_loop_seconds'] for r in dr if finite(r.get('node_loop_seconds')) and r['node_loop_seconds']>0]
 ot=[r['node_loop_seconds'] for r in op if finite(r.get('node_loop_seconds')) and r['node_loop_seconds']>0]
 dm,om=median(dt),median(ot);ratio=dm/om if dm is not None and om is not None else None
 dk=median([r['update_kernel_seconds'] for r in dr if finite(r.get('update_kernel_seconds')) and r['update_kernel_seconds']>0]);ok=median([r['update_kernel_seconds'] for r in op if finite(r.get('update_kernel_seconds')) and r['update_kernel_seconds']>0]);kernel_ratio=dk/ok if dk is not None and ok is not None else None
 da=median([r['allocated_bytes'] for r in dr if finite(r.get('allocated_bytes'))]);oa=median([r['allocated_bytes'] for r in op if finite(r.get('allocated_bytes'))]);allocation_ratio=da/oa if da is not None and oa is not None and oa>0 else None
 requested=row.get('requested_repeats',row.get('repeats'));completed=row.get('completed_repeats',min(len(dt),len(ot)))
 exact={k:bool(checks) and all(c.get(k) is True for c in checks) for k in BOOL_CHECKS}
 maxima={k:max([c[k] for c in checks if finite(c.get(k))],default=None) for k in NUM_CHECKS}
 registered_pass=(row.get('status')=='PASS' and bool(checks) and all(c.get('status')=='PASS' for c in checks))
 full_repeats=requested==7 and completed==requested and len(dt)==len(ot)==len(dr)==len(op)==len(checks)==requested
 good_numbers=bool(checks) and all(finite(c.get(k)) and c[k]>=0 for c in checks for k in NUM_CHECKS)
 numeric_ok=good_numbers and maxima['maximum_Green_absolute_error']<=2e-8 and maxima['maximum_logdet_absolute_error']<=2e-8 and maxima['maximum_measurement_relative_error']<=2e-8 and maxima['maximum_reported_deltaG']<=1e-7
 pair_ids_d=[r.get('repeat') for r in dr];pair_ids_o=[r.get('repeat') for r in op]
 pairs_valid=len(set(pair_ids_d))==len(dr) and len(set(pair_ids_o))==len(op) and set(pair_ids_d)==set(pair_ids_o)==set(range(1,8)) and {c.get('repeat') for c in checks}==set(range(1,8))
 reported=row.get('ratio_of_median_node_times');ratio_consistent=ratio is not None and finite(reported) and abs(ratio-reported)<=1e-12*max(1.,abs(reported))
 physics={'t':1.,'U':2.,'mu_centered':0.,'beta':2.,'dtau':.05,'lambda':.5,'delta_lambda':.125,'n_stab':5,'boundary':'pbc'}
 cut_ok=(row.get('L')==64 and row.get('nA') in (48,64)) if exploration else row.get('nA')==row.get('L',0)//2
 physics_ok=all(row.get(k)==v for k,v in physics.items()) and cut_ok
 size_ok=row.get('L') in ([64] if exploration else EXPECTED_SIZES)
 eligible=registered_pass and full_repeats and all(exact.values()) and numeric_ok and pairs_valid and ratio_consistent and physics_ok and size_ok and row.get('L') not in duplicates and row.get('two_replica_update_rounds',0)==12
 notes=[]
 if row.get('reason'):notes.append(row['reason'])
 if row.get('L') in duplicates:notes.append('Duplicate size records preserved; not eligible for a headline claim.')
 if not size_ok:notes.append('Size is outside the declared campaign grid; not eligible for its performance claim.')
 if registered_pass and not full_repeats:notes.append('A pilot or incomplete repeat set; not a production timing claim.')
 if row.get('status')=='PASS' and not all(exact.values()):notes.append('Missing or failed exact replay checks.')
 if row.get('status')=='PASS' and not physics_ok:notes.append('Physics parameters do not match the declared campaign protocol.')
 if ratio is not None and not ratio_consistent:notes.append('Reported ratio does not match medians recomputed from raw repeats.')
 interval=row.get('paired_bootstrap_95pct_ratio_interval',[None,None]);interval=(interval+[None,None])[:2]
 return dict(row_index=index,L=row.get('L'),nA=row.get('nA',row.get('L',0)//2),status=row.get('status','UNRESOLVED'),requested_repeats=requested,completed_repeats=completed,two_replica_rounds=row.get('two_replica_update_rounds'),median_direct_node_seconds=dm,median_optimized_node_seconds=om,node_loop_speedup=ratio,reported_node_loop_speedup=reported,ratio_recomputed_consistent=ratio_consistent,bootstrap_ratio_low=interval[0],bootstrap_ratio_high=interval[1],update_kernel_speedup=kernel_ratio,reported_update_kernel_speedup=row.get('ratio_of_median_update_kernel_times'),median_direct_kernel_seconds=dk,median_optimized_kernel_seconds=ok,median_direct_measurement_seconds=median([r['measurement_seconds'] for r in dr if finite(r.get('measurement_seconds'))]),median_optimized_measurement_seconds=median([r['measurement_seconds'] for r in op if finite(r.get('measurement_seconds'))]),median_direct_allocated_bytes=da,median_optimized_allocated_bytes=oa,allocation_reduction_ratio=allocation_ratio,median_direct_gc_seconds=median([r['gc_seconds'] for r in dr if finite(r.get('gc_seconds'))]),median_optimized_gc_seconds=median([r['gc_seconds'] for r in op if finite(r.get('gc_seconds'))]),initialization_seconds_excluded=row.get('initialization_seconds_excluded'),common_warmup_seconds_excluded=row.get('common_warmup_seconds_excluded'),exact_replay_checks=exact,correctness_maxima=maxima,production_claim_eligible=eligible,notes=notes)
def plot_results(rows,out):
 try:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
 except ImportError as ex:return {'status':'UNAVAILABLE','reason':str(ex)}
 fig,(ax,bx)=plt.subplots(1,2,figsize=(9.2,3.7),layout='constrained')
 good=[r for r in rows if finite(r['median_direct_node_seconds']) and finite(r['median_optimized_node_seconds'])]
 for key,label,col in [('median_direct_node_seconds','Direct determinant recomputation','#B34A3C'),('median_optimized_node_seconds','Cached rank-one update','#245E9B')]:
  gr=sorted(good,key=lambda r:r['L']);ax.plot([r['L'] for r in gr],[r[key] for r in gr],'o-',label=label,color=col,linewidth=1.5,markersize=5)
 ax.set(xscale='log',yscale='log',xlabel='Chain length L (subsystem nA=L/2)',ylabel='Median matched node-loop time (s)');ax.grid(alpha=.2,which='both');ax.legend(fontsize=8)
 for r in rows:
  if finite(r['node_loop_speedup']):
   color='#245E9B' if r['production_claim_eligible'] else '#888888';x=r['L'];y=r['node_loop_speedup'];lo=r['bootstrap_ratio_low'];hi=r['bootstrap_ratio_high']
   if finite(lo) and finite(hi) and lo<=y<=hi:bx.errorbar([x],[y],yerr=[[y-lo],[hi-y]],fmt='o',color=color,capsize=3)
   else:bx.plot(x,y,'o',color=color)
   bx.annotate(f'{y:.2f}×',(x,y),xytext=(0,6),textcoords='offset points',ha='center',fontsize=8)
  else:bx.annotate(f'L={r["L"]}: unresolved',(0.03,.95-.07*rows.index(r)),xycoords='axes fraction',fontsize=8,color='#666666')
 bx.axhline(1,color='#777777',linestyle='--',linewidth=1);bx.set(xscale='log',xlabel='Chain length L',ylabel='Direct / optimized node-loop time');bx.grid(alpha=.2);bx.set_title('All registered sizes retained',fontsize=10)
 fig.suptitle('Matched Rényi node loop: same local proposals, fields and measurements',fontsize=11)
 png=out/'timing_scaling.png';pdf=out/'timing_scaling.pdf';fig.savefig(png,dpi=220);fig.savefig(pdf);plt.close(fig)
 return {'status':'GENERATED','png':png.name,'pdf':pdf.name,'interpretation':'Error bars reflect timing jitter from paired repeat resampling, not uncertainty across hardware.'}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('receipt',type=Path);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--force',action='store_true');ap.add_argument('--no-plot',action='store_true');ap.add_argument('--campaign-mode',choices=['primary','cut-exploration'],default='primary');a=ap.parse_args();exploration=a.campaign_mode=='cut-exploration'
 receipt=a.receipt.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 targets=['timing_summary.json','timing_summary.csv','timing_repeats.csv','README.md']
 if not a.force and any((out/n).exists() for n in targets):raise SystemExit('Output already exists; choose a new directory or use --force for an intentional regeneration.')
 data=tomllib.loads(receipt.read_text());rawrows=data.get('results',[])
 declared_cuts=dict(zip(data.get('registered_sizes',[]),data.get('subsystem_sizes',[]))) if exploration else {}
 if exploration:rawrows=[{**r,'nA':r.get('nA',declared_cuts.get(r.get('L')))} for r in rawrows]
 sizes=[r.get('L') for r in rawrows];duplicates={L for L in sizes if sizes.count(L)>1};rows=[process_row(r,i,duplicates,exploration) for i,r in enumerate(rawrows,1)]
 registered=data.get('registered_sizes',[]);missing=[L for L in registered if L not in sizes]
 for L in missing:rows.append(dict(row_index=None,L=L,nA=declared_cuts.get(L,L//2),status='UNRESOLVED',requested_repeats=None,completed_repeats=0,two_replica_rounds=None,median_direct_node_seconds=None,median_optimized_node_seconds=None,node_loop_speedup=None,reported_node_loop_speedup=None,ratio_recomputed_consistent=False,bootstrap_ratio_low=None,bootstrap_ratio_high=None,update_kernel_speedup=None,exact_replay_checks={k:False for k in BOOL_CHECKS},correctness_maxima={k:None for k in NUM_CHECKS},production_claim_eligible=False,notes=['Registered size missing from receipt; retained as unresolved.']))
 # The selection rule is largest completed registered size, never maximum observed speedup.
 protocol_metadata_checks={'declared_campaign_grid':registered==([64] if exploration else EXPECTED_SIZES),'one_Julia_thread':data.get('julia_threads')==1,'one_BLAS_thread':data.get('blas_threads')==1,'frozen_benchmark_source':data.get('runner_sha256')==(LARGE_CUT_BENCHMARK_SHA256 if exploration else FROZEN_BENCHMARK_SHA256),'source_hashes_present':all(isinstance(data.get('source_hashes',{}).get(p),str) and len(data['source_hashes'][p])==64 for p in ['Project.toml','Manifest.toml','src/HubbardRenyiDQMC.jl','src/dqmc_kernel.jl','src/grover_rank1_update.jl']),'corrected_kernel_source_hashes':all(data.get('source_hashes',{}).get(p)==digest for p,digest in CORRECTED_SOURCE_HASHES.items()),'registered_seeds':data.get('seeds')==[940031,940032]}
 if exploration:protocol_metadata_checks.update({'exploration_class_declared':data.get('campaign_class')=='post_primary_subsystem_exploration','declared_exploratory_cut':data.get('subsystem_sizes') in ([48],[64]) and all(r.get('nA')==declared_cuts.get(r.get('L')) for r in rawrows),'original_primary_runner_identified':data.get('primary_campaign_runner_sha256')==FROZEN_BENCHMARK_SHA256})
 for r in rows:
  if r['production_claim_eligible'] and not all(protocol_metadata_checks.values()):r['notes'].append('Global protocol metadata checks failed; production claim suppressed.')
  r['production_claim_eligible']=r['production_claim_eligible'] and all(protocol_metadata_checks.values()) and r['L'] in registered
 eligible=[r for r in rows if r['production_claim_eligible']];headline=max(eligible,key=lambda r:r['L']) if eligible and not exploration else None
 metadata={k:data.get(k) for k in ['created_utc','finished_utc','comparison','scope','primary_statistic','timing_uncertainty','source_root','julia_version','host','cpu','julia_threads','blas_threads','blas_config','slurm_job_id','seeds','registered_sizes','soft_budget_seconds','runner_sha256','source_hashes','compilation_and_initial_probe_seconds_excluded','campaign_class','subsystem_sizes','cut_selection','primary_campaign_runner_sha256','selection_disclosure','entropy_scope']}
 globalmax={k:max([r['correctness_maxima'][k] for r in rows if finite(r['correctness_maxima'][k])],default=None) for k in NUM_CHECKS}
 status_counts={s:sum(r['status']==s for r in rows) for s in sorted(set(r['status'] for r in rows))}
 cv=None
 if headline:
  cv=f"Optimized Rényi entropy sampling with cached rank-one Grover updates, achieving {headline['node_loop_speedup']:.1f}× faster matched node-loop execution at L={headline['L']} (nA={headline['nA']}) versus direct determinant recomputation on ISAAC."
 try:figure={'status':'DISABLED'} if a.no_plot else plot_results(rows,out)
 except Exception as ex:figure={'status':'UNAVAILABLE','reason':f'{type(ex).__name__}: {ex}'}
 summary=dict(generated_utc=datetime.now(timezone.utc).isoformat(),campaign_mode=a.campaign_mode,receipt_path=str(receipt),receipt_sha256=hashlib.sha256(receipt.read_bytes()).hexdigest(),analysis_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),metadata=metadata,protocol_metadata_checks=protocol_metadata_checks,status_counts=status_counts,rows=rows,all_correctness_maxima=globalmax,headline_selection_rule='Exploration never replaces the primary headline; all explored cuts must be disclosed.' if exploration else 'Largest completed passing size in the full preregistered [8,16,32,64,128] grid; never select maximum observed speedup.',headline_case=headline,exploratory_comparisons=eligible if exploration else [],cv_candidate=cv,figure=figure,limitations=['One lambda=.5 node loop, not the complete final-entropy lambda grid.','Timing is host/workload specific; repeat bootstrap does not represent hardware-independent uncertainty.','This measures Grover determinant caching and removal of redundant candidate matrix work, not an isolated ordinary equal/unequal-time Green-function speedup.','Initialization, common warmup, compilation, deep copies and explicit pre-run GC are excluded; allocation-triggered GC inside the loop is included.','Unresolved/failed/pilot cases remain in the public summary and cannot supply a production CV number.']+(['Post-primary cut-size exploration, designed after seeing the original scaling results.','The nA=L cut is full-system finite-temperature thermal Renyi purity; it is not pure-state bipartite entanglement.'] if exploration else []))
 (out/'timing_summary.json').write_text(json.dumps(clean(summary),indent=2,ensure_ascii=False)+'\n')
 columns=['row_index','L','nA','status','requested_repeats','completed_repeats','two_replica_rounds','median_direct_node_seconds','median_optimized_node_seconds','node_loop_speedup','bootstrap_ratio_low','bootstrap_ratio_high','median_direct_kernel_seconds','median_optimized_kernel_seconds','update_kernel_speedup','median_direct_measurement_seconds','median_optimized_measurement_seconds','median_direct_allocated_bytes','median_optimized_allocated_bytes','allocation_reduction_ratio','median_direct_gc_seconds','median_optimized_gc_seconds','production_claim_eligible','notes']
 columns+=BOOL_CHECKS+NUM_CHECKS
 with (out/'timing_summary.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');w.writeheader()
  for r in rows:w.writerow({**r,**r['exact_replay_checks'],**r['correctness_maxima'],'notes':' '.join(r['notes'])})
 repeatcols=['row_index','L','nA','case_status','variant','repeat','execution_order','node_loop_seconds','update_kernel_seconds','measurement_seconds','other_loop_seconds','allocated_bytes','gc_seconds']
 with (out/'timing_repeats.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=repeatcols,extrasaction='ignore');w.writeheader()
  for i,r in enumerate(rawrows,1):
   for t in r.get('raw_timings',[]):w.writerow({'row_index':i,'L':r.get('L'),'nA':r.get('nA'),'case_status':r.get('status'),**t})
 lines=['# Matched Rényi node-loop timing','', 'This record compares the same locally lazy two-replica workflow **with and without rank-one Grover caching**. The primary metric is the ratio of median elapsed node-loop times; update-kernel results remain separate.','']
 if exploration:
  lines+=['**Post-primary cut-size exploration.** This campaign was designed after seeing the original half-chain results and does not replace that preregistered grid. Both exploratory cuts must remain visible in the combined public record. A full-system cut (nA=L) describes finite-temperature thermal Rényi purity, not pure-state bipartite entanglement. No primary headline or automatic CV candidate is selected from this exploration.','']
  for r in eligible:lines += [f"Completed exploratory case **L={r['L']}, nA={r['nA']}**: **{r['node_loop_speedup']:.3f}×** node-loop speedup ({fmt(r['median_direct_node_seconds'])} s direct / {fmt(r['median_optimized_node_seconds'])} s optimized).",'']
 elif headline:lines += [f"The largest completed passing registered case, **L={headline['L']}, nA={headline['nA']}**, ran **{headline['node_loop_speedup']:.3f}× faster**: {fmt(headline['median_direct_node_seconds'])} s direct versus {fmt(headline['median_optimized_node_seconds'])} s optimized. All sizes are shown below.",'']
 else:lines += ['No production CV number is certified by this receipt: a complete, passing registered comparison is required.','']
 lines += ['| L | nA | Status | Paired repeats | Direct median (s) | Optimized median (s) | Node-loop speedup | Timing-jitter 95% interval | Update-kernel speedup |','|---:|---:|---|---:|---:|---:|---:|---|---:|']
 for r in rows:lines.append(f"| {r['L']} | {r['nA']} | {r['status']} | {r['completed_repeats']}/{r['requested_repeats'] or '—'} | {fmt(r['median_direct_node_seconds'])} | {fmt(r['median_optimized_node_seconds'])} | {fmt(r['node_loop_speedup'],5)}× | [{fmt(r['bootstrap_ratio_low'],5)}, {fmt(r['bootstrap_ratio_high'],5)}] | {fmt(r['update_kernel_speedup'],5)}× |")
 notes=[f"- L={r['L']}: {' '.join(r['notes'])}" for r in rows if r['notes']]
 if notes:lines+=['','Retained case notes:','',*notes]
 failed_metadata=[k for k,v in protocol_metadata_checks.items() if not v]
 if failed_metadata:lines+=['','Production metadata checks requiring attention: '+', '.join(failed_metadata)+'.']
 lines+=['','## Allocations','','These are median cumulative allocated bytes inside each timed node loop, including allocation churn. They do not measure peak resident memory. The reduction ratio is direct median bytes divided by optimized median bytes.','','| L | nA | Direct median allocated bytes | Optimized median allocated bytes | Allocation reduction |','|---:|---:|---:|---:|---:|']
 for r in rows:lines.append(f"| {r['L']} | {r['nA']} | {fmt(r.get('median_direct_allocated_bytes'),10)} | {fmt(r.get('median_optimized_allocated_bytes'),10)} | {fmt(r.get('allocation_reduction_ratio'),5)}× |")
 lines += ['', '## Scope and fairness','', 'Fixed physics: nearest-neighbor periodic Hubbard chain, t=1, U=2, centered mu=0, beta=2, dtau=0.05, lambda=0.5, nA=L/2 and stabilization period 5. Both paths use the corrected 50% independent idle decision before every local proposal. Both start from clones of identical warmed states and RNG streams; execution order alternates. The production protocol uses seven repeats of twelve two-replica rounds. The same forward/backward and endpoint ratio measurements and sweep-boundary snapshots are included.','', 'Initialization, common warmup, compilation, deep copies and explicit GC are outside the primary timer. Allocation-triggered GC inside the measured loop is included. Raw allocation, GC and component timings are preserved in [timing_repeats.csv](timing_repeats.csv). The bootstrap interval measures repeat-to-repeat timing jitter on this host, not performance variation across machines.','', 'The direct baseline rebuilds subsystem Grover matrices and refactors determinants for proposals. The optimized path maintains their determinant and inverse, uses rank-one updates, avoids redundant full candidate-state copying/rebuilding, and reconstructs at registered stabilization boundaries with conditioning safeguards. This is a Grover-update measurement; it does not establish an isolated speedup for ordinary equal/unequal-time Green functions.','', '## Correctness','', 'Every accepted timing case must reproduce sweep-boundary auxiliary fields, sweep acceptance rates, determinant signs and next RNG draws exactly. Full Green matrices, determinant logs and all ratio observations are checked numerically. Maximum discrepancies across every retained repeat are:','']
 if exploration:lines=[line.replace('nA=L/2 and stabilization period 5','nA as declared in the table and stabilization period 5') for line in lines]
 for key,value in globalmax.items():lines.append(f'- `{key}`: {fmt(value,8)}')
 lines += ['', 'These checks establish matched sampled trajectories at the recorded sweep boundaries. They are not an independent proof of the final entropy’s physical accuracy, which belongs to the external scientific benchmark campaign. Failed and unresolved sizes remain visible.','', '## Provenance','',f"- Recorded UTC: {metadata.get('created_utc')} to {metadata.get('finished_utc')}",f"- Host / Slurm job: `{metadata.get('host')}` / `{metadata.get('slurm_job_id')}`",f"- CPU: `{metadata.get('cpu')}`",f"- Julia / threads / BLAS threads: `{metadata.get('julia_version')}` / `{metadata.get('julia_threads')}` / `{metadata.get('blas_threads')}`",f"- Seeds: `{metadata.get('seeds')}`",f"- Receipt SHA256: `{summary['receipt_sha256']}`",f"- Benchmark runner SHA256: `{metadata.get('runner_sha256')}`",f"- Analysis script SHA256: `{summary['analysis_script_sha256']}`",'', 'Source hashes:','']
 for name,digest in (metadata.get('source_hashes') or {}).items():lines.append(f'- `{name}`: `{digest}`')
 lines += ['',f"BLAS configuration: `{metadata.get('blas_config')}`",'', 'Full metadata and every row are in [timing_summary.json](timing_summary.json); compact values are in [timing_summary.csv](timing_summary.csv).']
 if figure['status']=='GENERATED':lines += ['', '![Matched timing and scaling](timing_scaling.png)']
 else:lines += ['',f"Optional figure status: {figure['status']}. {figure.get('reason','')}"]
 if cv:lines += ['', '## Scoped CV wording','',cv,'', 'This wording identifies the tested algorithm and workload. The largest completed registered size was selected by size, not by searching for the largest observed ratio. The full table must accompany the public evidence record.']
 (out/'README.md').write_text('\n'.join(lines)+'\n')
 print(json.dumps({'status_counts':status_counts,'headline':None if headline is None else {'L':headline['L'],'node_loop_speedup':headline['node_loop_speedup']},'figure':figure,'output_dir':str(out)},indent=2))
if __name__=='__main__':main()

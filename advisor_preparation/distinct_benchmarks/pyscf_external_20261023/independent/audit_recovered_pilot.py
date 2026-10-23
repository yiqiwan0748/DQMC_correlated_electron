#!/usr/bin/env python3
"""Independent immutable-input/linear-algebra/resource audit; no physics solver."""
import argparse,hashlib,json,math,re,tomllib
from pathlib import Path
import numpy as np

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return tomllib.loads(p.read_text()) if p.suffix=='.toml' else json.loads(p.read_text())
def coo(x):
 a=np.zeros(x['shape']);a[x['row'],x['col']]=x['value'];return a

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--receipt-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();p=a.receipt_root
 assert not a.output.exists();checks=[]
 def need(n,b):checks.append({'name':n,'passed':bool(b)});assert b,n
 def close(n,x,y,tol=1e-10):
  e=float(np.max(np.abs(np.asarray(x)-np.asarray(y))));checks.append({'name':n,'error':e,'limit':tol,'passed':math.isfinite(e) and e<=tol});assert checks[-1]['passed'],n
 before=load(p/'transfer_before.json');after=load(p/'transfer_after.json');need('remote after PASS',after['status']=='PASS')
 aft={r['local']:r for r in after['files']}
 for row in before['files']:
  q=p/row['local'];need('bytes/hash '+row['local'],q.stat().st_size==row['bytes'] and sha(q)==row['sha256']);r=aft[row['local']];need('remote stable '+row['local'],r['sha256']==row['sha256'] and r['bytes']==row['bytes'] and r['mtime_ns']==row['before']['mtime_ns'] and r['inode']==row['before']['inode'])
 for role in ('repair_source','original_source'):
  f=load(p/role/'source_freeze.json')
  for n,h in f['files'].items():need('frozen '+role+'/'+n,sha(p/role/n)==h)
 original=load(p/'original_results/execution.json');need('original failure retained',original['status']=='UNRESOLVED' and [r['exit_code'] for r in original['processes']]==[0,0,1])
 cfg=load(p/'original_source/manifest.toml');comp=load(p/'results/comparison.json');replay=load(p/'independent_replay/comparison.json')
 replay_max=0.
 def equivalent(x,y):
  nonlocal replay_max
  if isinstance(x,float):
   z=abs(x-y);replay_max=max(replay_max,z);return math.isfinite(z) and z<=1e-12
  if isinstance(x,dict):return x.keys()==y.keys() and all(equivalent(v,y[k]) for k,v in x.items())
  if isinstance(x,list):return len(x)==len(y) and all(equivalent(v,w) for v,w in zip(x,y))
  return x==y
 for key in comp:
  if key not in ('input_sha256','resource_warm_postprocessing_seconds'):need('frozen replay metadata/gates and roundoff agreement '+key,equivalent(comp[key],replay[key]))
 need('replay input hashes',sorted(comp['input_sha256'].values())==sorted(replay['input_sha256'].values()))
 need('zero credits',comp['completion_count']==0 and comp['pilot_excluded'] and [r['status'] for r in comp['rows']]==['NOT_RUN','NOT_RUN'])
 need('all comparator gates',len(comp['checks'])==529 and all(r['passed'] for r in comp['checks']))
 repair=load(p/'results/repair_receipt.json');need('strict bool-only recovery',repair['numpy_bool_leaves_converted']==2 and repair['other_unsupported_types_converted']==0 and repair['original_source_and_raw_unchanged'] and repair['comparator_exit_code']==0 and not repair['physics_executed'])
 t=load(p/'original_results/target.toml');q=load(p/'original_results/reference.json');need('pilot fixture grid',[r['id'] for r in t['cases']]==[r['id'] for r in q['cases']]==cfg['pilot_case_ids'])
 summaries=[]
 for tr,qr in zip(t['cases'],q['cases']):
  pairs=[]
  for tag,raw in [('target',tr),('reference',qr)]:
   coords=list(zip(raw['basis_up'],raw['basis_down']));ix=sorted(range(len(coords)),key=coords.__getitem__);E=np.asarray(raw['energies']);V=np.asarray(raw['vectors'])[ix];H=coo(raw['H'])[np.ix_(ix,ix)];C=coo(raw['C'])[np.ix_(ix,ix)];S=coo(raw['S2'])[np.ix_(ix,ix)]
   close(tag+'/'+raw['id']+'/eigen residual',H@V,V*E);close(tag+'/'+raw['id']+'/orthogonality',V.T@V,np.eye(len(E)))
   d=raw['ground_dimension'];A=V.T@C@V[:,:d];B=V.T@C.T@V[:,:d]
   close(tag+'/'+raw['id']+'/directed amplitudes',A,raw['transition_amplitudes']);close(tag+'/'+raw['id']+'/reverse amplitudes',B,raw['reverse_transition_amplitudes'])
   groups=[];start=0
   for j in range(1,len(E)+1):
    if j==len(E) or E[j]-E[start]>cfg['gates']['cluster_diameter_max']:groups.append(np.arange(start,j));start=j
   need(tag+'/'+raw['id']+'/eight complete manifolds',len(groups)>=8)
   pairs.append((E,V,H,C,S,A,B,groups))
   if raw['id']=='DIR_N2_U1':
    close(tag+'/witness energy',E[0],(1-math.sqrt(50))/2);close(tag+'/witness raising total',np.sum(A*A)/d,2);close(tag+'/witness lowering total',np.sum(B*B)/d,0)
  x,y=pairs
  for idx,name in [(0,'energies'),(2,'H'),(3,'C'),(4,'S2')]:close(tr['id']+'/'+name,x[idx],y[idx],1e-8)
  for k,(I,J) in enumerate(zip(x[7][:8],y[7][:8])):
   need(tr['id']+'/manifold dimension '+str(k),len(I)==len(J));close(tr['id']+'/projector '+str(k),x[1][:,I]@x[1][:,I].T,y[1][:,J]@y[1][:,J].T,1e-8)
   for idx,name in [(5,'raising'),(6,'lowering')]:close(tr['id']+'/'+name+str(k),np.sum(x[idx][I]**2)/tr['ground_dimension'],np.sum(y[idx][J]**2)/qr['ground_dimension'],1e-8)
  summaries.append({'case_id':tr['id'],'manifolds':8,'ground_energy':float(x[0][0]),'ground_dimension':tr['ground_dimension']})
 resource=load(p/'results/resource_assessment.json');proc=load(p/'results/postprocessing_process.json');pol=cfg['resource_policy'];cold={r['name']:r['finished_unix']-r['started_unix'] for r in original['processes']};cold['comparator']+=proc['finished_unix']-proc['started_unix'];project={};peaks={}
 for name,raw in [('target',t),('reference',q)]:
  tm=raw['resource_repeat_cases'][0]['timing_seconds'];di=tm['full_diagonalization'];other=sum(tm.values())-di;R=784/64;project[name]=pol['runtime_safety_factor']*(cold[name]+2*pol['particle_complexity_guard']*sum(tm.values())+3*(pol['particle_complexity_guard']*R**2*other+R**3*di))+pol['fixed_io_allowance_per_process_seconds']
 w=comp['resource_warm_postprocessing_seconds'];project['comparator']=pol['runtime_safety_factor']*(cold['comparator']+2*pol['particle_complexity_guard']*w+3*(784/64)**3*w)+pol['fixed_io_allowance_per_process_seconds']
 def rss(f):return int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',f.read_text())[1])/1024
 peaks={k:rss(p/'original_results'/f'{k}.resources.txt') for k in ('target','reference','comparator')};peaks['comparator']=max(peaks['comparator'],rss(p/'results/postprocessing.resources.txt'))
 extra=pol['bytes_per_retained_dimension_squared']*sum(s['dimension']**2 for s in cfg['cases'])/2**20;projectpeak={k:pol['memory_safety_factor']*(v+extra) for k,v in peaks.items()}
 for n in project:close('independent resource wall '+n,project[n],resource['projected_process_seconds'][n],1e-8);close('independent resource peak '+n,projectpeak[n],resource['projected_peak_mib'][n],1e-8)
 need('same resource policy',resource['thresholds']==pol);need('engineering limits only',sum(project.values())<=2100 and max(projectpeak.values())<=2800 and not resource['production_admitted'] and resource['scientific_completions']==0)
 inputs={str(f.relative_to(p)):sha(f) for f in sorted(p.rglob('*')) if f.is_file() and f!=a.output and f.name!='remote_text_receipts.tar'}
 out={'schema':'PYSCF_RECOVERED_PILOT_INDEPENDENT_AUDIT_V1','audit_status':'PASS','integrity':'PASS','catalogue_ids':['P23','P24'],'catalogue_status':{'P23':'NOT_RUN','P24':'NOT_RUN'},'scientific_completions':0,'physics_executed':False,'check_count':len(checks),'checks':checks,'pilot_cases':summaries,'frozen_comparator_checks':len(comp['checks']),'local_vs_ISAAC_replay_max_float_difference':replay_max,'replay_roundoff_allowance':1e-12,'original_execution_retained':'UNRESOLVED','only_numpy_bool_serialization_repaired':True,'projected_total_seconds':sum(project.values()),'maximum_projected_peak_mib':max(projectpeak.values()),'production_admitted':False,'input_sha256':inputs,'auditor_sha256':sha(__file__)}
 a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({k:out[k] for k in ['audit_status','check_count','scientific_completions','projected_total_seconds','maximum_projected_peak_mib','auditor_sha256']}))
if __name__=='__main__':main()

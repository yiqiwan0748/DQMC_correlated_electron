#!/usr/bin/env python3
"""Exact32-state HS scan analysis. No user implementation or MC simulation."""
import json,math,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
M,U,beta,mu=5,4.,.5,.3
dtau=beta/M;lam=math.acosh(math.exp(dtau*U/2));x=math.exp(beta*mu)
def spin_sum(state):return sum(1 if (state>>i)&1 else -1 for i in range(M))
def weight(state):
 s=spin_sum(state)
 return (1+x*math.exp(lam*s))*(1+x*math.exp(-lam*s))
weights=[weight(s) for s in range(1<<M)];Z=sum(weights);pi=[w/Z for w in weights]
def transition(proposal_probability,reverse=False):
 matrix=[]
 for state in range(1<<M):
  distribution={state:1.}
  for i in (range(M-1,-1,-1) if reverse else range(M)):
   following={}
   for source,p in distribution.items():
    target=source^(1<<i);accept=min(1.,weights[target]/weights[source]);flip=proposal_probability*accept
    following[target]=following.get(target,0)+p*flip
    following[source]=following.get(source,0)+p*(1-flip)
   distribution={s:p for s,p in following.items() if p>0}
  matrix.append([distribution.get(s,0.) for s in range(1<<M)])
 return matrix
def components(matrix):
 adj={i:{j for j,p in enumerate(row) if p>0} for i,row in enumerate(matrix)}
 def reach(start):
  found={start};pending=[start]
  while pending:
   for node in adj[pending.pop()]-found:found.add(node);pending.append(node)
  return found
 reachability={i:reach(i) for i in adj};pending=set(adj);groups=[]
 while pending:
  node=min(pending);group=sorted(j for j in reachability[node] if node in reachability[j]);groups.append(group);pending-=set(group)
 return groups
def compose(first,second):
 return [[sum(first[i][k]*second[k][j] for k in range(1<<M)) for j in range(1<<M)] for i in range(1<<M)]
original=transition(1.);lazy=transition(.5)
original_cycle=compose(original,transition(1.,reverse=True));lazy_cycle=compose(lazy,transition(.5,reverse=True))
classes=components(original);lazy_classes=components(lazy)
cycle_classes=components(original_cycle);lazy_cycle_classes=components(lazy_cycle)
closed=[s for s in range(1<<M) if all(((s>>i)&1) != ((s>>(i+1))&1) for i in range(M-1))]
assert len(closed)==2 and sorted(map(len,classes))==[2,30] and len(lazy_classes)==1
assert sorted(map(len,cycle_classes))==[1,1,30] and len(lazy_cycle_classes)==1

def observables(state):
 w=weights[state];D=x*x/w;N=1+(x*x-1)/w;E=U*(D-N/2+.25)
 return E,N,D

def mean_for(states):
 z=sum(weights[s] for s in states)
 return {name:sum(weights[s]*observables(s)[i] for s in states)/z for i,name in enumerate(('Ecenter','N','D'))}
all_states=list(range(1<<M));open_states=[s for s in all_states if s not in closed]
full,opened,trapped=mean_for(all_states),mean_for(open_states),mean_for(closed)

def scale_mix(number_trapped):return {name:(4-number_trapped)*opened[name]+number_trapped*trapped[name] for name in full}
report={'status':'PASS_EXACT_DIAGNOSIS_AND_LAZY_ATOMIC_TRANSITION_PROOF','method':'Independent analytical HS determinant weights and complete scan transition enumeration; no user sampler imported','model':{'time_slices':M,'U':U,'beta':beta,'mu_centered':mu,'dtau':dtau,'lambda':lam,'physical_sites':4},'single_site_weight':'(1+exp(beta*mu+lambda*sum_s))*(1+exp(beta*mu-lambda*sum_s))','original_forward_scan_components':classes,'original_component_sizes':list(map(len,classes)),'closed_alternating_states':[{'integer':s,'temporal_fields':[1 if (s>>i)&1 else -1 for i in range(M)],'target_probability':pi[s]} for s in closed],'closed_class_total_target_probability':sum(pi[s] for s in closed),'uniform_initial_closed_class_probability':len(closed)/(1<<M),'single_site_exact':full,'single_site_open_class':opened,'single_site_closed_class':trapped,'four_site_exact':{k:4*v for k,v in full.items()},'four_site_zero_trapped':scale_mix(0),'four_site_one_trapped':scale_mix(1),'four_site_uniform_initial_class_mixture_limit':{name:4*((1-len(closed)/(1<<M))*opened[name]+len(closed)/(1<<M)*trapped[name]) for name in full},'lazy_proposal_probability':.5,'lazy_component_count':len(lazy_classes),'lazy_minimum_full_scan_self_loop':min(lazy[i][i] for i in range(1<<M)),'original_stationarity_max_error':max(abs(sum(pi[i]*original[i][j] for i in range(1<<M))-pi[j]) for j in range(1<<M)),'lazy_stationarity_max_error':max(abs(sum(pi[i]*lazy[i][j] for i in range(1<<M))-pi[j]) for j in range(1<<M)),'scope':'Proves this atomic positive-weight scan reducibility and its lazy repair. It is not a global ergodicity proof for arbitrary determinant-zero/sign-changing many-body supports.'}
report['transition_conventions']={'forward_scan':'Visits temporal fields0,1,2,3,4 exactly once','forward_reverse_cycle':'Forward scan followed by visits4,3,2,1,0; matches update_replica! measurement cycle'}
report['original_forward_reverse_cycle_component_sizes']=list(map(len,cycle_classes))
report['original_forward_reverse_cycle_alternating_self_transition']={str(s):original_cycle[s][s] for s in closed}
report['lazy_forward_reverse_cycle_component_count']=len(lazy_cycle_classes)
report['lazy_minimum_forward_reverse_cycle_self_loop']=min(lazy_cycle[i][i] for i in range(1<<M))
report['lazy_forward_reverse_cycle_stationarity_max_error']=max(abs(sum(pi[i]*lazy_cycle[i][j] for i in range(1<<M))-pi[j]) for j in range(1<<M))
(ROOT/'exact_atomic_scan_graph.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('original_forward_scan_components',)},indent=2))

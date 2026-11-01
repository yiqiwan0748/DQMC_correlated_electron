#!/usr/bin/env python3
"""Define exactly 100 external scientific benchmark comparisons. No simulations."""
from pathlib import Path
import json,csv,collections
P=Path(__file__).resolve().parent
source_quspin='https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html'
source_paper='https://doi.org/10.21468/SciPostPhys.7.2.020'
source_peschel='https://arxiv.org/abs/cond-mat/0212631'
source_hubbard='https://doi.org/10.1098/rspa.1963.0204'
fixtures=[(4,0.,'apbc'),(4,2.,'pbc'),(4,2.,'apbc'),(4,4.,'pbc'),(4,4.,'apbc'),(4,8.,'pbc'),(6,2.,'pbc'),(6,4.,'apbc')]
obs=[('energy_uncentered',[],None),('double_occupancy_total',[],None),('spin_z_correlation',[0,1],None),('spin_z_correlation',[0,2],None),('onebody_up',[0,1],None),('onebody_up',[0,2],None),('S2',[],1),('S2',[],0)]
cases=[]
for L,U,bc in fixtures:
 for ob,sites,na in obs:
  cases.append(dict(id=f'QG-{len(cases)+1:02d}',family='QuSpin_ground',L=L,t=1.,U=U,mu_centered=0.,boundary=bc,Nup=L//2,Ndown=L//2,ensemble='canonical_ground',observable=ob,sites=sites,nA=(L//2 if na==0 else na or 0),target_methods=['ED','DMRG'],oracle='QuSpin exact fixed-sector ground state; unique-ground-state gate',source_urls=[source_quspin,source_paper],absolute_tolerance_ed=1e-8,absolute_tolerance_dmrg=(1e-7 if ob.startswith('energy') else 2e-5),runtime_tier='T1',status='NOT_RUN',independence='Separate well-known package implements Hamiltonian, diagonalization and partial trace; target ED/DMRG use their own basis/MPO'))
for fidx,(U,beta) in enumerate([(2.,1.),(4.,2.)]):
 for oidx,(ob,sites,na) in enumerate(obs):
  cases.append(dict(id=f'QT-{fidx*8+oidx+1:02d}',family='QuSpin_thermal_trotter',L=4,t=1.,U=U,mu_centered=0.,boundary='pbc',Nup=2,Ndown=2,beta=beta,dtau=.1,ensemble='grand_canonical_trotter',observable=ob,sites=sites,nA=(2 if na==0 else na or 0),target_methods=['thermal_ED','DQMC'],oracle='QuSpin full grand Fock space; normalized [exp(-dtau V) exp(-dtau K)]^(beta/dtau)',source_urls=[source_quspin,source_paper,'https://arxiv.org/abs/1307.1486'],absolute_tolerance_ed=1e-8,equivalence_margin=(.05 if ob in ['S2','energy_uncentered'] else .03),runtime_tier='T3',status='NOT_RUN',independence='Independent package finite-dtau thermal oracle; no comparison to a pure ground state'))
for L,beta,bc in [(8,.5,'pbc'),(8,4.,'apbc'),(6,2.,'pbc'),(6,1.,'apbc')]:
 for ob in ['energy_uncentered','onebody_up','S2']:
  cases.append(dict(id=f'FREE-{sum(c["family"]=="analytic_free" for c in cases)+1:02d}',family='analytic_free',L=L,t=1.,U=0.,mu_centered=0.,boundary=bc,beta=beta,dtau=.1,ensemble='grand_canonical_thermal',observable=ob,sites=([0,1] if ob=='onebody_up' else []),nA=L//2 if ob=='S2' else 0,target_methods=['DQMC'],oracle='Exact Fourier modes epsilon_k=-2t cos(k), Fermi weights; restricted-correlation Gaussian purity',source_urls=[source_peschel,'https://arxiv.org/abs/0906.1663'],absolute_tolerance_ed=1e-9,runtime_tier='T0',status='NOT_RUN',independence='Closed analytic one-particle dispersion and Gaussian reduced-density formula; no Monte Carlo uncertainty at U=0'))
for U,beta,mu in [(2.,1.,0.),(4.,.5,.3)]:
 for ob in ['energy_centered','density_total','double_occupancy_total','S2']:
  cases.append(dict(id=f'ATOM-{sum(c["family"]=="analytic_atom" for c in cases)+1:02d}',family='analytic_atom',L=4,t=0.,U=U,mu_centered=mu,boundary='pbc',beta=beta,dtau=.1,ensemble='grand_canonical_thermal',observable=ob,sites=[],nA=2 if ob=='S2' else 0,target_methods=['thermal_ED','DQMC'],oracle='Four local occupation states with exact centered Boltzmann weights; independent-site factorization',source_urls=[source_hubbard],absolute_tolerance_ed=1e-10,equivalence_margin=.05 if ob=='S2' else .03,runtime_tier='T2',status='NOT_RUN',independence='Analytical diagonal four-state Hubbard atom, derived explicitly in guide; interactions commute exactly so no Trotter bias'))
assert len(cases)==100
payload=dict(schema_version=2,description='100 externally referenced benchmark comparisons, not 100 independent numerical methods. Parameter-family repeats are explicit. Internal regressions do not count toward this 100.',family_counts=dict(collections.Counter(c['family'] for c in cases)),status='PROPOSED_NOT_EXECUTED',cases=cases)
(P/'validation_100.json').write_text(json.dumps(payload,indent=2)+'\n')
def val(x):
 if isinstance(x,bool):return str(x).lower()
 if isinstance(x,(float,int)):return str(x)
 if isinstance(x,list):return '['+', '.join(val(y) for y in x)+']'
 return json.dumps(x)
with (P/'benchmark_manifest.toml').open('w') as out:
 out.write('schema_version = 2\nfamilywise_alpha = 0.01\nscalar_comparisons = 100\n\n')
 for c in cases:
  out.write('[[cases]]\n')
  for k,v in c.items():out.write(k+' = '+val(v)+'\n')
  out.write('\n')
fields=['id','family','observable','L','t','U','mu_centered','boundary','ensemble','beta','dtau','nA','sites','target_methods','oracle','source_urls','status']
with (P/'validation_100.csv').open('w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=fields,extrasaction='ignore');w.writeheader()
 for c in cases:w.writerow({k:json.dumps(v) if isinstance(v,list) else v for k,v in c.items()})
print(json.dumps({'cases':len(cases),'family_counts':payload['family_counts']}))

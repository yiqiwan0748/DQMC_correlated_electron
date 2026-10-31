"""Prospective P13 independent TeNPy core; no scientific imports or automatic run.

A future reviewed wrapper must supply the actual source/runtime-qualified
NumPy/npc/MPS/ExactDiag/FermiHubbardChain objects before calling produce.
No user Hamiltonian/vector/result is accepted as an argument.
"""
import copy
import hashlib
import inspect
from pathlib import Path
from protocol_helpers import L,CUT,U_INITIAL,U_FINAL,ARMS,coordinates,regroup_sign

def pack(array,np):
 a=np.asarray(array)
 if a.dtype.kind not in 'fciub' or not np.all(np.isfinite(a)):
  raise RuntimeError('Only finite numeric actual arrays may be exported')
 return {'shape':list(a.shape),'order':'F','real':a.real.astype(np.float64).ravel(order='F').tolist(),'imag':a.imag.astype(np.float64).ravel(order='F').tolist()}

def transpose_actual(array,names):
 labels=list(array.get_leg_labels())
 if len(labels)!=len(names)or set(labels)!=set(names):
  raise RuntimeError('Unexpected actual tensor-axis labels')
 return array.to_ndarray().transpose([labels.index(name)for name in names])

def operators(site):
 names={'Id':'Id','Nup':'Nu','Ndn':'Nd','Ntot':'Ntot','Sz':'Sz','Cdagup':'Cdu','Cdagdn':'Cdd','Cup':'Cu','Cdn':'Cd','F':'JW'}
 return {key:transpose_actual(site.get_op(name),['p','p*'])for key,name in names.items()}

def local_layout(op,np):
 occ={}
 for j in range(4):
  u,d=op['Nup'][j,j],op['Ndn'][j,j]
  if u not in(0,1)or d not in(0,1):raise RuntimeError('Nonbinary actual local occupations')
  k=(int(u.real),int(d.real))
  if k in occ:raise RuntimeError('Repeated actual occupation label')
  occ[k]=j
 if set(occ)!={(0,0),(1,0),(0,1),(1,1)}:raise RuntimeError('Incomplete local Fock basis')
 vac=np.zeros(4);vac[occ[0,0]]=1
 vectors={(0,0):vac,(1,0):op['Cdagup']@vac,(0,1):op['Cdagdn']@vac,(1,1):op['Cdagup']@op['Cdagdn']@vac}
 phases={k:v[occ[k]]for k,v in vectors.items()}
 if any(x not in(-1,1)for x in phases.values()):raise RuntimeError('Local phase not exact+/-1')
 for k,v in vectors.items():
  expected=np.zeros(4,dtype=np.complex128);expected[occ[k]]=phases[k]
  if not np.array_equal(v,expected):raise RuntimeError('Local canonical creation action failed')
 cu,cd=op['Cup'],op['Cdn'];eye=np.eye(4);zero=np.zeros((4,4))
 checks=[op['Id']==eye,op['Cdagup']==cu.conj().T,op['Cdagdn']==cd.conj().T,cu@cu.conj().T+cu.conj().T@cu==eye,cd@cd.conj().T+cd.conj().T@cd==eye,cu@cd+cd@cu==zero,cu@cd.conj().T+cd.conj().T@cu==zero,op['Nup']==cu.conj().T@cu,op['Ndn']==cd.conj().T@cd,op['Ntot']==op['Nup']+op['Ndn'],op['Sz']==(op['Nup']-op['Ndn'])/2,op['F']==eye-2*op['Ntot']+4*(op['Nup']@op['Ndn'])]
 if not all(np.all(x)for x in checks):raise RuntimeError('Actual CAR/number/JW qualification failed')
 return occ,phases

def canonical_projection(ed,layouts,np):
 selected=np.flatnonzero(ed._mask)
 if selected.shape!=(400,)or ed.full_H.shape!=(400,400):
  raise RuntimeError('Actual balanced charge projection must be exactly400')
 indices=[];phases=[];digits=[]
 for u,d in coordinates():
  occ=[((u>>j)&1,(d>>j)&1)for j in range(L)]
  rowdigits=tuple(layouts[j][0][occ[j]]for j in range(L))
  full=int(ed._pipe.map_incoming_flat(rowdigits))
  pos=int(np.searchsorted(selected,full))
  if pos>=400 or int(selected[pos])!=full:
   raise RuntimeError('Canonical state not in actual selected charge sector')
  phase=regroup_sign(u,d)
  for j in range(L):phase*=layouts[j][1][occ[j]]
  indices.append(pos);phases.append(phase);digits.append(list(rowdigits))
 if sorted(indices)!=list(range(400)):raise RuntimeError('Canonical projection is not a complete permutation')
 return np.asarray(indices,dtype=np.int64),np.asarray(phases,dtype=np.complex128),digits,selected

def export_mpo(actual,np):
 left,right=actual.get_IdL(0),actual.get_IdR(L-1)
 if left is None or right is None:raise RuntimeError('Missing actual finite-MPO identity boundary')
 return {'axes':['left','right','bra','ket'],'sites':[pack(transpose_actual(actual.get_W(j),['wL','wR','p','p*']),np)for j in range(L)],'left_boundary_index':int(left),'right_boundary_index':int(right)}

def restore_for_registered_mps(ed,vector,MPS,npc,cutoff):
 # Exact source-derived LegPipe restoration from ExactDiag.full_to_mps,
 # followed by explicit MPS.from_full cutoff instead of an invented chi_max.
 full=npc.zeros([ed._pipe],vector.dtype,vector.qtotal)
 full[ed._mask]=vector
 full.iset_leg_labels(['('+'.'.join(ed._labels_p)+')'])
 full=full.split_legs([0])
 return MPS.from_full(ed._sites,full,form='B',cutoff=cutoff,normalize=True,bc='finite',unit_cell_width=ed.model.lat.mps_unit_cell_width)

def export_after_entropy_mps(psi,np):
 scalar=float(psi.norm);forms=[tuple(x)for x in psi.form]
 left=np.asarray(psi.get_SL(0));right=np.asarray(psi.get_SR(L-1))
 if scalar!=1.0 or forms!=[(0.0,1.0)]*L or not np.array_equal(left,np.ones(1))or not np.array_equal(right,np.ones(1)):
  raise RuntimeError('Post-entropy MPS not normalized all-B/unit-boundary form')
 theta=transpose_actual(psi.get_theta(0,n=L),['vL']+[f'p{j}'for j in range(L)]+['vR'])
 if theta.shape[0]!=1 or theta.shape[-1]!=1:raise RuntimeError('Nontrivial MPS endpoint')
 return {'axes':['left','right','physical'],'sites':[pack(transpose_actual(psi.get_B(j,form='B'),['vL','vR','p']),np)for j in range(L)],'complete_site_amplitudes':pack((scalar*theta[0,...,0]).ravel(order='F'),np),'external_norm_scalar':scalar,'actual_canonical_forms':[list(x)for x in forms],'actual_left_schmidt':left.tolist(),'actual_right_schmidt':right.tolist(),'scalar_placement':'UNSCALED actual B tensors; actual psi.norm multiplied at END of complete contraction','export_stage':'AFTER_ACTUAL_TENPY_ENTROPY_CALL'}

def callable_witness(name,function):
 actual=getattr(function,'__func__',function)
 module=inspect.getmodule(actual)
 try:source=inspect.getsourcefile(actual)
 except (TypeError,OSError):source=None
 try:line=inspect.getsourcelines(actual)[1] if source else None
 except (TypeError,OSError):line=None
 module_path=getattr(module,'__file__',None)if module is not None else None
 path=source or module_path
 result={'name':name,'callable':callable(actual),'callable_type':type(actual).__module__+'.'+type(actual).__qualname__,'module':getattr(actual,'__module__',None),'qualname':getattr(actual,'__qualname__',getattr(actual,'__name__',None)),'python_source_available':source is not None,'source_line':line}
 if path is not None:
  real=Path(path).resolve(strict=True);data=real.read_bytes()
  result.update(file=str(real),bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
 else:result['file_resolution_scope']='No standalone Python source path; qualified array-extension/provider records remain mandatory.'
 return result

def api_witnesses(np,npc,MPS,ExactDiag,FermiHubbardChain):
 entries=[('FermiHubbardChain.init_terms',FermiHubbardChain.init_terms),('ExactDiag.__init__',ExactDiag.__init__),('ExactDiag.build_full_H_from_mpo',ExactDiag.build_full_H_from_mpo),('ExactDiag.full_diagonalization',ExactDiag.full_diagonalization),('ExactDiag.groundstate',ExactDiag.groundstate),('ExactDiag.exp_H',ExactDiag.exp_H),('ExactDiag.mps_to_full',ExactDiag.mps_to_full),('MPS.from_full',MPS.from_full),('MPS.entanglement_entropy',MPS.entanglement_entropy),('MPS.get_theta',MPS.get_theta),('MPS.get_B',MPS.get_B),('LegPipe.map_incoming_flat',npc.LegPipe.map_incoming_flat),('npc.zeros',npc.zeros),('npc.tensordot',npc.tensordot),('numpy.exp',np.exp)]
 result=[callable_witness(name,function)for name,function in entries]
 extension=np._core._multiarray_umath
 path=Path(extension.__file__).resolve(strict=True);data=path.read_bytes()
 result.append({'name':'actual_NumPy_array_extension','module':extension.__name__,'file':str(path),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'scope':'Compiled array provider identity; no invented Python source for ufunc.'})
 return result

def qualified_stage(runtime_guard,stage,np,npc,MPS,ExactDiag,FermiHubbardChain):
 witness=api_witnesses(np,npc,MPS,ExactDiag,FermiHubbardChain)
 receipt=runtime_guard(stage,witness)
 if not isinstance(receipt,dict)or receipt.get('status')!='PASS':
  raise RuntimeError('Fresh P13 actual imported API/provider/runtime qualification failed')
 return receipt

def raw_pipe_mapping(pipe):
 return {'nlegs':int(pipe.nlegs),'ind_len':int(pipe.ind_len),'subshape':[int(x)for x in pipe.subshape],'qconj':int(pipe.qconj),'slices':pipe.slices.tolist(),'q_map':pipe.q_map.tolist(),'qflat':pipe.to_qflat().tolist(),'chinfo_mod':pipe.chinfo.mod.tolist(),'chinfo_names':list(pipe.chinfo.names),'incoming_legs':[{'ind_len':int(leg.ind_len),'qconj':int(leg.qconj),'slices':leg.slices.tolist(),'qflat':leg.to_qflat().tolist()}for leg in pipe.legs]}

def produce(np,npc,MPS,ExactDiag,FermiHubbardChain,*,pilot=False,on_record,runtime_guard):
 runtime_before=qualified_stage(runtime_guard,'before',np,npc,MPS,ExactDiag,FermiHubbardChain)
 on_record('runtime_before',copy.deepcopy(runtime_before))
 models=[];objects=[];mappings=[];all_local=[]
 for U in(U_INITIAL,U_FINAL):
  model=FermiHubbardChain({'L':L,'t':1.0,'U':U,'mu':0.0,'bc_MPS':'finite','bc_x':'open','cons_N':'N','cons_Sz':'Sz'})
  sites=model.lat.mps_sites();local=[operators(s)for s in sites];layouts=[local_layout(x,np)for x in local]
  # Derive the actual total charge from a balanced product occupation pattern.
  charge=np.zeros(sites[0].leg.chinfo.qnumber,dtype=np.int64)
  for j,s in enumerate(sites):charge+=s.leg.to_qflat()[layouts[j][0][(1,0)if j%2==0 else(0,1)]]
  charge=sites[0].leg.chinfo.make_valid(charge)
  ed=ExactDiag(model,charge_sector=charge,max_size=16777216)
  ed.build_full_H_from_mpo()
  if ed.full_H is None or ed.full_H.shape!=(400,400):raise RuntimeError('No actual complete400x400 charge Hamiltonian')
  indices,phases,digits,selected=canonical_projection(ed,layouts,np)
  matrix=ed.full_H.to_ndarray();canonical=phases.conj()[:,None]*matrix[np.ix_(indices,indices)]*phases[None,:]
  if np.any(canonical.imag!=0):raise RuntimeError('Unexpected complex physical Hamiltonian')
  ed.full_diagonalization()
  values=np.asarray(ed.E);vectors=ed.V.to_ndarray()
  if values.shape!=(400,)or vectors.shape!=(400,400)or np.any(np.diff(values)<0):raise RuntimeError('Complete sorted400eigensystem required')
  cv=phases.conj()[:,None]*vectors[indices,:]
  if np.any(cv.imag!=0):raise RuntimeError('Unexpected complex real-Hubbard eigensystem')
  models.append({'U':U,'represented_mpo':export_mpo(model.H_MPO,np),'Hamiltonian':pack(canonical.real,np),'eigenvalues':values.tolist(),'eigenvectors':pack(cv.real,np),'actual_charge_sector':charge.tolist(),'actual_charge_names':list(sites[0].leg.chinfo.names),'actual_charge_mask':selected.tolist(),'actual_canonical_to_selected_permutation':indices.tolist(),'actual_canonical_local_digits':digits,'actual_canonical_phases':pack(phases,np),'actual_pipe_charge_table':ed._pipe.to_qflat().tolist(),'actual_pipe_mapping':raw_pipe_mapping(ed._pipe),'eigensolver':'Actual TeNPy ExactDiag.full_diagonalization on projected400x400full_H'})
  on_record('model',copy.deepcopy(models[-1]))
  objects.append(ed);mappings.append((indices,phases));all_local.append(local)
 initial_ed,final_ed=objects
 # A LegCharge equality alone is insufficient inside the400-dimensional block:
 # all states have the same total charge. Require the physical maps themselves.
 for key in ['actual_charge_sector','actual_charge_names','actual_charge_mask','actual_canonical_to_selected_permutation','actual_canonical_local_digits','actual_canonical_phases','actual_pipe_charge_table','actual_pipe_mapping']:
  if models[0][key]!=models[1][key]:raise RuntimeError('Initial/final physical basis mapping differs: '+key)
 for initial_local,final_local in zip(all_local[0],all_local[1]):
  if set(initial_local)!=set(final_local)or not all(np.array_equal(initial_local[k],final_local[k])for k in initial_local):
   raise RuntimeError('Initial/final actual local physical operators differ')
 _,initial=initial_ed.groundstate()
 # Exact equality of incoming physical charge legs is required before action.
 initial.legs[0].test_equal(final_ed.full_H.legs[0])
 initial_array=initial.to_ndarray();initial_norm=float(np.linalg.norm(initial_array))
 if not np.isfinite(initial_norm)or initial_norm<=0:raise RuntimeError('Invalid initial norm')
 initial=initial/initial_norm
 initial_indices,initial_phases=mappings[0]
 own_initial=initial_phases.conj()*initial.to_ndarray()[initial_indices]
 indices,phases=mappings[1]
 samples=[]
 for k in range(5 if pilot else 33):
  time=k/8
  actual_U=final_ed.exp_H(time)
  if actual_U.shape!=(400,400)or list(actual_U.get_leg_labels())!=['ps','ps*']:
   raise RuntimeError('Actual ExactDiag.exp_H shape/axes differ')
  raw=npc.tensordot(actual_U,initial,axes=[['ps*'],['ps']])
  raw_array=raw.to_ndarray();raw_norm=float(np.linalg.norm(raw_array))
  if not np.isfinite(raw_norm)or raw_norm<=0:raise RuntimeError('Invalid evolved raw norm')
  normalized=raw/raw_norm
  for arm,settings in ARMS.items():
   budget=settings['tenpy_budget'];cutoff=settings['tenpy_cutoff']
   psi=restore_for_registered_mps(final_ed,normalized.copy(),MPS,npc,cutoff)
   value=float(psi.entanglement_entropy(n=2,bonds=[CUT])[0])
   maxbond=int(max(psi.chi))
   if maxbond>64 or maxbond>=budget:raise RuntimeError('Registered external bond budget exhausted')
   state=export_after_entropy_mps(psi,np)
   schmidt=np.asarray(psi.get_SL(CUT))
   samples.append({'time_numerator':k,'time_denominator':8,'time':time,'arm':arm,'external_bond_budget':budget,'from_full_cutoff':cutoff,'actual_raw_evolved_state':pack(phases.conj()*raw_array[indices],np),'raw_norm':raw_norm,'actual_normalized_MPS_input_sector_state':pack(phases.conj()*normalized.to_ndarray()[indices],np),'library_S2':value,'library_schmidt_probability_sum':float(np.sum(schmidt**2)),'observed_max_bond':maxbond,'post_call_schmidt_values_new_observer':schmidt.tolist(),'state':state})
   on_record('sample',copy.deepcopy(samples[-1]))
 runtime_after=qualified_stage(runtime_guard,'after',np,npc,MPS,ExactDiag,FermiHubbardChain)
 on_record('runtime_after',copy.deepcopy(runtime_after))
 return {'schema':'P13_TENPY_EXACTDIAG_ENTROPY_PRODUCER_CORE_V2','runtime_before':runtime_before,'runtime_after':runtime_after,'status':'PRODUCER_COMPLETE_NOT_A_SCIENTIFIC_PASS','scientific_credit':0,'catalogue_id':'P13','pilot_only':pilot,'coordinates':[list(x)for x in coordinates()],'local_operators_initial_and_final':[[{key:pack(value,np)for key,value in table.items()}for table in local]for local in all_local],'models':models,'actual_initial_state':pack(own_initial,np),'original_initial_eigenvector_norm':initial_norm,'initial_final_physical_basis_mapping_equal':True,'samples':samples,'time_discretization_error':'ZERO_TROTTER_DIRECT_COMPLETE_SPECTRAL_EXPONENTIAL;floatingeigensystem/phase/arithmeticerrorseparate','attribution':'Independent actual TeNPy FermiHubbardChain/ExactDiag.exp_H/MPS.from_full/Schmidtentropy; NEW registration,basis,configuration,complete-export and certificate adapters. No target matrices/vectors/results are received.'}

# P13 prospective numerical core. NO CLI or automatic execution.
# A future separately reviewed source/runtime/admission wrapper must qualify
# actual packages/providers before loading and invoking this module.
module P13NativeProducerCore
using LinearAlgebra, SHA
using ITensors, ITensorMPS
module ActualUserP13
include(joinpath(@__DIR__,"source_snapshot","compute_s2_itensor_dmrg.jl"))
end
const L=6
const CUT=3
const UI=8.0
const UF=4.0
const ARMS=(("coarse",128,1e-24),("primary",256,1e-28))
pack(x)=Dict("shape"=>collect(size(x)),"order"=>"F","real"=>vec(Float64.(real.(x))),"imag"=>vec(Float64.(imag.(x))))
coords()=[(u,d) for u in 0:63 for d in 0:63 if count_ones(u)==3 && count_ones(d)==3]
function regroup_sign(u,d)
 # Count down_j/up_k crossings directly. An empty suffix contributes zero.
 inversions=0
 for j in 0:5, k in (j+1):5
  inversions+=((d>>j)&1)*((u>>k)&1)
 end
 isodd(inversions) ? -1.0 : 1.0
end
function local_tables(sites)
 [Dict(name=>Array(op(name,s),prime(s),dag(s)) for name in ("Id","Nup","Ndn","Ntot","Sz","Cdagup","Cdagdn","Cup","Cdn","F")) for s in sites]
end
function actual_layout(tables)
 out=Any[]
 for op in tables
  occ=Dict{Tuple{Int,Int},Int}()
  for j in 1:4
   u,d=op["Nup"][j,j],op["Ndn"][j,j]
   u in (0,1)&&d in (0,1)||error("Nonbinary actual local occupation")
   k=(Int(real(u)),Int(real(d)));haskey(occ,k)&&error("Repeated occupation label");occ[k]=j
  end
  Set(keys(occ))==Set(((0,0),(1,0),(0,1),(1,1)))||error("Incomplete occupation map")
  vac=zeros(ComplexF64,4);vac[occ[(0,0)]]=1
  vectors=Dict((0,0)=>vac,(1,0)=>op["Cdagup"]*vac,(0,1)=>op["Cdagdn"]*vac,(1,1)=>op["Cdagup"]*op["Cdagdn"]*vac)
  phases=Dict(k=>v[occ[k]] for (k,v) in vectors)
  all(p in(-1,1) for p in values(phases))||error("Actual local phase must be exact+/-1")
  for (k,v) in vectors
   ideal=zeros(ComplexF64,4);ideal[occ[k]]=phases[k];v==ideal||error("Local state not exact canonical Fock state")
  end
  id=Matrix{Float64}(I,4,4);z=zeros(4,4);cu=op["Cup"];cd=op["Cdn"]
  checks=(op["Id"]==id,op["Cdagup"]==cu',op["Cdagdn"]==cd',cu*cu'+cu'*cu==id,cd*cd'+cd'*cd==id,cu*cd+cd*cu==z,cu*cd'+cd'*cu==z,op["Nup"]==cu'*cu,op["Ndn"]==cd'*cd,op["Ntot"]==op["Nup"]+op["Ndn"],op["Sz"]==(op["Nup"]-op["Ndn"])/2,op["F"]==id-2*op["Ntot"]+4*op["Nup"]*op["Ndn"])
  all(checks)||error("Actual local CAR/number/JW identity failed")
  push!(out,(occ,phases))
 end
 out
end
function basis_map(layouts)
 digits=zeros(Int,400,L);phases=ones(ComplexF64,400)
 for (row,(u,d)) in enumerate(coords())
  phases[row]=regroup_sign(u,d)
  for j in 1:L
   key=((u>>(j-1))&1,(d>>(j-1))&1)
   digits[row,j]=layouts[j][1][key];phases[row]*=layouts[j][2][key]
  end
 end
 digits,phases
end
function mpo_arrays(H,sites)
 out=Any[]
 for j in 1:L
  left=j>1 ? (linkind(H,j-1),) : ();right=j<L ? (linkind(H,j),) : ()
  a=Array(H[j],left...,right...,prime(sites[j]),dag(sites[j]))
  dl=j>1 ? dim(only(left)) : 1;dr=j<L ? dim(only(right)) : 1
  push!(out,reshape(a,dl,dr,4,4))
 end
 out
end
function sector_from_actual_mpo(tensors,digits,phases)
 # Only160000balanced-sector entries; never form a4096x4096 matrix.
 matrix=zeros(Float64,400,400)
 for row in 1:400,col in 1:400
  x=ComplexF64[1]
  for j in 1:L
   a=tensors[j];size(a,1)==length(x)||error("Actual MPO link mismatch")
   y=zeros(ComplexF64,size(a,2))
   for r in eachindex(y),l in eachindex(x)
    y[r]+=x[l]*a[l,r,digits[row,j],digits[col,j]]
   end
   x=y
  end
  length(x)==1||error("Nontrivial actual MPO endpoint")
  value=conj(phases[row])*x[1]*phases[col]
  iszero(imag(value))||error("Unexpected complex real-Hubbard matrix entry")
  matrix[row,col]=real(value)
 end
 matrix
end
function mps_state_export(psi,sites)
 tensors=Any[]
 for j in 1:L
  left=j>1 ? (linkind(psi,j-1),) : ();right=j<L ? (linkind(psi,j),) : ()
  a=Array(psi[j],left...,right...,sites[j]);dl=j>1 ? dim(only(left)) : 1;dr=j<L ? dim(only(right)) : 1
  push!(tensors,pack(reshape(a,dl,dr,4)))
 end
 complete=psi[1]
 for j in 2:L;complete*=psi[j];end
 Dict("axes"=>["left","right","physical"],"sites"=>tensors,"complete_site_amplitudes"=>pack(vec(Array(complete,sites...))),"external_norm_scalar"=>1.0,"export_stage"=>"AFTER_ACTUAL_USER_ENTROPY_CALL")
end
function full_site_tensor(vector,digits,phases)
 a=zeros(ComplexF64,ntuple(_->4,L))
 for row in 1:400;a[Tuple(digits[row,:])...]=phases[row]*vector[row];end
 a
end
function method_witness(name,f,types)
 m=which(f,types);raw=String(m.file)
 record=Dict{String,Any}("name"=>name,"module"=>string(m.module),"signature"=>string(m.sig),"raw_file"=>raw,"line"=>Int(m.line),"current_file_readable"=>isfile(raw))
 if isfile(raw)
  path=realpath(raw);record["file"]=path;record["bytes"]=filesize(path);record["sha256"]=bytes2hex(sha256(read(path)))
 else
  record["source_resolution_scope"]="Raw compiled-image source location; fresh runtime guard must bind the Julia image and independently resolve against its source closure."
 end
 record
end
function api_witnesses()
 [method_witness("actual_user_hubbard_mpo",ActualUserP13.hubbard_mpo,(Vector{Index{Int}},)),
  method_witness("actual_user_renyi2_bond_entropy",ActualUserP13.renyi2_bond_entropy,(MPS,Int)),
  method_witness("actual_complex_ITensor_constructor",ITensor,Tuple{Array{ComplexF64,6},Vararg{Index{Int},6}}),
  method_witness("actual_dense_MPS_constructor",MPS,(ITensor,Vector{Index{Int}})),
  method_witness("actual_complete_real_eigen",eigen,(Symmetric{Float64,Matrix{Float64}},)),
  method_witness("actual_complex_phase_cis",cis,(Float64,)),
  method_witness("actual_MPS_normalize",normalize!,(MPS,))]
end
function qualified_stage(runtime_guard,stage)
 witness=api_witnesses();receipt=runtime_guard(stage,witness)
 receipt isa AbstractDict && get(receipt,"status",nothing)=="PASS"||error("Fresh P13 native runtime qualification failed")
 receipt
end
function produce(;pilot::Bool=false,on_record::Function,runtime_guard::Function)
 # This is a numerical core, not an admission boundary. Root wrapper is required.
 string(VERSION)=="1.12.1"||error("Unreviewed Julia version")
 string(pkgversion(ITensors))=="0.9.30"||error("Unreviewed ITensors version")
 string(pkgversion(ITensorMPS))=="0.4.1"||error("Unreviewed ITensorMPS version")
 BLAS.set_num_threads(1)
 runtime_before=qualified_stage(runtime_guard,"before")
 on_record("runtime_before",deepcopy(runtime_before))
 sites=siteinds("Electron",L;conserve_qns=false)
 tables=local_tables(sites);layouts=actual_layout(tables);digits,phases=basis_map(layouts)
 models=Any[];spectra=Any[]
 for U in(UI,UF)
  H=ActualUserP13.hubbard_mpo(sites;L=L,t=1.0,U=U,mu=0.0,boundary="obc")
  tensors=mpo_arrays(H,sites);matrix=sector_from_actual_mpo(tensors,digits,phases)
  d=eigen(Symmetric(matrix,:U));issorted(d.values)||error("Complete eigensystem not sorted")
  push!(spectra,d)
  push!(models,Dict("U"=>U,"represented_mpo"=>Dict("axes"=>["left","right","bra","ket"],"sites"=>[pack(x) for x in tensors]),"Hamiltonian"=>pack(matrix),"eigenvalues"=>d.values,"eigenvectors"=>pack(d.vectors),"eigensolver"=>"Julia eigen(Symmetric(actual_complete400matrix,:U)); raw full matrix retained"))
  on_record("model",deepcopy(models[end]))
 end
 initial=ComplexF64.(spectra[1].vectors[:,1]);initial_norm=norm(initial);initial/=initial_norm
 final=spectra[2];coefficient=final.vectors'*initial;samples=Any[]
 for k in 0:(pilot ? 4 : 32)
  time=k/8;raw=final.vectors*(cis.(-time.*final.values).*coefficient)
  raw_norm=norm(raw);isfinite(raw_norm)&&raw_norm>0||error("Invalid evolved raw norm")
  normalized=raw/raw_norm
  a=full_site_tensor(normalized,digits,phases)
  for (arm,budget,cutoff) in ARMS
   psi=MPS(ITensor(a,sites...),sites;cutoff=cutoff,maxdim=budget,orthocenter=CUT)
   pre_entropy_MPS_norm=norm(psi)
   isfinite(pre_entropy_MPS_norm)&&pre_entropy_MPS_norm>0||error("Invalid MPS conversion norm")
   normalize!(psi)
   normalized_MPS_norm=norm(psi)
   S2,normcheck=ActualUserP13.renyi2_bond_entropy(psi,CUT)
   maxbond=maxlinkdim(psi);maxbond<=64&&maxbond<budget||error("Predeclared bond budget exhausted")
   # User routine has already orthogonalized psi at the registered cut.
   left=(linkind(psi,CUT-1),sites[CUT]);_,S,_=svd(psi[CUT],left)
   schmidt=[real(S[n,n]) for n in 1:dim(S,1)]
   push!(samples,Dict("time_numerator"=>k,"time_denominator"=>8,"time"=>time,"arm"=>arm,"maxdim"=>budget,"cutoff"=>cutoff,"actual_raw_evolved_state"=>pack(raw),"raw_norm"=>raw_norm,"actual_normalized_MPS_input_sector_state"=>pack(normalized),"pre_entropy_MPS_norm_before_normalize"=>pre_entropy_MPS_norm,"pre_entropy_MPS_norm_after_normalize"=>normalized_MPS_norm,"actual_normalize_call"=>"LinearAlgebra.normalize!(psi) before actual user entropy","library_S2"=>S2,"library_normcheck"=>normcheck,"observed_max_bond"=>maxbond,"post_call_schmidt_values_new_observer"=>schmidt,"state"=>mps_state_export(psi,sites)))
   on_record("sample",deepcopy(samples[end]))
  end
 end
 runtime_after=qualified_stage(runtime_guard,"after")
 on_record("runtime_after",deepcopy(runtime_after))
 Dict("schema"=>"P13_NATIVE_MPO_ENTROPY_DENSE_ADAPTER_CORE_V2","runtime_before"=>runtime_before,"runtime_after"=>runtime_after,"status"=>"PRODUCER_COMPLETE_NOT_A_SCIENTIFIC_PASS","scientific_credit"=>0,"catalogue_id"=>"P13","pilot_only"=>pilot,"coordinates"=>[collect(x) for x in coords()],"local_operators"=>[Dict(k=>pack(v) for (k,v) in x) for x in tables],"models"=>models,"actual_initial_state"=>pack(initial),"original_initial_eigenvector_norm"=>initial_norm,"samples"=>samples,"attribution"=>"Unchanged actual user OBC Hubbard MPO and renyi2_bond_entropy; ITensors/ITensorMPS tensor/SVD routines; NEW balanced-matrix,dense-evolution,complex-state/MPS and certificate adapters. No original native realtime evolver or DMRG optimizer is claimed.")
end
end

#!/usr/bin/env julia
# Actual user Hamiltonians/thermal/RDM APIs plus native DMRG and explicit new adapters.
using LinearAlgebra, Random, TOML, SHA
using ITensors, ITensorMPS
BLAS.set_num_threads(1)
length(ARGS)%2==0 || error("Expected option/value pairs")
const OPT=Dict(ARGS[i]=>ARGS[i+1] for i in 1:2:length(ARGS))
const CONFIG=abspath(get(OPT,"--manifest",joinpath(@__DIR__,"manifest.toml")))
const CFG=TOML.parsefile(CONFIG)
const OUTPUT=abspath(OPT["--output"])
ispath(OUTPUT) && error("Use a new output")
const SOURCE_ROOT=get(ENV,"RENYI_ROOT",normpath(joinpath(@__DIR__,"../../validated_package/HubbardRenyiDQMC")))
sha(p)=bytes2hex(sha256(read(p)))
for (rel,key) in (("analysis/compute_s2_ed.jl","ed_source_sha256"),("analysis/audit_u4_l4_thermal_ed.jl","thermal_source_sha256"),
                  ("analysis/compute_s2_itensor_dmrg.jl","dmrg_source_sha256"),("Project.toml","project_sha256"),("Manifest.toml","environment_manifest_sha256"))
    sha(joinpath(SOURCE_ROOT,rel))==CFG[key] || error("Unpinned source/environment: $key")
end
string(VERSION)==CFG["julia_version"] || error("Julia version")
string(pkgversion(ITensors))==CFG["itensors_version"] || error("ITensors version")
string(pkgversion(ITensorMPS))==CFG["itensormps_version"] || error("ITensorMPS version")
include(joinpath(SOURCE_ROOT,"analysis/audit_u4_l4_thermal_ed.jl"))
const A=U4L4ThermalEDAudit
const E=A.GroundStateED
module UserDMRG
include(joinpath(Main.SOURCE_ROOT,"analysis/compute_s2_itensor_dmrg.jl"))
end

function materialize(space)
    n=E.dimension(space);H=zeros(n,n);unit=zeros(n);col=zeros(n)
    for j in 1:n
        fill!(unit,0.);unit[j]=1.;E.apply_hamiltonian!(col,unit,space);H[:,j]=col
    end
    H
end
function site_map(space,L)
    indices=Int[];signs=Float64[]
    for u in space.up_masks,d in space.dn_masks
        row=1;crossings=0
        for j in 0:L-1
            row+=(Int((u>>j)&1)+2Int((d>>j)&1))*4^j
            crossings+=Int((d>>j)&1)*count_ones(u>>(j+1))
        end
        push!(indices,row);push!(signs,isodd(crossings) ? -1. : 1.)
    end
    indices,signs
end
function mpo_sector(H,sites,indices,signs)
    tensor=H[1]
    for j in 2:length(H);tensor*=H[j];end
    dense=reshape(Array(tensor,prime.(sites)...,dag.(sites)...),4^length(sites),:)
    mapped=(signs*signs').*dense[indices,indices]
    real.(mapped),maximum(abs,imag.(mapped))
end
function mps_row(label,kind,psi,H,sites,indices,signs;extra=Dict())
    tensor=psi[1]
    for j in 2:length(psi);tensor*=psi[j];end
    amps=vec(Array(tensor,sites...));coeff=signs.*amps[indices]
    n=real(inner(psi,psi));energy=real(inner(psi',H,psi))/n
    E2=real(inner(H,psi,H,psi))/n;variance=E2-energy^2
    merge(Dict{String,Any}("label"=>label,"kind"=>kind,"energy"=>energy,"energy_second_moment"=>E2,"variance"=>variance,
        "norm"=>n,"Nup"=>sum(real.(expect(psi,"Nup"))),"Ndown"=>sum(real.(expect(psi,"Ndn"))),
        "coeff_real"=>real.(coeff),"coeff_imag"=>imag.(coeff),
        "outside_sector_weight"=>real(sum(abs2,amps)-sum(abs2,coeff)),"max_bond"=>maxlinkdim(psi)),extra)
end
function annihilation(L,spin,j)
    n=1<<L;c=zeros(n^2,n^2);bit=1<<j
    for u in 0:n-1,d in 0:n-1
        occupied=spin==:up ? (u&bit)!=0 : (d&bit)!=0
        occupied || continue
        parity=spin==:up ? count_ones(u&(bit-1)) : count_ones(u)+count_ones(d&(bit-1))
        unew=spin==:up ? xor(u,bit) : u;dnew=spin==:dn ? xor(d,bit) : d
        c[unew*n+dnew+1,u*n+d+1]=isodd(parity) ? -1. : 1.
    end
    c
end
function compute()
    s=CFG["X09"];L=s["L"]
    space=E.build_space(L,s["t"],s["U"],0.,s["Nup"],s["Ndown"];boundary=s["boundary"])
    HH=materialize(space);dec=eigen(Symmetric(HH));indices,signs=site_map(space,L)
    sites=siteinds("Electron",L;conserve_qns=true)
    H=UserDMRG.hubbard_mpo(sites;L,t=s["t"],U=s["U"],mu=0.,boundary=s["boundary"])
    rows=Dict{String,Any}[]
    for (name,seed) in zip(s["starts"],s["seeds"])
        Random.seed!(seed)
        state=name=="charge_density_wave" ? ["UpDn","Emp","UpDn","Emp"] : name=="spin_domains" ? ["Up","Up","Dn","Dn"] : ["Up","Dn","Up","Dn"]
        name=="shuffled" && shuffle!(state)
        initial=random_mps(sites,state;linkdims=s["initial_link_dimension"])
        obs=DMRGObserver()
        reported,psi=dmrg(H,initial;nsweeps=length(s["maxdim_schedule"]),maxdim=s["maxdim_schedule"],cutoff=s["cutoff"],
            eigsolve_tol=s["eigsolve_tol"],noise=s["noise"],observer=obs,outputlevel=1)
        push!(rows,mps_row(name,"native_dmrg",psi,H,sites,indices,signs;extra=Dict("seed"=>seed,
            "reported_energy"=>reported,"truncation_errors"=>truncerrors(obs))))
    end
    dense_sites=siteinds("Electron",L;conserve_qns=false)
    Hd=UserDMRG.hubbard_mpo(dense_sites;L,t=s["t"],U=s["U"],mu=0.,boundary=s["boundary"])
    excited=findfirst(>(dec.values[1]+CFG["gates"]["ground_cluster_abs"]),dec.values)
    for (label,k) in zip(s["challenge_states"],(excited,length(dec.values)))
        amps=zeros(4^L);amps[indices]=signs.*dec.vectors[:,k]
        tensor=ITensor(reshape(amps,ntuple(_->4,L)),dense_sites...)
        psi=MPS(tensor,dense_sites;cutoff=0.,maxdim=4^L)
        push!(rows,mps_row(label,"manufactured_excited_mps",psi,Hd,dense_sites,indices,signs;
            extra=Dict("source_eigenvalue"=>dec.values[k],"source_eigen_residual"=>norm(HH*dec.vectors[:,k]-dec.values[k]*dec.vectors[:,k]))))
    end
    mpo_qn,mpo_qn_imag=mpo_sector(H,sites,indices,signs)
    mpo_dense,mpo_dense_imag=mpo_sector(Hd,dense_sites,indices,signs)
    X09=Dict("coordinates"=>[[Int(u),Int(d)] for u in space.up_masks for d in space.dn_masks],
        "Hamiltonian"=>[collect(HH[i,:]) for i in axes(HH,1)],
        "MPO_qn_Hamiltonian"=>[collect(mpo_qn[i,:]) for i in axes(HH,1)],
        "MPO_dense_Hamiltonian"=>[collect(mpo_dense[i,:]) for i in axes(HH,1)],
        "MPO_qn_max_imaginary"=>mpo_qn_imag,"MPO_dense_max_imaginary"=>mpo_dense_imag,
        "source_eigen_residual"=>norm(HH*dec.vectors-dec.vectors*Diagonal(dec.values)),"rows"=>rows)
    # X12 uses actual native Gibbs density and prefix/full RDM routines.
    s=CFG["X12"];L=s["L"];n=1<<L
    K,V,masks=A.full_fock_matrices(L,s["t"],s["U"],0.;boundary=s["boundary"])
    Nu=[Matrix(Diagonal([Float64((u>>j)&1) for u in masks for d in masks])) for j in 0:L-1]
    Nd=[Matrix(Diagonal([Float64((d>>j)&1) for u in masks for d in masks])) for j in 0:L-1]
    N=sum(Nu+Nd);Hphys=Matrix(K+Diagonal(V))+s["U"]/2*N-s["U"]*L/4*I
    ids=[Int(u)*n+Int(d)+1 for u in masks for d in masks if count_ones(u)+count_ones(d)==s["N_total"]]
    hc=Hphys[ids,ids];d=eigen(Symmetric(hc));keep=findall(x->x-d.values[1]<=CFG["gates"]["ground_cluster_abs"],d.values)
    vg=zeros(n^2,length(keep));vg[ids,:]=d.vectors[:,keep];mixed=vg*vg'/length(keep)
    pure_space=E.build_space(L,s["t"],s["U"],0.,s["pure_Nup"],s["pure_Ndown"];boundary=s["boundary"])
    pure_E,pu_sector,it,res=E.ground_state_lanczos(pure_space;maxiter=100,tol=1e-12)
    pureids=[Int(u)*n+Int(d)+1 for u in pure_space.up_masks for d in pure_space.dn_masks]
    pu=zeros(n^2);pu[pureids]=pu_sector
    local_minus=[annihilation(L,:dn,j)'*annihilation(L,:up,j) for j in 0:L-1]
    Sminus=sum(local_minus)
    lowering_norm=norm(Sminus*pu);pd=Sminus*pu/lowering_norm
    theta=pi/s["rotation_theta_pi_denominator"];phi=pi/s["rotation_phi_pi_denominator"]
    rot=[cos(theta) -cis(-phi)*sin(theta);cis(phi)*sin(theta) cos(theta)]
    vr=vg*rot;mixrot=vr*vr'/2;pr=cos(theta)*pu+cis(phi)*sin(theta)*pd
    sz=[(Nu[j]-Nd[j])/2 for j in 1:L]
    sx=[(m+m')/2 for m in local_minus];sy=[(m'-m)/(2im) for m in local_minus]
    function measure(label,rho;beta=-1.,epsilon=-1.)
        rdms=Dict{String,Any}[]
        for nA in s["subsystem_sizes"]
            S2,purity,rhoA=A.renyi2(rho,masks,L,nA)
            push!(rdms,Dict("nA"=>nA,"spectrum"=>reverse(eigvals(Hermitian(rhoA))),"trace"=>real(tr(rhoA)),"purity"=>purity,"S2"=>S2))
        end
        means=[real(tr(rho*z)) for z in sz];zz=real(tr(rho*sz[1]*sz[2]))
        spin_xyz=[[real(tr(rho*x)),real(tr(rho*y)),real(tr(rho*z))] for (x,y,z) in zip(sx,sy,sz)]
        Dict("label"=>label,"beta"=>beta,"excited_weight"=>epsilon,"energy"=>real(tr(rho*Hphys)),"trace"=>real(tr(rho)),
            "Nmean"=>real(tr(rho*N)),"spin_z_means"=>means,"spin_xyz_means"=>spin_xyz,"Czz"=>zz,"connected_Czz"=>zz-prod(means),
            "distance_to_equal_ground"=>sum(abs,eigvals(Hermitian(rho-mixed))),"rdms"=>rdms)
    end
    rr=[measure("pure_up",pu*pu'),measure("pure_rotated",pr*pr'),measure("equal_ground_mixture",mixed),measure("rotated_basis_ground_mixture",mixrot)]
    for beta in s["beta"]
        rhoc,_=A.exact_thermal_density(hc,beta);rho=zeros(n^2,n^2);rho[ids,ids]=rhoc
        epsilon=real(1-tr(rho*(vg*vg')))
        push!(rr,measure("canonical_gibbs",rho;beta,epsilon))
    end
    X12=Dict("canonical_dimension"=>length(ids),"canonical_eigenvalues"=>d.values,"ground_dimension"=>length(keep),
        "gap"=>d.values[last(keep)+1]-d.values[1],"eigen_residual"=>norm(hc*d.vectors-d.vectors*Diagonal(d.values)),
        "mixture_rotation_error"=>norm(mixrot-mixed),"spin_lowering_norm"=>lowering_norm,
        "pure_lanczos_residual"=>res,"pure_lanczos_energy"=>pure_E,"rows"=>rr)
    Dict("X09"=>X09,"X12"=>X12)
end
function main()
    began=time();result=Dict{String,Any}("schema"=>CFG["schema"],"status"=>"COMPLETED","manifest_sha256"=>sha(CONFIG),
        "adapter_sha256"=>sha(@__FILE__),"julia_version"=>string(VERSION),"itensors_version"=>string(pkgversion(ITensors)),
        "itensormps_version"=>string(pkgversion(ITensorMPS)),"attribution"=>"Existing user Hubbard MPO, ED Hamiltonian/Lanczos, Gibbs density and RDM routines; ITensorMPS supplies DMRG. New certification and ensemble adapters are explicit.")
    for k in ("ed_source_sha256","thermal_source_sha256","dmrg_source_sha256","project_sha256","environment_manifest_sha256");result[k]=CFG[k];end
    merge!(result,compute());result["elapsed_seconds"]=time()-began
    mkpath(dirname(OUTPUT));open(OUTPUT,"w") do io;TOML.print(io,result);end
end
main()

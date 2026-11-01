#!/usr/bin/env julia
# Additive dynamics/fermion/spin-response adapters; existing user ED/thermal APIs.
using LinearAlgebra, TOML, SHA
BLAS.set_num_threads(1)
length(ARGS)%2==0 || error("Expected --manifest PATH --output PATH")
const OPT=Dict(ARGS[i]=>ARGS[i+1] for i in 1:2:length(ARGS))
const CONFIG=abspath(get(OPT,"--manifest",joinpath(@__DIR__,"manifest.toml")))
const CFG=TOML.parsefile(CONFIG)
const OUTPUT=abspath(OPT["--output"])
ispath(OUTPUT) && error("Output exists")
string(VERSION)==CFG["julia_version"] || error("Unpinned Julia")
const THERMAL=haskey(ENV,"RENYI_ROOT") ? joinpath(ENV["RENYI_ROOT"],"analysis/audit_u4_l4_thermal_ed.jl") : normpath(joinpath(dirname(CONFIG),CFG["thermal_source_relative"]))
const ED=joinpath(dirname(THERMAL),"compute_s2_ed.jl")
bytes2hex(sha256(read(THERMAL)))==CFG["thermal_source_sha256"] || error("Thermal source mismatch")
bytes2hex(sha256(read(ED)))==CFG["ed_source_sha256"] || error("ED source mismatch")
include(THERMAL)
const A=U4L4ThermalEDAudit
const E=A.GroundStateED

function materialize(space)
    n=E.dimension(space);H=zeros(n,n);unit=zeros(n);column=zeros(n)
    for j in 1:n
        fill!(unit,0.);unit[j]=1.;E.apply_hamiltonian!(column,unit,space);H[:,j]=column
    end
    H
end
function blocks(L,t,U,up,down)
    space=E.build_space(L,t,U,0.,up,down;boundary="pbc");H=materialize(space);d=eigen(Symmetric(H))
    space,H,d,norm(H*d.vectors-d.vectors*Diagonal(d.values))
end
function full_indices(space,L)
    [Int(u)*(1<<L)+Int(d)+1 for u in space.up_masks for d in space.dn_masks]
end
function annihilation(L,spin,k)
    n=1<<L;c=zeros(ComplexF64,n^2,n^2)
    for u in 0:n-1,d in 0:n-1,j in 0:L-1
        bit=1<<j;occupied=spin=="up" ? (u&bit)!=0 : (d&bit)!=0
        occupied || continue
        lower=bit-1
        parity=spin=="up" ? count_ones(u&lower) : count_ones(u)+count_ones(d&lower)
        unew=spin=="up" ? xor(u,bit) : u;dnew=spin=="down" ? xor(d,bit) : d
        c[unew*n+dnew+1,u*n+d+1]+=(-1.)^parity*cis(-k*j)/sqrt(L)
    end
    c
end
function clustered(freq,weight,tolerance)
    order=sortperm(freq);groups=Dict{String,Any}[];i=1
    while i<=length(order)
        j=i+1;anchor=freq[order[i]]
        while j<=length(order) && freq[order[j]]-anchor<=tolerance;j+=1;end
        ids=order[i:j-1]
        push!(groups,Dict("omega"=>sum(freq[ids])/length(ids),"weight"=>sum(weight[ids]),
            "multiplicity"=>length(ids),"span"=>maximum(freq[ids])-minimum(freq[ids])))
        i=j
    end
    groups
end

function compute()
    L=CFG["L"];t=CFG["t"];U=CFG["U"];L==4 && CFG["boundary"]=="pbc" || error("Wrong frozen fixture")
    K,V,masks=A.full_fock_matrices(L,t,U,0.;boundary="pbc");n=1<<L
    Nu=[Matrix(Diagonal([Float64((u>>j)&1) for u in masks for d in masks])) for j in 0:L-1]
    Nd=[Matrix(Diagonal([Float64((d>>j)&1) for u in masks for d in masks])) for j in 0:L-1]
    D=[Nu[j]*Nd[j] for j in 1:L];NuT=sum(Nu);NdT=sum(Nd);N=NuT+NdT;DT=sum(D)
    H=Matrix(K+Diagonal(V)+U/2*N-U*L/4*I)
    diagnostics=Dict("H_hermiticity"=>norm(H-H'),"N_commutator"=>norm(H*N-N*H))
    expected(psi,op)=real(dot(psi,op*psi))
    # P08 calls the existing matrix-vector Hamiltonian API in the exact fixed-N sector.
    s=CFG["P08"];space,HH,d,_=blocks(L,t,U,s["Nup"],s["Ndown"])
    occupied=UInt(0)
    for j in s["initial_doublon_sites_one_based"];occupied|=UInt(1)<<(j-1);end
    iu=findfirst(==(occupied),space.up_masks);idn=findfirst(==(occupied),space.dn_masks)
    initial=zeros(E.dimension(space));initial[(iu-1)*length(space.dn_masks)+idn]=1.
    ids=full_indices(space,L);proj=d.vectors'*initial;rows=Dict{String,Any}[]
    for tm in range(s["time_start"],s["time_stop"];length=s["time_count"])
        psi=d.vectors*(exp.(-im*tm.*d.values).*proj)
        density=[expected(psi,(Nu[j]+Nd[j])[ids,ids]) for j in 1:L]
        push!(rows,Dict("time"=>tm,"norm"=>real(dot(psi,psi)),"energy"=>expected(psi,HH),
            "Nup"=>expected(psi,NuT[ids,ids]),"Ndown"=>expected(psi,NdT[ids,ids]),"density"=>density,
            "imbalance"=>sum((iseven(j) ? 1. : -1.)*density[j] for j in 1:L)/(s["Nup"]+s["Ndown"]),
            "double_occupancy_total"=>expected(psi,DT[ids,ids])))
    end
    P08=Dict("rows"=>rows,"initial_energy"=>expected(initial,HH),
        "Hamiltonian_eigen_residual"=>norm(HH*d.vectors-d.vectors*Diagonal(d.values)))
    # P09 additionally exercises the actual Lanczos ground-state API.
    s=CFG["P09"];space,HH,d,er=blocks(L,t,U,s["Nup"],s["Ndown"])
    E0,p0,iterations,reported_residual=E.ground_state_lanczos(space;maxiter=300,tol=1e-12)
    psi=zeros(n^2);psi[full_indices(space,L)]=p0
    ground=Dict("energy"=>E0,"norm"=>real(dot(psi,psi)),"gap_same_sector"=>d.values[2]-d.values[1],
        "residual"=>max(norm(HH*p0-E0*p0),reported_residual,er),
        "Nup"=>expected(psi,NuT),"Ndown"=>expected(psi,NdT),"Lanczos_iterations"=>iterations,
        "dense_ground_overlap_error"=>abs(1-abs2(dot(d.vectors[:,1],p0))),"dense_ground_energy_error"=>abs(E0-d.values[1]))
    rows=Dict{String,Any}[];mu=s["physical_mu_frequency"];G=H-mu*N;g0=E0-mu*(s["Nup"]+s["Ndown"])
    for spin in s["spin"]
        plus=(s["Nup"]+(spin=="up"),s["Ndown"]+(spin=="down"))
        minus=(s["Nup"]-(spin=="up"),s["Ndown"]-(spin=="down"))
        sp,_,dp,rp=blocks(L,t,U,plus...);sm,_,dm,rm=blocks(L,t,U,minus...)
        ip=full_indices(sp,L);imask=full_indices(sm,L)
        for m in s["momentum_indices"]
            k=2pi*m/L;c=annihilation(L,spin,k);add=c'*psi;remove=c*psi
            wa=abs2.(dp.vectors'*add[ip]);wr=abs2.(dm.vectors'*remove[imask])
            fa=dp.values.-E0.-mu;fr=E0.-dm.values.-mu
            Xa=(G-g0*I)*add;Xr=(G-g0*I)*remove
            push!(rows,Dict("spin"=>spin,"momentum_index"=>m,
                "addition"=>clustered(fa,wa,s["degeneracy_cluster_abs"]),"removal"=>clustered(fr,wr,s["degeneracy_cluster_abs"]),
                "addition_dimension"=>length(ip),"removal_dimension"=>length(imask),"addition_eigen_residual"=>rp,"removal_eigen_residual"=>rm,
                "occupation"=>real(dot(remove,remove)),"addition_norm"=>real(dot(add,add)),"sum_weight"=>sum(wa)+sum(wr),
                "first_moment"=>sum(wa.*fa)+sum(wr.*fr),"first_moment_direct"=>real(dot(add,Xa)-dot(remove,Xr)),
                "first_moment_expected"=>-2t*cos(k)+U*expected(psi,spin=="up" ? NdT : NuT)/L-mu,
                "second_moment"=>sum(wa.*fa.^2)+sum(wr.*fr.^2),"second_moment_direct"=>real(dot(Xa,Xa)+dot(Xr,Xr)),
                "anticommutator_matrix_error"=>norm(c*c'+c'*c-I)))
        end
    end
    cross=s["cross_spin_momentum_pair"];cu=annihilation(L,"up",2pi*cross[1]/L);cd=annihilation(L,"down",2pi*cross[2]/L)
    P09=Dict("ground"=>ground,"rows"=>rows,"operator_checks"=>Dict(
        "cross_spin_annihilation_anticommutator"=>norm(cu*cd+cd*cu),
        "cross_spin_mixed_anticommutator"=>norm(cu*cd'+cd'*cu)))
    # P10 calls the existing full-Fock Gibbs-density API at physical-to-centered mu.
    s=CFG["P10"];mu=s["physical_mu"];kin,pot,_=A.full_fock_matrices(L,t,U,mu-U/2;boundary="pbc")
    generator=kin+Diagonal(pot);G=H-mu*N
    rows=Dict{String,Any}[];thermal=Dict{String,Any}[]
    for beta in s["beta"]
        rho,d=A.exact_thermal_density(generator,beta);shift=d.values.-d.values[1];Z=sum(exp.(-beta.*shift))
        direct=exp(-beta.*(G-(minimum(d.values)-U*L/4)*I));direct./=tr(direct)
        means=[real(tr(rho*(Nu[j]-Nd[j])/2)) for j in 1:L]
        moment=real(tr(rho*(N/4-DT/2)))-sum(abs2,means)
        push!(thermal,Dict("beta"=>beta,"mu_physical"=>mu,"mu_centered"=>mu-U/2,"trace"=>real(tr(rho)),
            "rho_expm_error"=>norm(rho-direct),"Nmean"=>real(tr(rho*N)),
            "double_occupancy_total"=>real(tr(rho*DT)),"local_moment_sum"=>moment,"spin_mean_abs_max"=>maximum(abs,means)))
        for m in s["momentum_indices"]
            q=2pi*m/L;S=sum(cis(-q*(j-1))*(Nu[j]-Nd[j])/2 for j in 1:L)/sqrt(L)
            B=d.vectors'*S*d.vectors;weights=abs2.(B);mean=tr(rho*S);offset=abs2(mean)
            endpoint=real(tr(rho*S*S'))-offset
            for fraction in s["tau_fractions"]
                tau=beta*fraction;boltz=exp.(-(beta-tau).*shift.-tau.*shift')./Z
                value=sum(boltz.*weights)-offset
                push!(rows,Dict("beta"=>beta,"momentum_index"=>m,"tau_fraction"=>fraction,"tau"=>tau,"value"=>value,
                    "equal_time_direct"=>endpoint,"spin_mean_abs"=>abs(mean),"spectral_weight_min"=>minimum(weights),"unconnected_value"=>value+offset))
            end
        end
    end
    Dict("P08"=>P08,"P09"=>P09,"P10"=>Dict("rows"=>rows,"thermal"=>thermal),"diagnostics"=>diagnostics)
end
function main()
    began=time();out=Dict{String,Any}("schema"=>CFG["schema"],"status"=>"COMPLETED",
        "manifest_sha256"=>bytes2hex(sha256(read(CONFIG))),"adapter_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "thermal_source_sha256"=>bytes2hex(sha256(read(THERMAL))),"ed_source_sha256"=>bytes2hex(sha256(read(ED))),
        "julia_version"=>string(VERSION),"attribution"=>CFG["new_scope"])
    merge!(out,compute());out["elapsed_seconds"]=time()-began
    mkpath(dirname(OUTPUT));open(OUTPUT,"w") do io;TOML.print(io,out);end
end
main()

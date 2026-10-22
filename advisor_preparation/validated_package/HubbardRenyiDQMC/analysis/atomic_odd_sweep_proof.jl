#!/usr/bin/env julia
# Stdlib proof of the odd-Ltau atomic closed orbit. --kernel additionally diagnoses actual JDQMC Greens.
using LinearAlgebra, Random, Statistics, TOML, Dates, SHA
const SEEDS=[1009,2011,3011,4001,5003,6007,7001,8009]
const L=4;const M=5;const BETA=.5;const DTAU=.1;const U=4.;const MU=.3
const ALPHA=acosh(exp(DTAU*U/2))
if "--kernel" in ARGS
    include(joinpath(ENV["RENYI_ROOT"],"src/HubbardRenyiDQMC.jl"))
end
fieldsum(mask)=2count_ones(mask)-M
weight(mask)=2exp(BETA*MU)*(cosh(BETA*MU)+cosh(ALPHA*fieldsum(mask)))
function scalar_observables(mask)
    q=fieldsum(mask);nu=1/(1+exp(-BETA*MU-ALPHA*q));nd=1/(1+exp(-BETA*MU+ALPHA*q));D=nu*nd;N=nu+nd
    [U*(D-N/2+.25),N,D]
end
alternating(mask)=all(((mask>>(j-1))&1)!=((mask>>j)&1) for j in 1:M-1)
function transition(slice;lazy=.0)
    n=1<<M;P=zeros(n,n)
    for x in 0:n-1
        y=x ⊻ (1<<(slice-1));a=(1-lazy)*min(1.,weight(y)/weight(x))
        P[x+1,y+1]=a;P[x+1,x+1]=1-a
    end
    P
end
function fullcycle(;lazy=.0)
    P=Matrix{Float64}(I,1<<M,1<<M)
    for j in vcat(collect(1:M),collect(M:-1:1));P=P*transition(j;lazy);end
    P
end
function communicating_classes(P)
    n=size(P,1);reach=BitMatrix(P.>0);for i in 1:n;reach[i,i]=true;end
    for k in 1:n,i in 1:n,j in 1:n;reach[i,j]|=reach[i,k]&&reach[k,j];end
    left=Set(1:n);classes=Vector{Int}[]
    while !isempty(left)
        i=first(left);members=[j for j in left if reach[i,j]&&reach[j,i]];push!(classes,sort!(members));setdiff!(left,members)
    end
    [Dict("masks"=>c.-1,"size"=>length(c),"closed"=>all(!any(P[i,j]>0 for j in setdiff(1:n,c)) for i in c)) for c in classes]
end
function exact_proof()
    W=weight.(0:(1<<M)-1);pi=W./sum(W);O=hcat(scalar_observables.(0:(1<<M)-1)...)
    isclosed=alternating.(0:(1<<M)-1);ordinary=.!isclosed
    closedmean=O[:,isclosed]*pi[isclosed]/sum(pi[isclosed]);openmean=O[:,ordinary]*pi[ordinary]/sum(pi[ordinary]);target=O*pi
    Ps=fullcycle();Pl=fullcycle(lazy=.5)
    records=Any[]
    for seed in SEEDS
        s=rand(Xoshiro(seed),-1:2:1,L,M);masks=[sum(Int(s[i,j]==1)<<(j-1) for j in 1:M) for i in 1:L];frozen=[i for i in 1:L if alternating(masks[i])];k=length(frozen);prediction=(L-k).*openmean.+k.*closedmean
        push!(records,Dict("seed"=>seed,"initial_masks"=>masks,"frozen_sites_1based"=>frozen,"predicted_energy_centered"=>prediction[1],"predicted_density_total"=>prediction[2],"predicted_double_occupancy_total"=>prediction[3]))
    end
    Dict("created_utc"=>string(now(UTC)),"parameters"=>Dict("L"=>L,"Ltau"=>M,"U"=>U,"beta"=>BETA,"dtau"=>DTAU,"mu_centered"=>MU,"t"=>0.),"derivation"=>"At t=0 the per-site HS target is W(q)=2exp(beta*mu)[cosh(beta*mu)+cosh(alpha*q)]. Alternating odd-length fields remain |q|=1 under every systematic local flip, so each flip accepts with probability1 and cannot leave the orbit.","ordinary_fullcycle_classes"=>communicating_classes(Ps),"lazy_local_fullcycle_classes"=>communicating_classes(Pl),"ordinary_stationarity_max_error"=>maximum(abs,transpose(pi)*Ps-transpose(pi)),"lazy_stationarity_max_error"=>maximum(abs,transpose(pi)*Pl-transpose(pi)),"lazy_minimum_diagonal_probability"=>minimum(diag(Pl)),"correct_closed_class_mass"=>sum(pi[isclosed]),"uniform_initialization_closed_mass"=>sum(isclosed)/length(isclosed),"exact_target_totals"=>Dict("energy_centered"=>L*target[1],"density_total"=>L*target[2],"double_occupancy_total"=>L*target[3]),"seed_predictions"=>records)
end
function kernel_check()
    H=getfield(@__MODULE__,:HubbardRenyiDQMC);J=H.jdqmcf
    rows=Any[];trajectories=Dict{Int,Any}();warm=512;meas=1024
    for ns in [1,2],seed in SEEDS
        r=H.initializer(0.,U,MU,BETA,L,DTAU,seed;n_stab=ns);initial=copy(r.s);initial_alt=[i for i in 1:L if all(initial[i,j]==-initial[i,j+1] for j in 1:M-1)]
        maxdirect=0.;maxfresh=0.;maxprop=0.;maxdelta=0.;maxcurrent=0.;vals=zeros(3);alt_visits=zeros(Int,L);trace=Vector{UInt8}();bounds=[Inf,-Inf];checks=0
        for sw in 1:warm+meas
            r,delta,_=H.update_replica!(r,0.);maxdelta=max(maxdelta,delta)
            q=vec(sum(r.s,dims=2));Guref=Diagonal(1 ./(1 .+exp.(BETA*MU.+ALPHA.*q)));Gdref=Diagonal(1 ./(1 .+exp.(BETA*MU.-ALPHA.*q)))
            maxdirect=max(maxdirect,maximum(abs,r.Gup_00-Guref),maximum(abs,r.Gdn_00-Gdref));maxcurrent=max(maxcurrent,maximum(abs,r.Gup-Guref),maximum(abs,r.Gdn-Gdref))
            up=ones(L);dn=ones(L)
            for B in r.Bup;up.*=B.expmΔτV;end
            for B in r.Bdn;dn.*=B.expmΔτV;end
            maxprop=max(maxprop,maximum(abs,log.(up).-BETA*MU.-ALPHA.*q),maximum(abs,log.(dn).-BETA*MU.+ALPHA.*q))
            if sw%64==0
                fc=J.FermionGreensCalculator(r.Bup,BETA,DTAU,ns);fresh=zeros(L,L);J.calculate_equaltime_greens!(fresh,fc);maxfresh=max(maxfresh,maximum(abs,r.Gup_00-fresh));checks+=1
            end
            append!(trace,UInt8.(vec(r.s).+1))
            if sw>warm
                Cu=Matrix{Float64}(I,L,L)-transpose(r.Gup_00);Cd=Matrix{Float64}(I,L,L)-transpose(r.Gdn_00);D=sum(diag(Cu).*diag(Cd));N=tr(Cu+Cd);E=U*(D-N/2+L/4);vals.+=[E,N,D];bounds[1]=min(bounds[1],minimum([diag(Cu);diag(Cd)]));bounds[2]=max(bounds[2],maximum([diag(Cu);diag(Cd)]))
                for i in 1:L;all(r.s[i,j]==-r.s[i,j+1] for j in 1:M-1)&&(alt_visits[i]+=1);end
            end
        end
        digest=bytes2hex(sha256(trace));ns==1&&(trajectories[seed]=digest)
        push!(rows,Dict("n_stab"=>ns,"seed"=>seed,"initial_alternating_sites"=>initial_alt,"mean_energy_centered"=>vals[1]/meas,"mean_density_total"=>vals[2]/meas,"mean_double_occupancy_total"=>vals[3]/meas,"maximum_G00_vs_exact_field_error"=>maxdirect,"maximum_Gcurrent_vs_exact_field_error"=>maxcurrent,"maximum_G00_vs_fresh_stabilized_error"=>maxfresh,"fresh_checks"=>checks,"maximum_propagator_vs_field_log_error"=>maxprop,"maximum_reported_deltaG"=>maxdelta,"density_eigenvalue_min"=>bounds[1],"density_eigenvalue_max"=>bounds[2],"measurement_alternating_visits_per_site"=>alt_visits,"trajectory_sha256"=>digest,"matches_nstab1_trajectory"=>get(trajectories,seed,"")==digest))
    end
    rows
end
function main()
    output=get(ENV,"ATOMIC_DIAGNOSTIC_OUTPUT",joinpath(@__DIR__,"atomic_odd_sweep_diagnostic.toml"));isfile(output)&&error("Choose a new ATOMIC_DIAGNOSTIC_OUTPUT")
    r=exact_proof();"--kernel" in ARGS&&(r["kernel_checks"]=kernel_check());r["runner_sha256"]=bytes2hex(sha256(read(@__FILE__)));r["julia_version"]=string(VERSION);r["slurm_job_id"]=get(ENV,"SLURM_JOB_ID","LOCAL")
    open(output,"w") do io;TOML.print(io,r);end
    println("Original full-cycle class sizes: ",[c["size"] for c in r["ordinary_fullcycle_classes"]]);println("Lazy local full-cycle class sizes: ",[c["size"] for c in r["lazy_local_fullcycle_classes"]]);println("Stationarity maxerror original=",r["ordinary_stationarity_max_error"]," lazy=",r["lazy_stationarity_max_error"]);println("Wrote ",output)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end

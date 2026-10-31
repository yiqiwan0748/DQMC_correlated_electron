#!/usr/bin/env julia
# External scientific comparisons only. Internal regressions are a separate suite.
using LinearAlgebra, Random, Statistics, TOML, Dates, SHA, Printf
const ROOT=get(ENV,"RENYI_ROOT",normpath(joinpath(@__DIR__,"../source_snapshot/HubbardRenyiDQMC")))
using HubbardRenyiDQMC
const H=HubbardRenyiDQMC
module ED
include(joinpath(Main.ROOT,"analysis/compute_s2_ed.jl"))
end
module DM
include(joinpath(Main.ROOT,"analysis/compute_s2_itensor_dmrg.jl"))
end
module THERM
include(joinpath(Main.ROOT,"analysis/audit_u4_l4_thermal_ed.jl"))
end
module NODE
include(joinpath(Main.ROOT,"run/run_lambda_node.jl"))
end
using ITensors, ITensorMPS
BLAS.set_num_threads(1)
const CACHE=Dict{Any,Any}()
cached(f,k)=get!(f,CACHE,k)
function sector_expectation(space,psi,c)
    ob=c["observable"];L=space.L;nd=length(space.dn_masks);ans=0.
    if ob=="S2"
        return ED.renyi2_entropy(psi,space,c["nA"])
    end
    for (iu,u) in enumerate(space.up_masks),(id,d) in enumerate(space.dn_masks)
        k=(iu-1)*nd+id;p=abs2(psi[k]);nu=count_ones(u);nn=count_ones(d);D=count_ones(u&d)
        if ob=="double_occupancy_total"
            ans+=p*D
        elseif ob=="density_total"
            ans+=p*(nu+nn)
        elseif ob=="spin_z_correlation"
            i,j=c["sites"];si=Int((u>>i)&1)-Int((d>>i)&1);sj=Int((u>>j)&1)-Int((d>>j)&1);ans+=p*si*sj/4
        elseif ob=="onebody_up"
            i,j=c["sites"]
            if i==j
                ans+=p*Int((u>>i)&1)
            elseif ((u>>j)&1)==1 && ((u>>i)&1)==0
                mask=u ⊻ (UInt(1)<<i) ⊻ (UInt(1)<<j)
                iu2=findfirst(==(mask),space.up_masks)
                sign=ED.fermion_hop_sign(u,j,i)
                ans+=psi[(iu2-1)*nd+id]*psi[k]*sign
            end
        else
            error("Unsupported sector observable $ob")
        end
    end
    ans
end
function user_ed(c)
    f=cached((:ed,c["L"],c["t"],c["U"],c["boundary"],c["Nup"],c["Ndown"])) do
        # The target solver is uncentered and fixed N. Its mu is set to zero.
        space=ED.build_space(c["L"],c["t"],c["U"],0.,c["Nup"],c["Ndown"];boundary=c["boundary"])
        E,psi,it,reported=ED.ground_state_lanczos(space;tol=1e-11,maxiter=1000)
        hv=similar(psi);ED.apply_hamiltonian!(hv,psi,space);res=norm(hv-E*psi)
        res<=1e-9 || error("Explicit ED eigenpair residual $res exceeds 1e-9")
        (;space,E,psi,it,res)
    end
    v=c["observable"]=="energy_uncentered" ? f.E : sector_expectation(f.space,f.psi,c)
    (;value=v,residual=f.res,iterations=f.it)
end
function user_dmrg(c)
    f=cached((:dmrg,c["L"],c["t"],c["U"],c["boundary"])) do
        L=c["L"];sites=siteinds("Electron",L;conserve_qns=true)
        mpo=DM.hubbard_mpo(sites;L,t=c["t"],U=c["U"],mu=0.,boundary=c["boundary"])
        Random.seed!(817291);psi0=random_mps(sites,DM.half_filled_sz0_state(L);linkdims=12)
        E,psi=dmrg(mpo,psi0;nsweeps=16,maxdim=[16,32,64,128],cutoff=1e-12,eigsolve_tol=1e-12,outputlevel=0)
        normcheck=real(inner(psi,psi));E2=real(inner(mpo,psi,mpo,psi));Echeck=real(inner(psi',mpo,psi));variance=E2-Echeck^2
        abs(normcheck-1)<=1e-9 || error("MPS norm $normcheck")
        isfinite(variance) && -1e-9<=variance<=1e-7 || error("MPS energy variance $variance outside [-1e-9,1e-7]")
        abs(E-Echeck)<=1e-8 || error("Reported/contracted MPS energies differ")
        D=sum(expect(psi,"Nupdn"));N=sum(expect(psi,"Ntot"));z=correlation_matrix(psi,"Sz","Sz");cup=correlation_matrix(psi,"Cdagup","Cup")
        abs(N-L)<=1e-8 || error("MPS wrong particle sector")
        (;psi,E=Echeck,D,N,z,cup,variance,normcheck)
    end
    ob=c["observable"]
    v=if ob=="energy_uncentered"
        f.E
    elseif ob=="double_occupancy_total"
        f.D
    elseif ob=="spin_z_correlation"
        i,j=c["sites"].+1;real(f.z[i,j])
    elseif ob=="onebody_up"
        i,j=c["sites"].+1;real(f.cup[i,j])
    elseif ob=="S2"
        first(DM.renyi2_bond_entropy(f.psi,c["nA"]))
    else
        error("Unsupported DMRG observable $ob")
    end
    (;value=v,variance=f.variance,norm=f.normcheck)
end
function fock_observable(rho,masks,L,c)
    ob=c["observable"];d=length(masks);ans=0.
    if ob=="S2"
        return first(THERM.U4L4ThermalEDAudit.renyi2(rho,masks,L,c["nA"]))
    end
    for (iu,u) in enumerate(masks),(id,v) in enumerate(masks)
        k=(iu-1)*d+id;p=real(rho[k,k])
        if ob=="double_occupancy_total"
            ans+=p*count_ones(u&v)
        elseif ob=="density_total"
            ans+=p*(count_ones(u)+count_ones(v))
        elseif ob=="spin_z_correlation"
            i,j=c["sites"];a=Int((u>>i)&1)-Int((v>>i)&1);b=Int((u>>j)&1)-Int((v>>j)&1);ans+=p*a*b/4
        elseif ob=="onebody_up"
            i,j=c["sites"]
            if i==j
                ans+=p*Int((u>>i)&1)
            elseif ((u>>j)&1)==1 && ((u>>i)&1)==0
                u2=u ⊻ (UInt(1)<<i) ⊻ (UInt(1)<<j);iu2=Int(u2)+1
                ans+=real(rho[k,(iu2-1)*d+id])*ED.fermion_hop_sign(u,j,i)
            end
        else
            error("Unsupported thermal observable $ob")
        end
    end
    ans
end
function user_thermal(c)
    f=cached((:thermal,c["L"],c["t"],c["U"],c["mu_centered"],c["beta"],get(c,"dtau",.1),c["boundary"])) do
        a=THERM.U4L4ThermalEDAudit
        K,V,masks=a.full_fock_matrices(c["L"],c["t"],c["U"],c["mu_centered"];boundary=c["boundary"])
        if c["ensemble"]=="grand_canonical_trotter"
            rho,_=a.asymmetric_trotter_density(K,V,c["beta"],c["dtau"])
        else
            rho,_=a.exact_thermal_density(K+Diagonal(V),c["beta"])
        end
        (;K,V,masks,rho)
    end
    ob=c["observable"]
    if ob in ["energy_centered","energy_uncentered"]
        grand=real(tr((f.K+Diagonal(f.V))*f.rho))
        n=fock_observable(f.rho,f.masks,c["L"],merge(c,Dict("observable"=>"density_total")))
        v=grand+c["mu_centered"]*n
        ob=="energy_uncentered" && (v+=c["U"]*n/2-c["U"]*c["L"]/4)
    else
        v=fock_observable(f.rho,f.masks,c["L"],c)
    end
    (;value=v,trace_error=abs(real(tr(f.rho))-1))
end
function analytic_free(c)
    L=c["L"];shift=c["boundary"]=="apbc" ? .5 : 0.
    ks=[2pi*(n+shift)/L for n in 0:L-1];es=-2c["t"].*cos.(ks)
    fs=1 ./(1 .+exp.(c["beta"].*(es.-c["mu_centered"])))
    corr(i,j)=real(sum(exp.(im.*ks.*(i-j)).*fs)/L)
    ob=c["observable"]
    if ob=="energy_uncentered"
        return 2sum(es.*fs)
    elseif ob=="onebody_up"
        return corr(c["sites"]...)
    elseif ob=="S2"
        C=[corr(i,j) for i in 0:c["nA"]-1,j in 0:c["nA"]-1]
        nu=eigvals(Hermitian(C));return -2sum(log.(nu.^2 .+(1 .-nu).^2))
    end
    error("Unknown free observable")
end
function analytic_atom(c)
    U=c["U"];mu=c["mu_centered"];beta=c["beta"]
    e=[U/4,-U/4-mu,-U/4-mu,U/4-2mu];w=exp.(-beta.*(e.-minimum(e)));p=w./sum(w)
    N=c["L"]*(p[2]+p[3]+2p[4]);D=c["L"]*p[4]
    ob=c["observable"]
    if ob=="energy_centered"
        return U*(D-N/2+c["L"]/4) # internal centered energy; no -mu*N
    elseif ob=="density_total"
        return N
    elseif ob=="double_occupancy_total"
        return D
    elseif ob=="S2"
        return -c["nA"]*log(sum(p.^2))
    end
    error("Unknown atom observable")
end
function user_free(c)
    r=cached((:free,c["L"],c["t"],c["U"],c["mu_centered"],c["beta"],c["dtau"],c["boundary"])) do
        H.initializer(c["t"],0.,c["mu_centered"],c["beta"],c["L"],c["dtau"],19271;n_stab=2,boundary_condition=Symbol(c["boundary"]))
    end
    Cu=Matrix{Float64}(I,c["L"],c["L"])-transpose(r.Gup_00);Cd=Matrix{Float64}(I,c["L"],c["L"])-transpose(r.Gdn_00)
    ob=c["observable"]
    if ob=="energy_uncentered"
        K=H.hubbard_chain_kinetic_matrix(c["t"],c["L"],Symbol(c["boundary"]));v=real(tr(K*(Cu+Cd)))
    elseif ob=="onebody_up"
        i,j=c["sites"].+1;v=Cu[i,j]
    elseif ob=="S2"
        A=collect(1:c["nA"]);gu=H.grover_matrix(r.Gup_00,r.Gup_00,A);gd=H.grover_matrix(r.Gdn_00,r.Gdn_00,A)
        su=logabsdet(gu);sd=logabsdet(gd);su[2]>0 && sd[2]>0 || error("Nonpositive free Grover determinant");v=-su[1]-sd[1]
    end
    (;value=v)
end
# Monte Carlo arms use eight independent chains; the protocol is fixed before execution.
const MC_SEEDS=[1009,2011,3011,4001,5003,6007,7001,8009]
const NMEAS=parse(Int,get(ENV,"VALIDATION_MC_MEASUREMENTS","4096"))
const NBIN=64
const BURN=512
function measurements(r,c)
    L=c["L"];Cu=Matrix{Float64}(I,L,L)-transpose(r.Gup_00);Cd=Matrix{Float64}(I,L,L)-transpose(r.Gdn_00)
    n=diag(Cu)+diag(Cd);D=sum(diag(Cu).*diag(Cd));N=sum(n)
    K=H.hubbard_chain_kinetic_matrix(c["t"],L,Symbol(c["boundary"]));kin=real(tr(K*(Cu+Cd)))
    z(i,j)=real(((Cu[i,i]-Cd[i,i])*(Cu[j,j]-Cd[j,j])-Cu[i,j]*Cu[j,i]-Cd[i,j]*Cd[j,i]+(i==j ? Cu[i,i]+Cd[i,i] : 0))/4)
    [kin+c["U"]*D,kin+c["U"]*(D-N/2+L/4),N,D,z(1,2),z(1,3),Cu[1,2],Cu[1,3]]
end
function ratio_jackknife(nums,dens)
    B=length(nums);sn=sum(nums);sd=sum(dens);abs(sd)>1e-12 || error("Unresolved sign denominator")
    center=sn/sd;loo=[(sn-nums[i])/(sd-dens[i]) for i in 1:B]
    corrected=B*center-(B-1)*mean(loo)
    (;value=corrected,se=sqrt((B-1)/B*sum((loo.-mean(loo)).^2)))
end
function effective_bins(x)
    n=length(x);y=x.-mean(x);g0=sum(abs2,y)/n
    g0<=eps(Float64)*max(1.,maximum(abs2,x)) && return Float64(n)
    tau=1.
    for lag in 1:min(n÷4,16)
        r=dot(view(y,1:n-lag),view(y,1+lag:n))/((n-lag)*g0)
        r<=0 && break
        tau+=2r
    end
    n/tau
end
# Split-chain variance diagnostic on paired jackknife pseudovalues. This is an
# additional rejection gate, not a proof of equilibration or a confidence bound.
function split_rhat(x::AbstractMatrix)
    n,m=size(x);iseven(n) || error("Split Rhat requires an even number of blocks")
    h=n÷2;y=hcat(x[1:h,:],x[h+1:end,:])
    W=mean(var(y;dims=1,corrected=true));B=h*var(vec(mean(y;dims=1));corrected=true)
    tiny=1e-24*max(1.,maximum(abs2,y))
    W<=tiny && return B<=tiny ? 1. : Inf
    sqrt(max(1.,((h-1)/h*W+B/h)/W))
end
function scalar_mc(c)
    cached((:mc,c["L"],c["t"],c["U"],c["mu_centered"],c["beta"],c["dtau"],c["boundary"],NMEAS)) do
        NMEAS%NBIN==0 || error("VALIDATION_MC_MEASUREMENTS must be divisible by64")
        seedmeans=zeros(length(MC_SEEDS),8);blockses=zeros(length(MC_SEEDS),8);doubleblockses=similar(blockses);effective=similar(blockses)
        bins_all=Any[];pseudovalues=zeros(NBIN,length(MC_SEEDS),8)
        for (si,seed) in enumerate(MC_SEEDS)
            r=H.initializer(c["t"],c["U"],c["mu_centered"],c["beta"],c["L"],c["dtau"],seed;n_stab=2,boundary_condition=Symbol(c["boundary"]))
            for _ in 1:BURN;r,_,_=H.update_replica!(r,0.);end
            nums=zeros(NBIN,8);dens=zeros(NBIN);bs=NMEAS÷NBIN
            for b in 1:NBIN,j in 1:bs
                r,_,_=H.update_replica!(r,0.);sgn=real(r.sgndetGup*r.sgndetGdn);nums[b,:].+=sgn.*measurements(r,c);dens[b]+=sgn
            end
            for k in 1:8
                a=ratio_jackknife(nums[:,k],dens);seedmeans[si,k]=a.value;blockses[si,k]=a.se;effective[si,k]=effective_bins(nums[:,k].-a.value.*dens)
                sn=sum(nums[:,k]);sd=sum(dens);center=sn/sd
                pseudovalues[:,si,k]=[NBIN*center-(NBIN-1)*(sn-nums[b,k])/(sd-dens[b]) for b in 1:NBIN]
                a2=ratio_jackknife(vec(sum(reshape(nums[:,k],2,:),dims=1)),vec(sum(reshape(dens,2,:),dims=1)));doubleblockses[si,k]=a2.se
            end
            push!(bins_all,Dict("seed"=>seed,"numerator_bins"=>[collect(nums[b,:]) for b in 1:NBIN],"denominator_bins"=>dens))
        end
        rhats=[split_rhat(pseudovalues[:,:,k]) for k in 1:8]
        (;seedmeans,blockses,doubleblockses,effective,bins_all,rhats)
    end
end
function entropy_mc(c)
    cached((:mc_entropy,c["L"],c["t"],c["U"],c["mu_centered"],c["beta"],c["dtau"],c["boundary"],c["nA"],NMEAS)) do
        seedmeans=Float64[];raw=Any[];blockses=Float64[];doubleblockses=Float64[];effective=Float64[];Ninc=4
        pseudovalues=zeros(NBIN,length(MC_SEEDS))
        for (si,seed) in enumerate(MC_SEEDS)
            config=(L=c["L"],t=c["t"],U=c["U"],mu=c["mu_centered"],beta=c["beta"],dtau=c["dtau"],n_stab=2,N_sweeps=1,boundary_condition=Symbol(c["boundary"]),nA=c["nA"],Ninc=Ninc,S2_ED=NaN,N_burnin=BURN,N_measurements=NMEAS,N_bins=NBIN)
            results=[NODE.run_lambda_node_simulation(config,k,seed) for k in 0:Ninc-1]
            function entropy_estimate(drop=0;pairblocks=false)
                val=0.
                for r in results
                    nums=pairblocks ? vec(sum(reshape(r.avg_step_plus,2,:),dims=1)) : r.avg_step_plus
                    dens=pairblocks ? vec(sum(reshape(r.avg_phase,2,:),dims=1)) : r.avg_phase
                    numerator=sum(nums)-(drop>0 ? nums[drop] : 0.);denominator=sum(dens)-(drop>0 ? dens[drop] : 0.)
                    ratio=numerator/denominator;isfinite(ratio)&&ratio>0 || error("Nonpositive/unresolved entropy ratio");val-=log(ratio)
                end
                val
            end
            push!(effective,minimum(effective_bins(r.avg_step_plus.-(sum(r.avg_step_plus)/sum(r.avg_phase)).*r.avg_phase) for r in results))
            center=entropy_estimate();loo=[entropy_estimate(b) for b in 1:NBIN];corrected=NBIN*center-(NBIN-1)*mean(loo)
            pseudovalues[:,si]=NBIN*center .-(NBIN-1).*loo
            push!(seedmeans,corrected);push!(blockses,sqrt((NBIN-1)/NBIN*sum((loo.-mean(loo)).^2)))
            loo2=[entropy_estimate(b;pairblocks=true) for b in 1:NBIN÷2];push!(doubleblockses,sqrt((NBIN÷2-1)/(NBIN÷2)*sum((loo2.-mean(loo2)).^2)))
            push!(raw,Dict("seed"=>seed,"nodes"=>[Dict("lambda_index"=>r.lambda_index,"numerator_bins"=>r.avg_step_plus,"denominator_bins"=>r.avg_phase) for r in results]))
        end
        (;seedmeans,blockses,doubleblockses,effective,raw,rhat=split_rhat(pseudovalues))
    end
end
# Exact elementary Student t7 tail: integrate cos(theta)^6 under x=sqrt(7)tan(theta).
function t7critical(alpha)
    tail(q)=begin th=atan(q/sqrt(7));.5-(10th+7.5sin(2th)+1.5sin(4th)+sin(6th)/6)/(10pi) end
    lo=0.;hi=100.
    for _ in 1:100;m=(lo+hi)/2;tail(m)>alpha/2 ? (lo=m) : (hi=m);end
    (lo+hi)/2
end
function user_mc(c)
    ob=c["observable"]
    if ob=="S2"
        f=entropy_mc(c);vals=f.seedmeans;ses=f.blockses;ses2=f.doubleblockses;raw=f.raw;eff=sum(f.effective);rhat=f.rhat
    else
        f=scalar_mc(c);k=Dict("energy_uncentered"=>1,"energy_centered"=>2,"density_total"=>3,"double_occupancy_total"=>4,"spin_z_correlation"=>(c["sites"]==[0,1] ? 5 : 6),"onebody_up"=>(c["sites"]==[0,1] ? 7 : 8))[ob]
        vals=f.seedmeans[:,k];ses=f.blockses[:,k];ses2=f.doubleblockses[:,k];raw=f.bins_all;eff=sum(f.effective[:,k]);rhat=f.rhats[k]
    end
    se=std(vals)/sqrt(length(vals));change=maximum(abs.(ses2.-ses)./max.(ses,1e-12))
    (;value=mean(vals),standard_error=se,seed_values=collect(vals),blocking_max_relative_change=change,raw=raw,critical_value=t7critical(.01/400),independent_seeds=length(vals),measurements_per_seed=NMEAS,effective_bins=eff,split_rhat=rhat)
end
function cleanvalue(x)
    x isa AbstractFloat && !isfinite(x) ? string(x) : x
end
function report_arm(method,c,reference;run_mc=true)
    start=time()
    if method=="DQMC" && c["family"]!="analytic_free" && !run_mc
        return Dict("method"=>method,"status"=>"UNRESOLVED","reason"=>"Monte Carlo arm was explicitly disabled","runtime_seconds"=>0.)
    end
    try
        r=method=="ED" ? user_ed(c) : method=="DMRG" ? user_dmrg(c) : method=="thermal_ED" ? user_thermal(c) : c["family"]=="analytic_free" ? user_free(c) : user_mc(c)
        value=Float64(r.value);err=abs(value-reference);result=Dict{String,Any}("method"=>method,"value"=>value,"reference"=>reference,"absolute_error"=>err,"runtime_seconds"=>time()-start)
        if hasproperty(r,:standard_error)
            half=r.critical_value*r.standard_error;margin=c["equivalence_margin"]
            # Comparing average independent-chain error scales is a diagnostic only; strict all-chain stability avoids hiding one stuck chain.
            good=isfinite(value)&&isfinite(half)&&r.blocking_max_relative_change<=.5 && r.effective_bins>=128 && r.split_rhat<=1.05
            result["status"]=good && err+half<=margin ? "PASS" : "UNRESOLVED"
            result["standard_error"]=r.standard_error;result["simultaneous_half_width"]=half;result["critical_value"]=r.critical_value;result["equivalence_margin"]=margin;result["blocking_max_relative_change"]=r.blocking_max_relative_change;result["independent_seeds"]=r.independent_seeds;result["measurements_per_seed"]=NMEAS;result["seed_values"]=r.seed_values;result["estimated_effective_bins"]=r.effective_bins;result["raw_bins"]=r.raw
            result["split_rhat"]=cleanvalue(r.split_rhat);result["split_rhat_limit"]=1.05
            good && err-half>margin && (result["status"]="FAIL")
        else
            tol=method=="DMRG" ? c["absolute_tolerance_dmrg"] : c["absolute_tolerance_ed"]
            result["tolerance"]=tol;result["status"]=isfinite(value)&&err<=tol ? "PASS" : "FAIL"
            for key in [:residual,:variance,:norm,:trace_error]
                hasproperty(r,key) && (result[String(key)]=getproperty(r,key))
            end
        end
        result
    catch ex
        Dict("method"=>method,"status"=>"UNRESOLVED","reason"=>sprint(showerror,ex,catch_backtrace()),"runtime_seconds"=>time()-start)
    end
end
function main(args=ARGS)
    manifest=get(ENV,"VALIDATION_MANIFEST",joinpath(@__DIR__,"benchmark_manifest.toml"))
    oraclepath=get(ENV,"VALIDATION_ORACLES",joinpath(@__DIR__,"../external_oracles/oracle_values.toml"))
    output=get(ENV,"VALIDATION_OUTPUT",joinpath(@__DIR__,"external_validation_results.toml"))
    cases=TOML.parsefile(manifest)["cases"]
    if "--list" in args
        for c in cases;println(c["id"]," ",c["family"]," ",c["observable"]);end
        return
    end
    isfile(output) && error("Output already exists; choose a new VALIDATION_OUTPUT path")
    values=isfile(oraclepath) ? TOML.parsefile(oraclepath)["values"] : Dict{String,Any}()
    filterprefix=get(ENV,"VALIDATION_ID_PREFIX","");run_mc=get(ENV,"VALIDATION_RUN_DQMC","1")=="1"
    selected=filter(!isempty,strip.(split(get(ENV,"VALIDATION_ID_LIST",""),',')))
    result=Dict{String,Any}("created_utc"=>string(now(UTC)),"julia_version"=>string(VERSION),"manifest_sha256"=>bytes2hex(sha256(read(manifest))),"oracle_file"=>oraclepath,"oracle_sha256"=>isfile(oraclepath) ? bytes2hex(sha256(read(oraclepath))) : "MISSING","source_root"=>ROOT,"source_hashes"=>Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in ["Project.toml","Manifest.toml","src/HubbardRenyiDQMC.jl","src/dqmc_kernel.jl","src/grover_rank1_update.jl","analysis/compute_s2_ed.jl","analysis/compute_s2_itensor_dmrg.jl","analysis/audit_u4_l4_thermal_ed.jl"]),"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"),"host"=>get(ENV,"HOSTNAME","UNKNOWN"),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"dmrg_controls"=>Dict("seed"=>817291,"nsweeps"=>16,"maxdim_schedule"=>[16,32,64,128],"cutoff"=>1e-12,"eigsolve_tolerance"=>1e-12),"run_monte_carlo"=>run_mc,"mc_measurements_per_seed"=>NMEAS,"cases"=>Any[])
    result["runner_sha256"]=bytes2hex(sha256(read(@__FILE__)))
    result["mc_protocol_revision"]="lazy-local-proposals-v1; split-Rhat<=1.05; alpha=0.01/400"
    result["lambda_runner_sha256"]=bytes2hex(sha256(read(joinpath(ROOT,"run/run_lambda_node.jl"))))
    result["all_src_sha256"]=Dict(relpath(joinpath(d,f),ROOT)=>bytes2hex(sha256(read(joinpath(d,f)))) for (d,_,fs) in walkdir(joinpath(ROOT,"src")) for f in fs if endswith(f,".jl"))
    for c in cases
        startswith(c["id"],filterprefix)||continue
        isempty(selected) || c["id"] in selected || continue
        if c["family"]=="analytic_free"
            ref=analytic_free(c)
        elseif c["family"]=="analytic_atom"
            ref=analytic_atom(c)
        elseif haskey(values,c["id"])
            ref=Float64(values[c["id"]])
        else
            push!(result["cases"],Dict("id"=>c["id"],"status"=>"UNRESOLVED","reason"=>"Independent external oracle unavailable"));continue
        end
        arms=[report_arm(m,c,ref;run_mc) for m in c["target_methods"]]
        status=all(a["status"]=="PASS" for a in arms) ? "PASS" : any(a["status"]=="FAIL" for a in arms) ? "FAIL" : "UNRESOLVED"
        push!(result["cases"],Dict("id"=>c["id"],"family"=>c["family"],"observable"=>c["observable"],"status"=>status,"reference"=>ref,"arms"=>arms))
        println(c["id"]," ",status," ",join([a["method"]*":"*a["status"] for a in arms],", "));flush(stdout)
        open(output,"w") do io;TOML.print(io,result);end
    end
    result["finished_utc"]=string(now(UTC));result["pass_count"]=count(c["status"]=="PASS" for c in result["cases"]);result["fail_count"]=count(c["status"]=="FAIL" for c in result["cases"]);result["unresolved_count"]=count(c["status"]=="UNRESOLVED" for c in result["cases"])
    open(output,"w") do io;TOML.print(io,result);end
    println("Cases=",length(result["cases"])," PASS=",result["pass_count"]," FAIL=",result["fail_count"]," UNRESOLVED=",result["unresolved_count"])
    result["fail_count"]>0 && exit(1)
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

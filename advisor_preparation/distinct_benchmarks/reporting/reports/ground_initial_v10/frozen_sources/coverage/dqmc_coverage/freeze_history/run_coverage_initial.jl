#!/usr/bin/env julia
# Frozen estimator calibration; no source patch, fitted correction, or replicate rejection.
using LinearAlgebra, Statistics, Random, TOML, SHA, Dates
using HubbardRenyiDQMC
const H=HubbardRenyiDQMC
BLAS.set_num_threads(1)
const ROOT=realpath(ENV["RENYI_ROOT"])
const CFGFILE=joinpath(@__DIR__,"manifest.toml")
const CFG=TOML.parsefile(CFGFILE)
const SEEDSFILE=joinpath(@__DIR__,"seeds.toml")
const SEEDS=TOML.parsefile(SEEDSFILE)
const MODE=get(ENV,"COVERAGE_MODE","pilot")
MODE in ("pilot","final") || error("COVERAGE_MODE must be pilot/final")
const COUNT=CFG[MODE=="pilot" ? "pilot_replicates" : "final_replicates"]
const FIRST=parse(Int,get(ENV,"COVERAGE_FIRST","1"))
const LAST=parse(Int,get(ENV,"COVERAGE_LAST",string(COUNT)))
1<=FIRST<=LAST<=COUNT || error("Replicate range invalid")
const OUTPUT=ENV["COVERAGE_OUTPUT"]
ispath(OUTPUT)&&error("Choose a new output path")
realpath(pathof(H))==realpath(joinpath(ROOT,"src/HubbardRenyiDQMC.jl")) || error("Active package/source root mismatch")
sha(p)=bytes2hex(sha256(read(p)))
function source_receipt()
    frozen=TOML.parsefile(joinpath(@__DIR__,"source_hashes.toml"))["files"]
    actual=Dict(rel=>sha(joinpath(ROOT,rel)) for rel in keys(frozen))
    actual==frozen || error("Source hash differs from frozen corrected package")
    actual
end
const SOURCE_HASHES=source_receipt()
const ORACLEPATH=ENV["COVERAGE_ORACLE"]
sha(ORACLEPATH)==CFG["oracle_sha256"] || error("External oracle bytes differ")
const REFERENCE=TOML.parsefile(ORACLEPATH)["values"][CFG["oracle_case"]]
REFERENCE==CFG["oracle_reference"] || error("Oracle reference mismatch")
function ratio_jackknife(nums,dens)
    B=length(nums);sn=sum(nums);sd=sum(dens)
    abs(sd)>1e-12 || error("Unresolved sign denominator")
    center=sn/sd;loo=[(sn-nums[b])/(sd-dens[b]) for b in 1:B]
    corrected=B*center-(B-1)*mean(loo)
    (;value=corrected,se=sqrt((B-1)/B*sum((loo.-mean(loo)).^2)),pseudo=B*center .-(B-1).*loo)
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
function split_rhat(x)
    n,m=size(x);h=n÷2;y=hcat(x[1:h,:],x[h+1:end,:])
    W=mean(var(y;dims=1,corrected=true));B=h*var(vec(mean(y;dims=1));corrected=true)
    tiny=1e-24*max(1.,maximum(abs2,y))
    W<=tiny && return B<=tiny ? 1. : Inf
    sqrt(max(1.,((h-1)/h*W+B/h)/W))
end
function energy(r,c,K)
    L=c["L"];Cu=Matrix{Float64}(I,L,L)-transpose(r.Gup_00);Cd=Matrix{Float64}(I,L,L)-transpose(r.Gdn_00)
    real(tr(K*(Cu+Cd)))+c["U"]*sum(diag(Cu).*diag(Cd))
end
safe(x)=isfinite(x) ? x : string(x)
function replicate(index)
    started=time();c=CFG["fixture"];B=CFG["bins"];C=CFG["chains_per_replicate"]
    N=CFG["measured_sweeps"];N%B==0||error("Bin count");ns=N÷B
    means=zeros(C);ses=zeros(C);ses2=zeros(C);eff=zeros(C);pseudo=zeros(B,C);raw=Any[]
    seeds=SEEDS[MODE][string(index)];K=H.hubbard_chain_kinetic_matrix(c["t"],c["L"],Symbol(c["boundary"]))
    for (chain,seed) in enumerate(seeds)
        r=H.initializer(c["t"],c["U"],c["mu_centered"],c["beta"],c["L"],c["dtau"],seed;n_stab=c["n_stab"],boundary_condition=Symbol(c["boundary"]))
        for _ in 1:CFG["burn_sweeps"];r,_,_=H.update_replica!(r,0.);end
        nums=zeros(B);dens=zeros(B)
        for b in 1:B,_ in 1:ns
            r,_,_=H.update_replica!(r,0.);sgn=real(r.sgndetGup*r.sgndetGdn)
            nums[b]+=sgn*energy(r,c,K);dens[b]+=sgn
        end
        a=ratio_jackknife(nums,dens);means[chain]=a.value;ses[chain]=a.se;pseudo[:,chain]=a.pseudo
        a2=ratio_jackknife(vec(sum(reshape(nums,2,:),dims=1)),vec(sum(reshape(dens,2,:),dims=1)))
        ses2[chain]=a2.se;eff[chain]=effective_bins(nums.-a.value.*dens)
        push!(raw,Dict("chain"=>chain,"seed"=>seed,"numerator_bins"=>nums,"denominator_bins"=>dens))
    end
    value=mean(means);se=std(means)/sqrt(C);half=CFG["student_t7_critical"]*se
    change=maximum(abs.(ses2.-ses)./max.(ses,1e-12));rh=split_rhat(pseudo)
    finite=all(isfinite,[value,se,half,change,rh,sum(eff)])
    good=finite && change<=CFG["blocking_relative_change_max"] && rh<=CFG["split_rhat_max"] && sum(eff)>=CFG["effective_bins_min"]
    Dict("id"=>index,"mode"=>MODE,"quality_pass"=>good,"value"=>safe(value),"standard_error"=>safe(se),"half_width_95"=>safe(half),"covered"=>(finite&&abs(value-REFERENCE)<=half),"chain_means"=>means,"chain_block_se"=>ses,"chain_double_block_se"=>ses2,"chain_effective_bins"=>eff,"blocking_max_relative_change"=>safe(change),"split_rhat"=>safe(rh),"raw_bins"=>raw,"elapsed_seconds"=>time()-started)
end
function main()
    mkpath(dirname(OUTPUT));started=time()
    meta=Dict("schema"=>CFG["schema"],"catalogue_id"=>"X06","mode"=>MODE,"first_replicate"=>FIRST,"last_replicate"=>LAST,"manifest_sha256"=>sha(CFGFILE),"seeds_sha256"=>sha(SEEDSFILE),"driver_sha256"=>sha(@__FILE__),"oracle_sha256"=>sha(ORACLEPATH),"reference"=>REFERENCE,"source_hashes"=>SOURCE_HASHES,"julia_version"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),"job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"),"created_utc"=>string(now(UTC)))
    open(OUTPUT,"w") do io;TOML.print(io,Dict("metadata"=>meta));end
    for index in FIRST:LAST
        row=replicate(index)
        open(OUTPUT,"a") do io;write(io,"\n");TOML.print(io,Dict("replicates"=>[row]));end
        println("X06 ",MODE," replicate ",index," quality=",row["quality_pass"]," covered=",row["covered"]," seconds=",row["elapsed_seconds"]);flush(stdout)
    end
    open(OUTPUT,"a") do io
        write(io,"\n");TOML.print(io,Dict("completion"=>Dict("status"=>MODE=="pilot" ? "PILOT_ONLY" : "SHARD_COMPLETE","replicate_count"=>LAST-FIRST+1,"elapsed_seconds"=>time()-started)))
    end
end
main()

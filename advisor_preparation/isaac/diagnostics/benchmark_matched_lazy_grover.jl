#!/usr/bin/env julia
# Dated, matched WITH-versus-WITHOUT optimization benchmark. No archived sampler is edited.
using LinearAlgebra, Random, Statistics, TOML, Dates, SHA, Sockets, Printf
using HubbardRenyiDQMC
const H=HubbardRenyiDQMC
BLAS.set_num_threads(1)
const ROOT=normpath(get(ENV,"RENYI_ROOT",dirname(dirname(pathof(H)))))
const DEFAULT_SIZES=[8,16,32,64,128]
const FIXED_SEED=940031
const BETA=2.;const DTAU=.05;const U=2.;const MU=0.;const T=1.;const LAMBDA=.5;const D_LAMBDA=1/8;const NSTAB=5
const ROUNDS=parse(Int,get(ENV,"TIMING_ROUNDS","12"))
const REPEATS=parse(Int,get(ENV,"TIMING_REPEATS","7"))
const WARMUP_ROUNDS=2
const BUDGET_SECONDS=parse(Float64,get(ENV,"TIMING_BUDGET_SECONDS","600"))
function make_replicas(L)
    a=H.initializer(T,U,MU,BETA,L,DTAU,FIXED_SEED;n_stab=NSTAB,boundary_condition=:pbc)
    b=H.initializer(T,U,MU,BETA,L,DTAU,FIXED_SEED+1;n_stab=NSTAB,boundary_condition=:pbc)
    a,b
end
function measure_node(a,b,C)
    gu=H.grover_matrix(a.Gup_00,b.Gup_00,C);gd=H.grover_matrix(a.Gdn_00,b.Gdn_00,C)
    lu,su=logabsdet(gu);ld,sd=logabsdet(gd);su*sd>0||error("Nonpositive physical Grover estimator")
    logX=lu+ld;sgn=real(a.sgndetGup*a.sgndetGdn*b.sgndetGup*b.sgndetGdn)
    # Same node-level forward/backward and endpoint measurements for both paths.
    [sgn,sgn*exp(D_LAMBDA*logX),sgn*exp(-D_LAMBDA*logX),sgn*exp(-LAMBDA*logX),sgn*exp((1-LAMBDA)*logX)]
end
function node_loop(update,a,b,C,rounds)
    acceptance=zeros(rounds,2);observations=zeros(rounds,5)
    # Cheap sweep-boundary field snapshots allow exact replay checks. Their cost
    # is shared by both variants and is explicitly included in node-loop time.
    snapshots=Vector{Vector{Int8}}();kernel_seconds=0.;measurement_seconds=0.;delta_max=0.
    for k in 1:rounds
        stamp=time_ns();a,delta,aa=update(a,b,0.,C,LAMBDA);b,delta2,ab=update(b,a,0.,C,LAMBDA);kernel_seconds+=(time_ns()-stamp)/1e9
        delta_max=max(delta_max,delta,delta2);acceptance[k,:]=[aa,ab]
        stamp=time_ns();observations[k,:]=measure_node(a,b,C);measurement_seconds+=(time_ns()-stamp)/1e9
        push!(snapshots,vcat(Int8.(vec(a.s)),Int8.(vec(b.s))))
    end
    (;a,b,acceptance,observations,snapshots,kernel_seconds,measurement_seconds,delta_max)
end
function compare_runs(d,o)
    field_equal=d.snapshots==o.snapshots
    acceptance_equal=d.acceptance==o.acceptance
    signs_equal=all(getproperty(dr,key)==getproperty(orr,key) for (dr,orr) in [(d.a,o.a),(d.b,o.b)] for key in [:sgndetGup,:sgndetGdn])
    next_rng_equal=rand(copy(d.a.rng),UInt64)==rand(copy(o.a.rng),UInt64)&&rand(copy(d.b.rng),UInt64)==rand(copy(o.b.rng),UInt64)
    green_error=maximum(maximum(abs,getproperty(dr,key)-getproperty(orr,key)) for (dr,orr) in [(d.a,o.a),(d.b,o.b)] for key in [:Gup,:Gdn,:Gup_00,:Gdn_00,:Gup_τ0,:Gup_0τ,:Gdn_τ0,:Gdn_0τ])
    determinant_log_error=maximum(abs(getproperty(dr,key)-getproperty(orr,key)) for (dr,orr) in [(d.a,o.a),(d.b,o.b)] for key in [:logdetGup,:logdetGdn])
    measurement_error=maximum(abs,d.observations-o.observations)
    relative_measurement_error=maximum(abs.(d.observations-o.observations)./max.(1.,abs.(d.observations)))
    valid=field_equal&&acceptance_equal&&signs_equal&&next_rng_equal&&green_error<=2e-8&&determinant_log_error<=2e-8&&relative_measurement_error<=2e-8&&max(d.delta_max,o.delta_max)<=1e-7
    Dict("status"=>valid ? "PASS" : "FAIL","sweep_boundary_HS_fields_exact"=>field_equal,"sweep_acceptance_rates_exact"=>acceptance_equal,"determinant_signs_exact"=>signs_equal,"next_RNG_draws_exact"=>next_rng_equal,"maximum_Green_absolute_error"=>green_error,"maximum_logdet_absolute_error"=>determinant_log_error,"maximum_measurement_absolute_error"=>measurement_error,"maximum_measurement_relative_error"=>relative_measurement_error,"maximum_reported_deltaG"=>max(d.delta_max,o.delta_max))
end
function row_from_timing(label,r,trial,order)
    Dict("variant"=>label,"repeat"=>r,"execution_order"=>order,"node_loop_seconds"=>trial.time,"update_kernel_seconds"=>trial.value.kernel_seconds,"measurement_seconds"=>trial.value.measurement_seconds,"other_loop_seconds"=>max(0.,trial.time-trial.value.kernel_seconds-trial.value.measurement_seconds),"allocated_bytes"=>trial.bytes,"gc_seconds"=>trial.gctime)
end
function paired_bootstrap(d,o)
    rng=Xoshiro(41911);B=2000;n=length(d);ratios=Float64[]
    for _ in 1:B;idx=rand(rng,1:n,n);push!(ratios,median(d[idx])/median(o[idx]));end
    quantile(ratios,[.025,.975])
end
function benchmark_size(L;deadline=Inf,progress=nothing)
    C=collect(1:L÷2);init_start=time();basea,baseb=make_replicas(L);initialization_seconds=time()-init_start
    # Shared direct warm-up gives both variants exactly the same starting state.
    warm_start=time();w=node_loop(H.update_replica_fused!,basea,baseb,C,WARMUP_ROUNDS);basea,baseb=w.a,w.b;warmup_seconds=time()-warm_start
    raw=Any[];checks=Any[];d_times=Float64[];o_times=Float64[];dk=Float64[];ok=Float64[]
    for r in 1:REPEATS
        time()>deadline && break
        da,db=deepcopy(basea),deepcopy(baseb);oa,ob=deepcopy(basea),deepcopy(baseb)
        # Deep copies and explicit GC are outside the measured node workflow.
        direct=nothing;optimized=nothing
        if isodd(r)
            GC.gc();direct=@timed node_loop(H.update_replica_fused!,da,db,C,ROUNDS)
            GC.gc();optimized=@timed node_loop(H.update_replica_rank1_fused!,oa,ob,C,ROUNDS)
        else
            GC.gc();optimized=@timed node_loop(H.update_replica_rank1_fused!,oa,ob,C,ROUNDS)
            GC.gc();direct=@timed node_loop(H.update_replica_fused!,da,db,C,ROUNDS)
        end
        check=compare_runs(direct.value,optimized.value);check["repeat"]=r;push!(checks,check)
        push!(raw,row_from_timing("without_Grover_cache_direct_fused",r,direct,isodd(r) ? 1 : 2));push!(raw,row_from_timing("with_rank1_Grover_cache",r,optimized,isodd(r) ? 2 : 1))
        push!(d_times,direct.time);push!(o_times,optimized.time);push!(dk,direct.value.kernel_seconds);push!(ok,optimized.value.kernel_seconds)
        if progress!==nothing
            progress(Dict("L"=>L,"nA"=>length(C),"beta"=>BETA,"dtau"=>DTAU,"U"=>U,"lambda"=>LAMBDA,"n_stab"=>NSTAB,"two_replica_update_rounds"=>ROUNDS,"status"=>"RUNNING","completed_repeats"=>r,"requested_repeats"=>REPEATS,"raw_timings"=>deepcopy(raw),"correctness_checks"=>deepcopy(checks)))
        end
        @printf("L=%d nA=%d repeat=%d direct=%.6fs optimized=%.6fs ratio=%.4fx correctness=%s\n",L,length(C),r,direct.time,optimized.time,direct.time/optimized.time,check["status"]);flush(stdout)
    end
    isempty(d_times)&&return Dict("L"=>L,"nA"=>length(C),"status"=>"UNRESOLVED","reason"=>"Soft timing budget exhausted before any repeat","raw_timings"=>raw)
    valid=all(c["status"]=="PASS" for c in checks)
    completed=length(d_times)
    status=!valid ? "FAIL" : completed<REPEATS ? "UNRESOLVED" : "PASS"
    Dict("L"=>L,"nA"=>length(C),"beta"=>BETA,"dtau"=>DTAU,"Ltau"=>round(Int,BETA/DTAU),"U"=>U,"t"=>T,"mu_centered"=>MU,"boundary"=>"pbc","lambda"=>LAMBDA,"delta_lambda"=>D_LAMBDA,"n_stab"=>NSTAB,"two_replica_update_rounds"=>ROUNDS,"warmup_rounds_excluded"=>WARMUP_ROUNDS,"requested_repeats"=>REPEATS,"completed_repeats"=>completed,"initialization_seconds_excluded"=>initialization_seconds,"common_warmup_seconds_excluded"=>warmup_seconds,"status"=>status,"median_without_optimization_seconds"=>median(d_times),"median_with_optimization_seconds"=>median(o_times),"ratio_of_median_node_times"=>median(d_times)/median(o_times),"paired_bootstrap_95pct_ratio_interval"=>paired_bootstrap(d_times,o_times),"ratio_of_median_update_kernel_times"=>median(dk)/median(ok),"raw_timings"=>raw,"correctness_checks"=>checks)
end
function main()
    output=get(ENV,"TIMING_OUTPUT",joinpath(@__DIR__,"matched_lazy_grover_timing.toml"));isfile(output)&&error("Choose a new TIMING_OUTPUT")
    Threads.nthreads()==1||error("Run with JULIA_NUM_THREADS=1 for this registered benchmark")
    H.LOCAL_PROPOSAL_IDLE_PROBABILITY==.5||error("Both variants must use corrected per-local-proposal laziness")
    normpath(dirname(pathof(H)))==normpath(joinpath(ROOT,"src"))||error("Loaded package source differs from RENYI_ROOT")
    REPEATS>=3&&ROUNDS>=2||error("Use at least3 repeats and2 rounds")
    sizes=haskey(ENV,"TIMING_SIZES") ? parse.(Int,split(ENV["TIMING_SIZES"],',')) : DEFAULT_SIZES
    cpu=isfile("/proc/cpuinfo") ? first(filter(l->startswith(l,"model name"),split(read("/proc/cpuinfo",String),'\n'))) : Sys.CPU_NAME
    report=Dict{String,Any}("created_utc"=>string(now(UTC)),"comparison"=>"WITH versus WITHOUT rank-one Grover caching in the same locally lazy two-replica Renyi lambda-node workflow","scope"=>"Node-loop times include two replica updates, the same ratio measurements and sweep-boundary snapshots; initialization/common warmup/deepcopy/explicit GC/compilation are excluded and separately described. This is not a full lambda-grid entropy campaign.","primary_statistic"=>"ratio of median node-loop elapsed times at each preregistered size; all sizes retained","timing_uncertainty"=>"Paired bootstrap interval describes repeat-to-repeat timing jitter on this host, not hardware-independent performance uncertainty","source_root"=>ROOT,"julia_version"=>string(VERSION),"host"=>gethostname(),"cpu"=>cpu,"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"blas_config"=>string(BLAS.get_config()),"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"),"seeds"=>[FIXED_SEED,FIXED_SEED+1],"registered_sizes"=>sizes,"soft_budget_seconds"=>BUDGET_SECONDS,"runner_sha256"=>bytes2hex(sha256(read(@__FILE__))),"source_hashes"=>Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in ["Project.toml","Manifest.toml","src/HubbardRenyiDQMC.jl","src/dqmc_kernel.jl","src/grover_rank1_update.jl"]),"results"=>Any[])
    # Compile both execution paths before any reported timing, using the same types.
    compilation_start=time();a,b=make_replicas(4);node_loop(H.update_replica_fused!,deepcopy(a),deepcopy(b),[1,2],2);node_loop(H.update_replica_rank1_fused!,deepcopy(a),deepcopy(b),[1,2],2);report["compilation_and_initial_probe_seconds_excluded"]=time()-compilation_start
    deadline=time()+BUDGET_SECONDS
    for L in sizes
        try
            if time()>deadline
                push!(report["results"],Dict("L"=>L,"status"=>"UNRESOLVED","reason"=>"Registered soft budget exhausted before this size"))
                continue
            end
            function progress(partial)
                if !isempty(report["results"])&&get(last(report["results"]),"L",0)==L
                    report["results"][end]=partial
                else
                    push!(report["results"],partial)
                end
                open(output,"w") do io;TOML.print(io,report);end
            end
            result=benchmark_size(L;deadline,progress)
            if !isempty(report["results"])&&get(last(report["results"]),"L",0)==L
                report["results"][end]=result
            else
                push!(report["results"],result)
            end
        catch ex
            push!(report["results"],Dict("L"=>L,"status"=>"UNRESOLVED","reason"=>sprint(showerror,ex,catch_backtrace())))
        end
        open(output,"w") do io;TOML.print(io,report);end
    end
    report["finished_utc"]=string(now(UTC));open(output,"w") do io;TOML.print(io,report);end
    println("Wrote ",output)
    for row in report["results"];haskey(row,"ratio_of_median_node_times")&&@printf("L=%d nA=%d median node speedup %.4fx status=%s\n",row["L"],row["nA"],row["ratio_of_median_node_times"],row["status"]);end
    any(r["status"]!="PASS" for r in report["results"])&&exit(1)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end

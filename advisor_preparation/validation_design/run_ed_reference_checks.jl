#!/usr/bin/env julia
# Preliminary ED arms of the canonical external100 manifest. No JDQMC/ITensor/Pkg import.
using LinearAlgebra, Random, Statistics, TOML, Dates, SHA, Printf
const ROOT=get(ENV,"RENYI_ROOT",normpath(joinpath(@__DIR__,"../source_snapshot/HubbardRenyiDQMC")))
module ED
include(joinpath(Main.ROOT,"analysis/compute_s2_ed.jl"))
end
module THERM
include(joinpath(Main.ROOT,"analysis/audit_u4_l4_thermal_ed.jl"))
end
BLAS.set_num_threads(1)
const CACHE=Dict{Any,Any}()
cached(f,k)=get!(f,CACHE,k)
include(joinpath(@__DIR__,"ed_reference_helpers.jl"))
function main()
    manifest=get(ENV,"VALIDATION_MANIFEST",joinpath(@__DIR__,"benchmark_manifest.toml"))
    oraclepath=get(ENV,"VALIDATION_ORACLES",joinpath(@__DIR__,"../external_oracles/oracle_values.toml"))
    output=get(ENV,"VALIDATION_OUTPUT",joinpath(@__DIR__,"ed_arm_results.toml"))
    isfile(output)&&error("Output exists; choose a new VALIDATION_OUTPUT")
    cases=TOML.parsefile(manifest)["cases"];oracles=isfile(oraclepath) ? TOML.parsefile(oraclepath)["values"] : Dict{String,Any}()
    filehash(path)=bytes2hex(sha256(read(path)))
    result=Dict{String,Any}("scope"=>"Preliminary ED arms only; no complete multi-method case is certified", "created_utc"=>string(now(UTC)),"julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"),"host"=>get(ENV,"HOSTNAME","UNKNOWN"),"manifest_sha256"=>filehash(manifest),"oracle_file"=>oraclepath,"oracle_sha256"=>isfile(oraclepath) ? filehash(oraclepath) : "MISSING","runner_sha256"=>filehash(@__FILE__),"helper_sha256"=>filehash(joinpath(@__DIR__,"ed_reference_helpers.jl")),"canonical_consumer_sha256"=>filehash(joinpath(@__DIR__,"run_external_validation.jl")),"source_root"=>ROOT,"source_hashes"=>Dict(p=>filehash(joinpath(ROOT,p)) for p in ["analysis/compute_s2_ed.jl","analysis/audit_u4_l4_thermal_ed.jl"]),"cases"=>Any[])
    for c in cases
        id=c["id"];family=c["family"];available=family=="QuSpin_ground" ? "ED" : family in ["QuSpin_thermal_trotter","analytic_atom"] ? "thermal_ED" : ""
        ref=family=="analytic_atom" ? analytic_atom(c) : get(oracles,id,nothing)
        arms=Any[]
        for method in c["target_methods"]
            if method!=available
                push!(arms,Dict("method"=>method,"status"=>"UNRESOLVED","reason"=>"Requires full dependency-backed runner"));continue
            elseif ref===nothing
                push!(arms,Dict("method"=>method,"status"=>"UNRESOLVED","reason"=>"Independent oracle unavailable"));continue
            end
            start=time()
            try
                r=method=="ED" ? user_ed(c) : user_thermal(c);value=Float64(r.value);err=abs(value-ref);tol=c["absolute_tolerance_ed"]
                arm=Dict{String,Any}("method"=>method,"status"=>isfinite(value)&&err<=tol ? "PASS" : "FAIL","value"=>value,"reference"=>Float64(ref),"absolute_error"=>err,"tolerance"=>tol,"runtime_seconds"=>time()-start)
                for key in [:residual,:iterations,:trace_error];hasproperty(r,key)&&(arm[String(key)]=getproperty(r,key));end
                push!(arms,arm)
            catch ex
                push!(arms,Dict("method"=>method,"status"=>"UNRESOLVED","reason"=>sprint(showerror,ex,catch_backtrace()),"runtime_seconds"=>time()-start))
            end
        end
        status=any(a["status"]=="FAIL" for a in arms) ? "FAIL" : any(a["status"]=="PASS" for a in arms) ? "PARTIAL" : "UNRESOLVED"
        push!(result["cases"],Dict("id"=>id,"family"=>family,"observable"=>c["observable"],"status"=>status,"arms"=>arms))
        println(id," ",status," ",join([a["method"]*":"*a["status"] for a in arms],", "));flush(stdout)
        open(output,"w") do io;TOML.print(io,result);end
    end
    arms=[a for c in result["cases"] for a in c["arms"]]
    result["finished_utc"]=string(now(UTC));result["ed_arm_pass_count"]=count(a["status"]=="PASS" for a in arms);result["ed_arm_fail_count"]=count(a["status"]=="FAIL" for a in arms);result["unresolved_arm_count"]=count(a["status"]=="UNRESOLVED" for a in arms);result["complete_case_pass_count"]=0
    open(output,"w") do io;TOML.print(io,result);end
    println("ED arms PASS=",result["ed_arm_pass_count"]," FAIL=",result["ed_arm_fail_count"]," unresolved/pending arms=",result["unresolved_arm_count"]," complete cases certified=0")
    result["ed_arm_fail_count"]>0 && exit(1)
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

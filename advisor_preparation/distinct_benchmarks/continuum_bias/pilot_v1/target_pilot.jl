#!/usr/bin/env julia
# Only the excluded timing pilot. Actual source kernels/native node runner unchanged.
using LinearAlgebra, TOML, SHA, Dates
BLAS.set_num_threads(1)
const HERE=@__DIR__
const CFGFILE=joinpath(HERE,"manifest.toml")
const CFG=TOML.parsefile(CFGFILE)
const SEEDFILE=joinpath(HERE,CFG["seed_file"])
const SEEDS=TOML.parsefile(SEEDFILE)
const ROOT=realpath(ENV["RENYI_ROOT"])
const OUTPUT=ENV["X01_PILOT_OUTPUT"]
sha(p)=bytes2hex(sha256(read(p)))
CFG["phase"]=="PILOT_ONLY"&&!CFG["scientific_completion_allowed"]&&!CFG["production_campaign_implemented"]||error("Timing pilot only")
CFG["pilot_dtau"]==[.25,.03125]&&CFG["bundles_per_dtau"]==2&&CFG["measured_sweeps"]==1024&&CFG["bins"]==64||error("Unregistered pilot grid")
ispath(OUTPUT)&&error("Use a new output")
string(VERSION)==CFG["julia_version"]||error("Julia runtime mismatch")
const SOURCE_HASHES=Dict(rel=>sha(joinpath(ROOT,rel)) for rel in keys(CFG["source_hashes"]))
SOURCE_HASHES==CFG["source_hashes"]||error("Source/Project/Manifest pin differs before include")
module EnergySource
include(joinpath(Main.ROOT,"src/HubbardRenyiDQMC.jl"))
end
const H=EnergySource.HubbardRenyiDQMC
module NodeSource
include(joinpath(Main.ROOT,"run/run_lambda_node.jl"))
end
string(pkgversion(H.jdqmcf))==CFG["jdqmcframework_version"]||error("JDQMCFramework runtime version")

function save_new(path,data)
    ispath(path)&&error("Receipt exists")
    mkpath(dirname(path));open(path,"w") do io;TOML.print(io,data);end
end
function append_row(key,row,filename)
    save_new(joinpath(dirname(OUTPUT),filename),row)
    open(OUTPUT,"a") do io;write(io,"\n");TOML.print(io,Dict(key=>[row]));end
end
function validate_seeds()
    SEEDS["namespace"]==CFG["seed_namespace"]||error("Seed namespace")
    expected=Set((i,b) for i in eachindex(CFG["pilot_dtau"]) for b in 1:2)
    rows=SEEDS["bundles"];Set((x["dtau_index"],x["bundle"]) for x in rows)==expected&&length(rows)==4||error("Incomplete/duplicate seeds")
    actual=Int[]
    for r in rows
        r["dtau"]==CFG["pilot_dtau"][r["dtau_index"]]||error("Seed/time-step mismatch")
        push!(actual,r["energy_seed"])
        for k in CFG["lambda_indices"],replica in 1:2
            derived=NodeSource.derived_seed(r["entropy_root_seed"],k,replica)
            derived==r["node_replica_seeds"][k+1][replica]||error("Native derived seed mismatch")
            push!(actual,derived)
        end
    end
    length(actual)==36&&length(unique(actual))==36&&all(>(0),actual)||error("Actual initialization seeds collide")
    nothing
end
function energy_chain(s)
    start=time();dt=s["dtau"];L=CFG["L"]
    r=H.initializer(CFG["t"],CFG["U"],CFG["mu_centered"],CFG["beta"],L,dt,s["energy_seed"];n_stab=CFG["n_stab"],boundary_condition=Symbol(CFG["boundary"]))
    deltaG=0.;acceptance=0.
    for _ in 1:CFG["burn_sweeps"]
        r,deltaG,a=H.update_replica!(r,deltaG);acceptance+=a
    end
    after_burn=deltaG;B=CFG["bins"];bs=CFG["measured_sweeps"]÷B
    nums=zeros(B);dens=zeros(B);nnums=zeros(B);deltas=zeros(B)
    K=H.hubbard_chain_kinetic_matrix(CFG["t"],L,Symbol(CFG["boundary"]))
    for b in 1:B
        for _ in 1:bs
            # Pass the returned cumulative maximum into every following update.
            r,deltaG,a=H.update_replica!(r,deltaG);acceptance+=a
            cu=Matrix{Float64}(I,L,L)-transpose(r.Gup_00)
            cd=Matrix{Float64}(I,L,L)-transpose(r.Gdn_00)
            physical=real(tr(K*(cu+cd)))+CFG["U"]*sum(diag(cu).*diag(cd))
            phase=real(r.sgndetGup*r.sgndetGdn)
            nums[b]+=phase*physical;dens[b]+=phase;nnums[b]+=phase*sum(diag(cu)+diag(cd))
        end
        nums[b]/=bs;dens[b]/=bs;nnums[b]/=bs;deltas[b]=deltaG
    end
    Dict("status"=>"COMPLETED","dtau_index"=>s["dtau_index"],"dtau"=>dt,"bundle"=>s["bundle"],
        "seed"=>s["energy_seed"],"numerator_bins"=>nums,"denominator_bins"=>dens,"N_numerator_bins"=>nnums,
        "deltaG_after_burnin"=>after_burn,"cumulative_deltaG_bins"=>deltas,"cumulative_deltaG"=>deltaG,
        "acceptance"=>acceptance/(CFG["burn_sweeps"]+CFG["measured_sweeps"]),"elapsed_seconds"=>time()-start)
end
function entropy_node(s,k)
    config=(L=CFG["L"],t=CFG["t"],U=CFG["U"],mu=CFG["mu_centered"],beta=CFG["beta"],dtau=s["dtau"],
        n_stab=CFG["n_stab"],N_sweeps=CFG["sweeps_per_measurement"],boundary_condition=Symbol(CFG["boundary"]),
        nA=CFG["nA"],Ninc=CFG["lambda_increments"],S2_ED=NaN,N_burnin=CFG["burn_sweeps"],
        N_measurements=CFG["measured_sweeps"],N_bins=CFG["bins"])
    r=NodeSource.run_lambda_node_simulation(config,k,s["entropy_root_seed"])
    # r.deltaG is the native runner's cumulative maximum over both replicas/burn/measurements.
    Dict("status"=>"COMPLETED","dtau_index"=>s["dtau_index"],"dtau"=>s["dtau"],"bundle"=>s["bundle"],
        "lambda_index"=>r.lambda_index,"lambda_value"=>r.lambda_value,"delta_lambda"=>r.delta_lambda,
        "root_seed"=>r.root_seed,"seed1"=>r.seed1,"seed2"=>r.seed2,"cut_sites_one_based"=>r.C,
        "numerator_bins"=>r.avg_step_plus,"denominator_bins"=>r.avg_phase,
        "native_avg_den_bins"=>r.avg_den,"native_avg_num_bins"=>r.avg_num,
        "native_cumulative_deltaG"=>r.deltaG,"acceptance"=>r.acceptance,"elapsed_seconds"=>r.runtime,
        "native_plus_stats"=>Dict(string(k)=>v for (k,v) in pairs(r.plus_stats)))
end
function main()
    validate_seeds();start=time();errors=String[]
    meta=Dict("schema"=>CFG["schema"],"phase"=>"PILOT_ONLY","catalogue_id"=>"X01","scientific_completions"=>0,
        "manifest_sha256"=>sha(CFGFILE),"seeds_sha256"=>sha(SEEDFILE),"adapter_sha256"=>sha(@__FILE__),
        "source_hashes"=>SOURCE_HASHES,"julia_version"=>string(VERSION),"jdqmcframework_version"=>string(pkgversion(H.jdqmcf)),
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"job_id"=>get(ENV,"SLURM_JOB_ID","NONE"),
        "created_utc"=>string(now(UTC)),"source_root"=>ROOT,
        "sampling"=>Dict(k=>CFG[k] for k in ("pilot_dtau","bundles_per_dtau","burn_sweeps","measured_sweeps","bins","n_stab","sweeps_per_measurement","lambda_indices","lambda_values","delta_lambda","cut_sites_one_based")),
        "attribution"=>"Actual corrected update_replica! and unchanged native run_lambda_node_simulation; new scalar measurement, bin-retention and timing adapters only. No production continuum fit.")
    save_new(OUTPUT,Dict("metadata"=>meta))
    for s in SEEDS["bundles"]
        prefix="dt$(s["dtau_index"])_bundle$(s["bundle"])"
        row=try energy_chain(s) catch ex
            push!(errors,sprint(showerror,ex));Dict("status"=>"EXECUTION_ERROR","dtau_index"=>s["dtau_index"],"bundle"=>s["bundle"],"error"=>last(errors))
        end
        append_row("energy_rows",row,prefix*"_energy.toml")
        for k in CFG["lambda_indices"]
            row=try entropy_node(s,k) catch ex
                push!(errors,sprint(showerror,ex));Dict("status"=>"EXECUTION_ERROR","dtau_index"=>s["dtau_index"],"bundle"=>s["bundle"],"lambda_index"=>k,"error"=>last(errors))
            end
            append_row("entropy_node_rows",row,prefix*"_lambda$(k).toml")
        end
        println("X01 timing pilot finished bundle ",prefix);flush(stdout)
    end
    open(OUTPUT,"a") do io
        write(io,"\n");TOML.print(io,Dict("completion"=>Dict("status"=>isempty(errors) ? "PILOT_ONLY_COMPLETE" : "EXECUTION_ERROR",
            "scientific_completions"=>0,"elapsed_seconds"=>time()-start,"errors"=>errors)))
    end
    isempty(errors)||exit(1)
end
main()

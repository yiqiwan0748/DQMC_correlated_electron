using Profile
using Printf
using Dates

repo = pwd()

param = ".\\parameters\\L8_U2_beta32_nA4_Ninc4_quick.toml"
lambda_index = 2

warm_seed = 91000
profile_seed = 91001

warm_out = ".\\results\\raw\\profile\\current_lambda_node_L8_beta32_k2_warmup"
profile_out = ".\\results\\raw\\profile\\current_lambda_node_L8_beta32_k2_profiled"

profile_tree = ".\\results\\processed\\profile_current_lambda_node_L8_beta32_k2_tree.txt"
profile_flat = ".\\results\\processed\\profile_current_lambda_node_L8_beta32_k2_flat.txt"

function set_run_args(param, lambda_index, seed, outdir)
    empty!(ARGS)
    append!(
        ARGS,
        [
            param,
            "--lambda-index", string(lambda_index),
            "--root-seed", string(seed),
            "--output-dir", outdir,
        ],
    )
    return nothing
end

function has_summary_csv(outdir)
    return isdir(outdir) && any(endswith(name, "_summary.csv") for name in readdir(outdir))
end

function clean_outdir(outdir)
    rm(outdir; force = true, recursive = true)
    mkpath(dirname(outdir))
    return nothing
end

mkpath(".\\results\\raw\\profile")
mkpath(".\\results\\processed")

println("========================================")
println("Current Renyi lambda-node CPU profile")
println("========================================")
println("parameter file = ", param)
println("lambda index   = ", lambda_index)
println("warm seed      = ", warm_seed)
println("profile seed   = ", profile_seed)
println("started        = ", Dates.now())
println("========================================")

# Load run_lambda_node.jl and warm up. If the included file auto-runs main(),
# the warm-up output will appear. If it does not auto-run, call main() manually.
clean_outdir(warm_out)
set_run_args(param, lambda_index, warm_seed, warm_out)

println()
println("Including run/run_lambda_node.jl for warm-up...")
include(joinpath(repo, "run", "run_lambda_node.jl"))

if !has_summary_csv(warm_out)
    @assert isdefined(Main, :main)
    println()
    println("Calling main() for warm-up...")
    main()
end

println()
println("Warm-up complete.")

# Profile a second run in the same Julia process.
clean_outdir(profile_out)
set_run_args(param, lambda_index, profile_seed, profile_out)

println()
println("Starting profiled run...")
Profile.init(; delay = 0.001, n = 10_000_000)
Profile.clear()

t0 = time()
@profile main()
wall_seconds = time() - t0

println()
@printf("Profiled run wall_seconds = %.6f\n", wall_seconds)

open(profile_tree, "w") do io
    println(io, "Current Renyi lambda-node CPU profile")
    println(io, "format = tree")
    println(io, "parameter_file = ", param)
    println(io, "lambda_index = ", lambda_index)
    println(io, "seed = ", profile_seed)
    @printf(io, "wall_seconds = %.6f\n", wall_seconds)
    println(io)

    try
        Profile.print(io; format = :tree, C = false, combine = true)
    catch err
        println(io, "Profile.print tree mode failed; falling back to default print.")
        println(io, err)
        println(io)
        Profile.print(io)
    end
end

open(profile_flat, "w") do io
    println(io, "Current Renyi lambda-node CPU profile")
    println(io, "format = flat")
    println(io, "parameter_file = ", param)
    println(io, "lambda_index = ", lambda_index)
    println(io, "seed = ", profile_seed)
    @printf(io, "wall_seconds = %.6f\n", wall_seconds)
    println(io)

    try
        Profile.print(io; format = :flat, C = false, combine = true, sortedby = :count)
    catch err
        println(io, "Profile.print flat sorted mode failed; falling back to flat/default print.")
        println(io, err)
        println(io)

        try
            Profile.print(io; format = :flat, C = false, combine = true)
        catch err2
            println(io, "Profile.print flat mode failed; falling back to default print.")
            println(io, err2)
            println(io)
            Profile.print(io)
        end
    end
end

println()
println("Saved profile files:")
println(profile_tree)
println(profile_flat)
println()
println("finished = ", Dates.now())
#!/usr/bin/env julia

using LinearAlgebra
using Printf

const PROJECT_ROOT = dirname(dirname(@__FILE__))
const BASELINE_COMMIT = "7c7cca6"
const BASELINE_SOURCE =
    "archive/validated_old_pipeline/run_lambda_grid_step_validated.jl"
const BASELINE_BLOB = "6b45732cd6a1339cca113bfea1c9b1b165ae4af8"
const RUN_NAME = "L4_U4_beta16_nA2_Ninc16_untouched_7c7cca6"
const ROOT_SEEDS = [1234, 11234, 21234, 31234]
const NINC = 16
const N_BINS = 40
const BINS_HEADER =
    "run_name,lambda_index,lambda,root_seed,bin,phase,den,num,step_plus,step_minus"

module UntouchedBaseline end

function load_baseline!()
    spec = "$BASELINE_COMMIT:$BASELINE_SOURCE"
    blob = strip(read(`git -C $PROJECT_ROOT rev-parse $spec`, String))
    blob == BASELINE_BLOB || error("Unexpected baseline blob: $blob")
    source = read(`git -C $PROJECT_ROOT show $spec`, String)
    Base.include_string(UntouchedBaseline, source, "git:$spec")
    return nothing
end

function splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = xor(z, z >> 30) * 0xbf58476d1ce4e5b9
    z = xor(z, z >> 27) * 0x94d049bb133111eb
    return xor(z, z >> 31)
end

function derived_seed(root_seed::Integer, lambda_index::Integer, replica_id::Integer)
    x = UInt64(root_seed)
    x = xor(x, UInt64(lambda_index + 1) * 0x9e3779b97f4a7c15)
    x = xor(x, UInt64(replica_id) * 0xbf58476d1ce4e5b9)
    return Int(splitmix64(x) % UInt64(typemax(Int)))
end

function parse_cli(args)
    root_seed = nothing
    output_dir = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--root-seed"
            i += 1
            i <= length(args) || error("--root-seed requires a value")
            root_seed = parse(Int, args[i])
        elseif arg == "--output-dir"
            i += 1
            i <= length(args) || error("--output-dir requires a value")
            output_dir = abspath(args[i])
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end
    root_seed in ROOT_SEEDS || error("root seed must be one of $ROOT_SEEDS")
    output_dir === nothing && error("Missing --output-dir")
    return (; root_seed, output_dir)
end

function node_path(output_dir, lambda_index, root_seed)
    return joinpath(
        output_dir,
        @sprintf("%s_k%03d_seed%d_bins.csv", RUN_NAME, lambda_index, root_seed),
    )
end

function valid_node(path, lambda_index, root_seed)
    isfile(path) || return false
    lines = readlines(path)
    length(lines) == N_BINS + 1 || return false
    lines[1] == BINS_HEADER || return false
    for (bin, line) in enumerate(lines[2:end])
        fields = split(line, ','; keepempty = true)
        length(fields) == 10 || return false
        fields[1] == RUN_NAME || return false
        parse(Int, fields[2]) == lambda_index || return false
        parse(Int, fields[4]) == root_seed || return false
        parse(Int, fields[5]) == bin || return false
        all(isfinite, parse.(Float64, fields[6:9])) || return false
        isempty(fields[10]) || return false
    end
    return true
end

function write_node_atomic(path, result, lambda_index, root_seed)
    ispath(path) && error("Refusing to overwrite existing output: $path")
    temporary, io = mktemp(dirname(path))
    try
        println(io, BINS_HEADER)
        lambda_value = lambda_index / NINC
        for bin in 1:N_BINS
            println(io, join((
                RUN_NAME,
                lambda_index,
                @sprintf("%.12f", lambda_value),
                root_seed,
                bin,
                @sprintf("%.17g", result.avg_phase[bin]),
                @sprintf("%.17g", result.avg_den[bin]),
                @sprintf("%.17g", result.avg_num[bin]),
                @sprintf("%.17g", result.avg_step[bin]),
                "",
            ), ','))
        end
        close(io)
        mv(temporary, path)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary; force = true)
    end
    return nothing
end

function run_root(root_seed, output_dir)
    mkpath(output_dir)
    for lambda_index in 0:(NINC - 1)
        output = node_path(output_dir, lambda_index, root_seed)
        if ispath(output)
            valid_node(output, lambda_index, root_seed) ||
                error("Existing node is invalid: $output")
            println("Skipping valid node k=$lambda_index root=$root_seed")
            continue
        end

        lambda_value = lambda_index / NINC
        seed1 = derived_seed(root_seed, lambda_index, 1)
        seed2 = derived_seed(root_seed, lambda_index, 2)
        println("Starting k=$lambda_index root=$root_seed")
        flush(stdout)
        result = Base.invokelatest(
            UntouchedBaseline.run_power_increment_seed,
            lambda_value,
            seed1,
            seed2;
            Δλ = 1 / NINC,
            L = 4,
            t = 1.0,
            U = 4.0,
            μ = 0.0,
            β = 16.0,
            Δτ = 0.05,
            n_stab = 10,
            N_sweeps = 1,
            nA = 2,
            N_burnin = 8000,
            N_measurements = 4000,
            N_bins = N_BINS,
        )
        write_node_atomic(output, result, lambda_index, root_seed)
        println(@sprintf(
            "Finished k=%d root=%d runtime=%.3f seconds",
            lambda_index,
            root_seed,
            result.runtime,
        ))
        flush(stdout)
    end
end

function main(args = ARGS)
    cli = parse_cli(args)
    BLAS.set_num_threads(1)
    load_baseline!()
    run_root(cli.root_seed, cli.output_dir)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

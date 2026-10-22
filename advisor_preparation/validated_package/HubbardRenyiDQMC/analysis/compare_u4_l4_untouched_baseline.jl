#!/usr/bin/env julia

using Printf

module GridSummary
include(joinpath(@__DIR__, "summarize_grid.jl"))
end

const NINC = 16
const N_BINS = 40
const ROOT_SEEDS = [1234, 11234, 21234, 31234]
const BASELINE_SOURCE =
    "git:7c7cca6:6b45732cd6a1339cca113bfea1c9b1b165ae4af8"

function read_forward_bins(directory)
    rows = Dict{Tuple{Int, Int, Int}, NamedTuple}()
    files = sort(filter(
        name -> endswith(name, "_bins.csv"),
        readdir(directory; join = true),
    ))

    for file in files
        lines = readlines(file)
        isempty(lines) && error("Empty bins file: $file")
        lines[1] == GridSummary.BINS_HEADER || error("Wrong header: $file")
        for line in lines[2:end]
            fields = split(line, ','; keepempty = true)
            length(fields) == 10 || error("Malformed row: $file")
            k = parse(Int, fields[2])
            k < NINC || continue
            root_seed = parse(Int, fields[4])
            bin = parse(Int, fields[5])
            key = (k, root_seed, bin)
            haskey(rows, key) && error("Duplicate row: $key")
            rows[key] = (
                phase = parse(Float64, fields[6]),
                den = parse(Float64, fields[7]),
                num = parse(Float64, fields[8]),
                plus = parse(Float64, fields[9]),
            )
        end
    end

    expected = Set(
        (k, root_seed, bin)
        for k in 0:(NINC - 1)
        for root_seed in ROOT_SEEDS
        for bin in 1:N_BINS
    )
    actual = Set(keys(rows))
    actual == expected || error(
        "Forward-grid coverage mismatch: missing=$(length(setdiff(expected, actual))) " *
        "unexpected=$(length(setdiff(actual, expected)))",
    )
    return rows
end

function summarize_forward(rows)
    step_stats = Dict{Int, NamedTuple}()
    for k in 0:(NINC - 1)
        keys_k = [
            (k, root_seed, bin)
            for root_seed in ROOT_SEEDS
            for bin in 1:N_BINS
        ]
        phases = [rows[key].phase for key in keys_k]
        weights = [rows[key].plus for key in keys_k]
        step_stats[k] = GridSummary.ratio_statistics(weights, phases)
    end
    S2 = sum(step_stats[k].entropy for k in 0:(NINC - 1))
    error = sqrt(sum(step_stats[k].entropy_err^2 for k in 0:(NINC - 1)))
    return (; S2, error, step_stats)
end

function maximum_differences(baseline, optimized)
    maxima = Dict(field => 0.0 for field in (:phase, :den, :num, :plus))
    for key in keys(baseline), field in keys(maxima)
        difference = abs(getproperty(baseline[key], field) - getproperty(optimized[key], field))
        maxima[field] = max(maxima[field], difference)
    end
    return maxima
end

function read_audit(path)
    lines = readlines(path)
    length(lines) == 2 || error("Expected one audit row: $path")
    header = split(lines[1], ',')
    values = split(lines[2], ',')
    length(header) == length(values) || error("Malformed audit CSV: $path")
    row = Dict(header[i] => values[i] for i in eachindex(header))
    return (
        ground = parse(Float64, row["fixed_ground_S2"]),
        thermal = parse(Float64, row["thermal_exact_S2"]),
        trotter = parse(Float64, row["trotter_asymmetric_S2"]),
    )
end

function write_comparison(path, baseline, optimized, audit, maxima)
    ispath(path) && error("Refusing to overwrite existing output: $path")
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        println(
            io,
            "method,S2,error,delta_vs_trotter,z_vs_trotter,source," *
            "max_abs_phase_vs_optimized,max_abs_den_vs_optimized," *
            "max_abs_num_vs_optimized,max_abs_step_plus_vs_optimized",
        )

        function write_dqmc(method, result, source, differences)
            delta = result.S2 - audit.trotter
            println(io, join((
                method,
                @sprintf("%.15e", result.S2),
                @sprintf("%.15e", result.error),
                @sprintf("%.15e", delta),
                @sprintf("%.8f", delta / result.error),
                source,
                differences...,
            ), ','))
        end

        formatted_maxima = (
            @sprintf("%.15e", maxima[:phase]),
            @sprintf("%.15e", maxima[:den]),
            @sprintf("%.15e", maxima[:num]),
            @sprintf("%.15e", maxima[:plus]),
        )
        write_dqmc("DQMC_untouched", baseline, BASELINE_SOURCE, formatted_maxima)
        write_dqmc(
            "DQMC_optimized",
            optimized,
            "completed_U4_production_bins",
            ("0", "0", "0", "0"),
        )
        println(io, "ED_ground,$(@sprintf("%.15e", audit.ground)),,,,fixed_sector_ED,,,,")
        println(io, "ED_thermal_exact,$(@sprintf("%.15e", audit.thermal)),,,,grand_canonical_ED,,,,")
        println(io, "ED_trotter_exact,$(@sprintf("%.15e", audit.trotter)),,,,asymmetric_320_slice_ED,,,,")
        close(io)
        mv(temporary, path)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary; force = true)
    end
end

function main(args = ARGS)
    length(args) == 4 || error(
        "Usage: compare_u4_l4_untouched_baseline.jl " *
        "BASELINE_RAW OPTIMIZED_RAW AUDIT_CSV OUTPUT_CSV",
    )
    baseline_rows = read_forward_bins(args[1])
    optimized_rows = read_forward_bins(args[2])
    baseline = summarize_forward(baseline_rows)
    optimized = summarize_forward(optimized_rows)
    audit = read_audit(args[3])
    maxima = maximum_differences(baseline_rows, optimized_rows)

    abs(baseline.S2 - optimized.S2) <= 1e-10 ||
        error("Untouched and optimized S2 disagree")
    maximum(values(maxima)) <= 1e-9 ||
        error("Untouched and optimized bins disagree")

    write_comparison(args[4], baseline, optimized, audit, maxima)
    @printf("untouched S2 = %.12f +/- %.12f\n", baseline.S2, baseline.error)
    @printf("optimized S2 = %.12f +/- %.12f\n", optimized.S2, optimized.error)
    @printf("difference   = %+.3e\n", baseline.S2 - optimized.S2)
    @printf("max bin diff = %.3e\n", maximum(values(maxima)))
    println("Wrote comparison: ", abspath(args[4]))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

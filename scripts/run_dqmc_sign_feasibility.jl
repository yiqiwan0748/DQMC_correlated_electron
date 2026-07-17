using AltermagnetDQMCTrial
using LinearAlgebra: BLAS

BLAS.set_num_threads(1)

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTPUT_DIRECTORY = joinpath(
    ROOT,
    "results",
    "dqmc_sign_feasibility",
)
const SUMMARY_FILENAME = "sign_feasibility.csv"
const METADATA_FILENAME = "metadata.txt"

const TARGET_MU = -1.897
const TARGET_FILLING = 1.46021875
const ALL_PHASES = (
    :u0_control,
    :mu_calibration,
    :fixed_mu_diagnostics,
)

function feasibility_chain(;
    phase,
    case_name,
    L,
    beta,
    U,
    mu,
    seed,
    burnin_sweeps,
    measurement_sweeps,
    sweeps_between_measurements,
    block_size,
    target_filling=missing,
)
    return (
        phase=phase,
        case=case_name,
        L=L,
        beta=beta,
        dtau=0.1,
        U=U,
        mu=mu,
        stabilization_period=5,
        burnin_sweeps=burnin_sweeps,
        measurement_sweeps=measurement_sweeps,
        sweeps_between_measurements=sweeps_between_measurements,
        block_size=block_size,
        seed=seed,
        target_filling=target_filling,
    )
end

const FEASIBILITY_CHAINS = (
    feasibility_chain(
        phase=:u0_control,
        case_name="u0_target_L4_beta4",
        L=4, beta=4.0, U=0.0, mu=TARGET_MU, seed=6002,
        burnin_sweeps=100, measurement_sweeps=400,
        sweeps_between_measurements=1, block_size=50,
    ),
    feasibility_chain(
        phase=:mu_calibration,
        case_name="calibrate_target_mu_m2p497",
        L=4, beta=4.0, U=2.0, mu=-2.497, seed=6111,
        burnin_sweeps=200, measurement_sweeps=800,
        sweeps_between_measurements=2, block_size=100,
        target_filling=TARGET_FILLING,
    ),
    feasibility_chain(
        phase=:mu_calibration,
        case_name="calibrate_target_mu_m2p197",
        L=4, beta=4.0, U=2.0, mu=-2.197, seed=6112,
        burnin_sweeps=200, measurement_sweeps=800,
        sweeps_between_measurements=2, block_size=100,
        target_filling=TARGET_FILLING,
    ),
    feasibility_chain(
        phase=:mu_calibration,
        case_name="calibrate_target_mu_m1p897",
        L=4, beta=4.0, U=2.0, mu=TARGET_MU, seed=6113,
        burnin_sweeps=200, measurement_sweeps=800,
        sweeps_between_measurements=2, block_size=100,
        target_filling=TARGET_FILLING,
    ),
    feasibility_chain(
        phase=:mu_calibration,
        case_name="calibrate_target_mu_m1p597",
        L=4, beta=4.0, U=2.0, mu=-1.597, seed=6114,
        burnin_sweeps=200, measurement_sweeps=800,
        sweeps_between_measurements=2, block_size=100,
        target_filling=TARGET_FILLING,
    ),
    feasibility_chain(
        phase=:mu_calibration,
        case_name="calibrate_target_mu_m1p297",
        L=4, beta=4.0, U=2.0, mu=-1.297, seed=6115,
        burnin_sweeps=200, measurement_sweeps=800,
        sweeps_between_measurements=2, block_size=100,
        target_filling=TARGET_FILLING,
    ),
    feasibility_chain(
        phase=:fixed_mu_diagnostics,
        case_name="fixed_target_L4_beta4",
        L=4, beta=4.0, U=2.0, mu=TARGET_MU, seed=6211,
        burnin_sweeps=400, measurement_sweeps=1600,
        sweeps_between_measurements=2, block_size=100,
    ),
    feasibility_chain(
        phase=:fixed_mu_diagnostics,
        case_name="fixed_target_L4_beta4",
        L=4, beta=4.0, U=2.0, mu=TARGET_MU, seed=6212,
        burnin_sweeps=400, measurement_sweeps=1600,
        sweeps_between_measurements=2, block_size=100,
    ),
    feasibility_chain(
        phase=:fixed_mu_diagnostics,
        case_name="fixed_target_L6_beta4",
        L=6, beta=4.0, U=2.0, mu=TARGET_MU, seed=6231,
        burnin_sweeps=400, measurement_sweeps=1600,
        sweeps_between_measurements=2, block_size=100,
    ),
    feasibility_chain(
        phase=:fixed_mu_diagnostics,
        case_name="fixed_target_L6_beta4",
        L=6, beta=4.0, U=2.0, mu=TARGET_MU, seed=6232,
        burnin_sweeps=400, measurement_sweeps=1600,
        sweeps_between_measurements=2, block_size=100,
    ),
)

function phase_symbol(value::AbstractString)
    phase = Symbol(replace(lowercase(value), '-' => '_'))
    phase in ALL_PHASES || error(
        "unknown phase '$value'; choose u0-control, mu-calibration, or " *
        "fixed-mu-diagnostics",
    )
    return phase
end

function parse_runner_options(arguments)
    run_requested = false
    overwrite = false
    no_write = false
    help_requested = false
    selected_phases = Symbol[]
    selected_case = nothing
    selected_seed = nothing
    output_directory = DEFAULT_OUTPUT_DIRECTORY
    for argument in arguments
        if argument == "--run"
            run_requested = true
        elseif argument == "--overwrite"
            overwrite = true
        elseif argument == "--no-write"
            no_write = true
        elseif argument == "--help" || argument == "-h"
            help_requested = true
        elseif startswith(argument, "--phase=")
            push!(selected_phases, phase_symbol(split(argument, '='; limit=2)[2]))
        elseif startswith(argument, "--case=")
            selected_case = split(argument, '='; limit=2)[2]
        elseif startswith(argument, "--seed=")
            selected_seed = parse(Int, split(argument, '='; limit=2)[2])
        elseif startswith(argument, "--output-dir=")
            output_directory = abspath(split(argument, '='; limit=2)[2])
        else
            error("unknown option '$argument'; use --help for usage")
        end
    end
    overwrite && no_write && error("--overwrite and --no-write are incompatible")
    phases = isempty(selected_phases) ? collect(ALL_PHASES) : unique(selected_phases)
    return (
        run_requested=run_requested,
        overwrite=overwrite,
        no_write=no_write,
        help_requested=help_requested,
        phases=phases,
        selected_case=selected_case,
        selected_seed=selected_seed,
        output_directory=output_directory,
    )
end

function selected_chains(options)
    selected = [
        chain for chain in FEASIBILITY_CHAINS
        if chain.phase in options.phases &&
           (isnothing(options.selected_case) || chain.case == options.selected_case) &&
           (isnothing(options.selected_seed) || chain.seed == options.selected_seed)
    ]
    isempty(selected) && error("the phase/case/seed selection matched no chains")
    return selected
end

function print_usage(io::IO=stdout)
    println(io, "Usage: julia --project=. scripts/run_dqmc_sign_feasibility.jl [options]")
    println(io)
    println(io, "Default behavior is dry-run: print the full selected chain table only.")
    println(io, "  --run                         explicitly construct and run selected chains")
    println(io, "  --phase=PHASE                 repeatable phase selector")
    println(io, "  --case=NAME                   exact case-name selector")
    println(io, "  --seed=INTEGER                exact seed selector")
    println(io, "  --output-dir=PATH             output directory used only with --run")
    println(io, "  --no-write                    run but write no CSV or metadata")
    println(io, "  --overwrite                   allow replacement of existing output files")
    println(io, "  --help                        show this message")
end

table_value(value) = ismissing(value) ? "-" : string(value)

function print_chain_table(chains; io::IO=stdout)
    columns = (
        :phase, :case, :L, :beta, :dtau, :Ltau, :U, :mu, :burnin,
        :measurements, :between, :block_size, :seed, :target_filling,
    )
    println(io, join(string.(columns), '\t'))
    for chain in chains
        values = (
            chain.phase,
            chain.case,
            chain.L,
            chain.beta,
            chain.dtau,
            round(Int, chain.beta / chain.dtau),
            chain.U,
            chain.mu,
            chain.burnin_sweeps,
            chain.measurement_sweeps,
            chain.sweeps_between_measurements,
            chain.block_size,
            chain.seed,
            chain.target_filling,
        )
        println(io, join(table_value.(values), '\t'))
    end
end

function csv_value(value)
    ismissing(value) && return ""
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return '"' * replace(text, '"' => "\"\"") * '"'
    end
    return text
end

function write_summary(path, summaries)
    columns = propertynames(first(summaries))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for summary in summaries
            println(io, join(
                (csv_value(getproperty(summary, column)) for column in columns),
                ',',
            ))
        end
    end
end

function write_metadata(path, chains, criteria)
    open(path, "w") do io
        println(io, "Milestone 6A DQMC sign-feasibility diagnostics")
        println(io)
        println(io, "No average-sign conclusion is encoded in this configuration file.")
        println(io, "Each CSV row must be interpreted through all diagnostic flags.")
        println(io, "average_sign_gate=average_sign - 2*SE > threshold (strict)")
        println(io, "average_sign_threshold=$(criteria.average_sign_threshold)")
        println(io, "maximum_average_sign_standard_error=$(criteria.maximum_average_sign_standard_error)")
        println(io, "minimum_block_count=$(criteria.minimum_block_count)")
        println(io, "maximum_block_drift_zscore=$(criteria.maximum_block_drift_zscore)")
        println(io, "minimum_sign_signal_to_noise=$(criteria.minimum_sign_signal_to_noise)")
        println(io, "minimum_absolute_average_sign=$(criteria.minimum_absolute_average_sign)")
        println(io, "maximum_greens_recompute_error=$(criteria.maximum_greens_recompute_error)")
        println(io, "u0_filling_tolerance=$(criteria.u0_filling_tolerance)")
        println(io)
        println(io, "Selected chain table")
        println(io, "--------------------")
        print_chain_table(chains; io=io)
    end
end

function check_output_collisions(options)
    options.no_write && return nothing
    summary_path = joinpath(options.output_directory, SUMMARY_FILENAME)
    metadata_path = joinpath(options.output_directory, METADATA_FILENAME)
    collisions = filter(isfile, (summary_path, metadata_path))
    isempty(collisions) || options.overwrite || error(
        "refusing to overwrite existing output(s): $(join(collisions, ", ")); " *
        "pass --overwrite explicitly",
    )
    return nothing
end

function run_selected_chains(chains, criteria)
    return map(chains) do chain
        parameters = DQMCParameters(
            L=chain.L,
            beta=chain.beta,
            dtau=chain.dtau,
            U=chain.U,
            mu=chain.mu,
            stabilization_period=chain.stabilization_period,
            burnin_sweeps=chain.burnin_sweeps,
            measurement_sweeps=chain.measurement_sweeps,
            sweeps_between_measurements=chain.sweeps_between_measurements,
            seed=chain.seed,
        )
        summary = run_dqmc_sign_feasibility(
            parameters;
            case_name=chain.case,
            phase=chain.phase,
            block_size=chain.block_size,
            criteria=criteria,
        )
        println(
            "CHAIN_COMPLETE phase=$(chain.phase) case=$(chain.case) " *
            "seed=$(chain.seed) runtime_seconds=$(summary.runtime_seconds)",
        )
        return merge((target_filling=chain.target_filling,), summary)
    end
end

function main(arguments=ARGS)
    options = parse_runner_options(arguments)
    if options.help_requested
        print_usage()
        return nothing
    end
    chains = selected_chains(options)
    println("mode=", options.run_requested ? "run" : "dry-run")
    println("selected_chain_count=", length(chains))
    print_chain_table(chains)

    if !options.run_requested
        println("dqmc_chains_constructed=0")
        println("outputs_written=false")
        return nothing
    end

    check_output_collisions(options)
    criteria = DEFAULT_SIGN_FEASIBILITY_CRITERIA
    summaries = run_selected_chains(chains, criteria)
    if options.no_write
        println("outputs_written=false")
        return nothing
    end
    mkpath(options.output_directory)
    summary_path = joinpath(options.output_directory, SUMMARY_FILENAME)
    metadata_path = joinpath(options.output_directory, METADATA_FILENAME)
    write_summary(summary_path, summaries)
    write_metadata(metadata_path, chains, criteria)
    println("Wrote ", summary_path)
    println("Wrote ", metadata_path)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

using AltermagnetDQMCTrial
using LinearAlgebra: BLAS

BLAS.set_num_threads(1)

const CONFIRMED_FIXED_DENSITY_MU = -2.5202
const TARGET_FILLING = 1.46021875
const WARMUP_SWEEPS = 1
const TIMED_SWEEPS = 10
const PROJECTED_COMPLETE_SWEEPS_PER_CHAIN = 3600
const PROJECTED_CHAIN_COUNT = 2

function resource_benchmark_case(;
    case_name,
    L,
    seed,
    execution_authorized=true,
)
    return (
        case=case_name,
        L=L,
        beta=4.0,
        dtau=0.1,
        U=2.0,
        mu=CONFIRMED_FIXED_DENSITY_MU,
        stabilization_period=5,
        seed=seed,
        target_filling=TARGET_FILLING,
        warmup_sweeps=WARMUP_SWEEPS,
        timed_sweeps=TIMED_SWEEPS,
        execution_authorized=execution_authorized,
        projected_complete_sweeps_per_chain=
            PROJECTED_COMPLETE_SWEEPS_PER_CHAIN,
        projected_chain_count=PROJECTED_CHAIN_COUNT,
    )
end

const RESOURCE_BENCHMARKS = (
    resource_benchmark_case(
        case_name="resource_L4_beta4",
        L=4,
        seed=6704,
    ),
    resource_benchmark_case(
        case_name="resource_L6_beta4",
        L=6,
        seed=6706,
    ),
    resource_benchmark_case(
        case_name="resource_L8_beta4",
        L=8,
        seed=6708,
    ),
    resource_benchmark_case(
        case_name="resource_L10_beta4",
        L=10,
        seed=6710,
    ),
    resource_benchmark_case(
        case_name="resource_L12_beta4",
        L=12,
        seed=6712,
    ),
    resource_benchmark_case(
        case_name="resource_L14_beta4",
        L=14,
        seed=6714,
    ),
    resource_benchmark_case(
        case_name="resource_L16_beta4",
        L=16,
        seed=6716,
        execution_authorized=true,
    ),
)

function parse_benchmark_options(arguments)
    run_requested = false
    require_single = false
    help_requested = false
    selected_case = nothing
    for argument in arguments
        if argument == "--run"
            run_requested = true
        elseif argument == "--require-single"
            require_single = true
        elseif argument == "--help" || argument == "-h"
            help_requested = true
        elseif startswith(argument, "--case=")
            selected_case = split(argument, '='; limit=2)[2]
        else
            error(string("unknown option '", argument, "'; use --help for usage"))
        end
    end
    return (
        run_requested=run_requested,
        require_single=require_single,
        help_requested=help_requested,
        selected_case=selected_case,
    )
end

function selected_benchmarks(options)
    selected = [
        case for case in RESOURCE_BENCHMARKS
        if isnothing(options.selected_case) || case.case == options.selected_case
    ]
    isempty(selected) && error("the case selector matched no resource benchmark")
    return selected
end

function validate_benchmark_selection(options, cases)
    if options.require_single
        single_case_guard_passed = length(cases) == 1
        println("single_case_guard_passed=", single_case_guard_passed)
        single_case_guard_passed || error(
            string(
                "--require-single matched ",
                length(cases),
                " resource benchmarks; expected exactly 1",
            ),
        )
    end
    if options.run_requested
        options.require_single || error("--run requires --require-single")
        isnothing(options.selected_case) &&
            error("--run requires an exact --case selector")
        only(cases).execution_authorized || error(
            string(
                "real execution of ",
                only(cases).case,
                " is not authorized by the current resource registration",
            ),
        )
    end
    return nothing
end

function print_usage(io::IO=stdout)
    println(
        io,
        "Usage: julia --project=. scripts/benchmark_dqmc_scaling.jl [options]",
    )
    println(io)
    println(io, "Default behavior is dry-run: print registered definitions only.")
    println(io, "  --run               execute one selected resource benchmark")
    println(io, "  --case=NAME         exact registered case selector")
    println(io, "  --require-single    require exactly one selected case")
    println(io, "  --help              show this message")
    println(io)
    println(io, "Real execution requires --run, --require-single, and --case.")
    println(
        io,
        "Capture peak RSS externally with /usr/bin/time -v around the Julia command.",
    )
end

function print_benchmark_table(cases; io::IO=stdout)
    columns = (
        :case,
        :L,
        :beta,
        :dtau,
        :Ltau,
        :U,
        :mu,
        :seed,
        :warmup_sweeps,
        :timed_complete_sweeps,
        :execution_authorized,
        :projected_complete_sweeps_per_chain,
        :projected_chain_count,
        :target_filling,
    )
    println(io, join(string.(columns), '\t'))
    for case in cases
        values = (
            case.case,
            case.L,
            case.beta,
            case.dtau,
            round(Int, case.beta / case.dtau),
            case.U,
            case.mu,
            case.seed,
            case.warmup_sweeps,
            case.timed_sweeps,
            case.execution_authorized,
            case.projected_complete_sweeps_per_chain,
            case.projected_chain_count,
            case.target_filling,
        )
        println(io, join(string.(values), '\t'))
    end
    return nothing
end

function projected_two_chain_seconds(
    elapsed_seconds::Real,
    timed_sweeps::Integer;
    projected_complete_sweeps_per_chain::Integer=
        PROJECTED_COMPLETE_SWEEPS_PER_CHAIN,
    projected_chain_count::Integer=PROJECTED_CHAIN_COUNT,
)
    isfinite(elapsed_seconds) && elapsed_seconds >= 0 ||
        throw(ArgumentError("elapsed seconds must be finite and nonnegative"))
    timed_sweeps > 0 ||
        throw(ArgumentError("timed sweeps must be positive"))
    projected_complete_sweeps_per_chain > 0 || throw(
        ArgumentError("projected complete sweeps per chain must be positive"),
    )
    projected_chain_count > 0 ||
        throw(ArgumentError("projected chain count must be positive"))
    seconds_per_sweep = elapsed_seconds / timed_sweeps
    return seconds_per_sweep *
           projected_complete_sweeps_per_chain *
           projected_chain_count
end

function run_resource_benchmark(case)
    parameters = DQMCParameters(
        L=case.L,
        beta=case.beta,
        dtau=case.dtau,
        U=case.U,
        mu=case.mu,
        stabilization_period=case.stabilization_period,
        burnin_sweeps=0,
        measurement_sweeps=1,
        sweeps_between_measurements=1,
        seed=case.seed,
    )
    state = initialize_dqmc_state(parameters)
    for _ in 1:case.warmup_sweeps
        sweep_all_time_slices!(state)
    end

    attempted_updates = 0
    accepted_updates = 0
    start_time = time_ns()
    for _ in 1:case.timed_sweeps
        sweep = sweep_all_time_slices!(state)
        attempted_updates += sweep.attempted
        accepted_updates += sweep.accepted
    end
    elapsed_seconds = (time_ns() - start_time) / 1.0e9
    seconds_per_complete_sweep = elapsed_seconds / case.timed_sweeps
    projected_seconds = projected_two_chain_seconds(
        elapsed_seconds,
        case.timed_sweeps;
        projected_complete_sweeps_per_chain=
            case.projected_complete_sweeps_per_chain,
        projected_chain_count=case.projected_chain_count,
    )
    return (
        case=case.case,
        L=case.L,
        elapsed_seconds=elapsed_seconds,
        seconds_per_complete_sweep=seconds_per_complete_sweep,
        projected_two_chain_seconds=projected_seconds,
        projected_two_chain_hours=projected_seconds / 3600,
        attempted_updates=attempted_updates,
        accepted_updates=accepted_updates,
        acceptance_rate=accepted_updates / attempted_updates,
        maximum_greens_recompute_error=
            state.maximum_greens_recompute_error,
    )
end

function print_benchmark_result(result; io::IO=stdout)
    println(io, "benchmark_case=", result.case)
    println(io, "L=", result.L)
    println(io, "timed_sweep_seconds=", result.elapsed_seconds)
    println(
        io,
        "seconds_per_complete_sweep=",
        result.seconds_per_complete_sweep,
    )
    println(
        io,
        "projected_complete_sweeps_per_chain=",
        PROJECTED_COMPLETE_SWEEPS_PER_CHAIN,
    )
    println(io, "projected_chain_count=", PROJECTED_CHAIN_COUNT)
    println(
        io,
        "projected_two_chain_seconds=",
        result.projected_two_chain_seconds,
    )
    println(io, "projected_two_chain_hours=", result.projected_two_chain_hours)
    println(io, "attempted_updates=", result.attempted_updates)
    println(io, "accepted_updates=", result.accepted_updates)
    println(io, "acceptance_rate=", result.acceptance_rate)
    println(
        io,
        "maximum_greens_recompute_error=",
        result.maximum_greens_recompute_error,
    )
    println(io, "peak_rss_source=external_/usr/bin/time_-v")
    println(io, "peak_rss_captured_by_script=false")
    return nothing
end

function main(arguments=ARGS)
    options = parse_benchmark_options(arguments)
    if options.help_requested
        print_usage()
        return nothing
    end

    cases = selected_benchmarks(options)
    println("RESOURCE_BENCHMARK_ONLY=true")
    println("PHYSICS_CLAIM_AUTHORIZED=false")
    println("no_physics_claim_reason=short_timing_run_without_equilibration")
    println("mode=", options.run_requested ? "run" : "dry-run")
    println("selected_benchmark_count=", length(cases))
    print_benchmark_table(cases)
    validate_benchmark_selection(options, cases)

    if !options.run_requested
        println("dqmc_states_constructed=0")
        println("benchmark_sweeps_executed=0")
        println("outputs_written=false")
        return nothing
    end

    result = run_resource_benchmark(only(cases))
    print_benchmark_result(result)
    println("outputs_written=false")
    return result
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

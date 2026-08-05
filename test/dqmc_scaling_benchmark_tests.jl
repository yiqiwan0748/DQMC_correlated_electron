using Test

module DQMCScalingBenchmark
include(joinpath(@__DIR__, "..", "scripts", "benchmark_dqmc_scaling.jl"))
end

@testset "DQMC resource benchmark registration and safety" begin
    benchmarks = DQMCScalingBenchmark.RESOURCE_BENCHMARKS

    @test DQMCScalingBenchmark.CONFIRMED_FIXED_DENSITY_MU == -2.5202
    @test DQMCScalingBenchmark.TARGET_FILLING == 1.46021875
    @test DQMCScalingBenchmark.WARMUP_SWEEPS == 1
    @test DQMCScalingBenchmark.TIMED_SWEEPS == 10
    @test DQMCScalingBenchmark.PROJECTED_COMPLETE_SWEEPS_PER_CHAIN == 3600
    @test DQMCScalingBenchmark.PROJECTED_CHAIN_COUNT == 2
    @test length(benchmarks) == 7
    @test [case.case for case in benchmarks] == [
        "resource_L4_beta4",
        "resource_L6_beta4",
        "resource_L8_beta4",
        "resource_L10_beta4",
        "resource_L12_beta4",
        "resource_L14_beta4",
        "resource_L16_beta4",
    ]
    @test [case.L for case in benchmarks] == [4, 6, 8, 10, 12, 14, 16]
    @test [case.seed for case in benchmarks] ==
          [6704, 6706, 6708, 6710, 6712, 6714, 6716]
    @test all(case.execution_authorized for case in benchmarks)

    for case in benchmarks
        @test case.beta == 4.0
        @test case.dtau == 0.1
        @test round(Int, case.beta / case.dtau) == 40
        @test case.U == 2.0
        @test case.mu == -2.5202
        @test case.stabilization_period == 5
        @test case.target_filling == 1.46021875
        @test case.warmup_sweeps == 1
        @test case.timed_sweeps == 10
        @test case.projected_complete_sweeps_per_chain == 3600
        @test case.projected_chain_count == 2
    end

    default_options =
        DQMCScalingBenchmark.parse_benchmark_options(String[])
    @test !default_options.run_requested
    @test !default_options.require_single
    @test isnothing(default_options.selected_case)

    dry_run_pipe = Pipe()
    redirect_stdout(dry_run_pipe) do
        @test isnothing(DQMCScalingBenchmark.main(String[]))
    end
    close(dry_run_pipe.in)
    dry_run_output = read(dry_run_pipe, String)
    @test occursin("RESOURCE_BENCHMARK_ONLY=true", dry_run_output)
    @test occursin("PHYSICS_CLAIM_AUTHORIZED=false", dry_run_output)
    @test occursin(
        "no_physics_claim_reason=short_timing_run_without_equilibration",
        dry_run_output,
    )
    @test occursin("mode=dry-run", dry_run_output)
    @test occursin("selected_benchmark_count=7", dry_run_output)
    @test occursin("dqmc_states_constructed=0", dry_run_output)
    @test occursin("benchmark_sweeps_executed=0", dry_run_output)
    @test occursin("outputs_written=false", dry_run_output)

    selected_arguments = [
        "--case=resource_L8_beta4",
        "--require-single",
    ]
    selected_options =
        DQMCScalingBenchmark.parse_benchmark_options(selected_arguments)
    selected =
        DQMCScalingBenchmark.selected_benchmarks(selected_options)
    @test length(selected) == 1
    @test only(selected).L == 8
    selected_pipe = Pipe()
    redirect_stdout(selected_pipe) do
        @test isnothing(DQMCScalingBenchmark.main(selected_arguments))
    end
    close(selected_pipe.in)
    selected_output = read(selected_pipe, String)
    @test occursin("single_case_guard_passed=true", selected_output)
    @test occursin("dqmc_states_constructed=0", selected_output)
    @test occursin("benchmark_sweeps_executed=0", selected_output)

    runnable_options = DQMCScalingBenchmark.parse_benchmark_options([
        "--run",
        "--case=resource_L4_beta4",
        "--require-single",
    ])
    runnable =
        DQMCScalingBenchmark.selected_benchmarks(runnable_options)
    redirect_stdout(devnull) do
        @test isnothing(
            DQMCScalingBenchmark.validate_benchmark_selection(
                runnable_options,
                runnable,
            ),
        )
    end

    unguarded_options = DQMCScalingBenchmark.parse_benchmark_options([
        "--run",
        "--case=resource_L4_beta4",
    ])
    unguarded =
        DQMCScalingBenchmark.selected_benchmarks(unguarded_options)
    @test_throws ErrorException DQMCScalingBenchmark.validate_benchmark_selection(
        unguarded_options,
        unguarded,
    )

    unselected_run_options =
        DQMCScalingBenchmark.parse_benchmark_options([
            "--run",
            "--require-single",
        ])
    redirect_stdout(devnull) do
        @test_throws ErrorException DQMCScalingBenchmark.validate_benchmark_selection(
            unselected_run_options,
            collect(benchmarks),
        )
    end

    extension_options = DQMCScalingBenchmark.parse_benchmark_options([
        "--case=resource_L16_beta4",
        "--require-single",
    ])
    extension = DQMCScalingBenchmark.selected_benchmarks(extension_options)
    @test length(extension) == 1
    @test only(extension).L == 16
    @test only(extension).seed == 6716
    @test only(extension).execution_authorized
    extension_run =
        DQMCScalingBenchmark.parse_benchmark_options([
            "--run",
            "--case=resource_L16_beta4",
            "--require-single",
        ])
    redirect_stdout(devnull) do
        @test isnothing(
            DQMCScalingBenchmark.validate_benchmark_selection(
                extension_run,
                extension,
            ),
        )
    end

    missing_options = DQMCScalingBenchmark.parse_benchmark_options([
        "--case=resource_L18_beta4",
    ])
    @test_throws ErrorException DQMCScalingBenchmark.selected_benchmarks(
        missing_options,
    )
    @test_throws ErrorException DQMCScalingBenchmark.parse_benchmark_options([
        "--output-dir=results",
    ])

    @test DQMCScalingBenchmark.projected_two_chain_seconds(10.0, 10) ==
          7200.0
    @test DQMCScalingBenchmark.projected_two_chain_seconds(
        2.5,
        10;
        projected_complete_sweeps_per_chain=100,
        projected_chain_count=3,
    ) == 75.0
    @test_throws ArgumentError DQMCScalingBenchmark.projected_two_chain_seconds(
        -1.0,
        10,
    )
    @test_throws ArgumentError DQMCScalingBenchmark.projected_two_chain_seconds(
        1.0,
        0,
    )

    usage = IOBuffer()
    DQMCScalingBenchmark.print_usage(usage)
    usage_text = String(take!(usage))
    @test occursin("/usr/bin/time -v", usage_text)
    @test occursin("Real execution requires", usage_text)

    result_output = IOBuffer()
    DQMCScalingBenchmark.print_benchmark_result(
        (
            case="resource_L4_beta4",
            L=4,
            elapsed_seconds=1.0,
            seconds_per_complete_sweep=0.1,
            projected_two_chain_seconds=720.0,
            projected_two_chain_hours=0.2,
            attempted_updates=1280,
            accepted_updates=640,
            acceptance_rate=0.5,
            maximum_greens_recompute_error=1.0e-12,
        );
        io=result_output,
    )
    result_text = String(take!(result_output))
    @test occursin("projected_complete_sweeps_per_chain=3600", result_text)
    @test occursin("projected_chain_count=2", result_text)
    @test occursin("projected_two_chain_seconds=720.0", result_text)
    @test occursin("maximum_greens_recompute_error=1.0e-12", result_text)
    @test occursin("peak_rss_source=external_/usr/bin/time_-v", result_text)
    @test occursin("peak_rss_captured_by_script=false", result_text)

    source = read(
        joinpath(@__DIR__, "..", "scripts", "benchmark_dqmc_scaling.jl"),
        String,
    )
    @test !occursin("open(", source)
    @test !occursin("mkpath", source)
end

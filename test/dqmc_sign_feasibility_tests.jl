using AltermagnetDQMCTrial
using Test

module SignFeasibilityRunner
include(joinpath(@__DIR__, "..", "scripts", "run_dqmc_sign_feasibility.jl"))
end

function _synthetic_thermodynamic_observation(
    density,
    kinetic_energy_per_cell,
    shifted_interaction_energy_per_cell,
    mu,
)
    internal_energy_per_cell =
        kinetic_energy_per_cell + shifted_interaction_energy_per_cell
    return (
        density=density,
        doping=(density - 2) / 4,
        kinetic_energy_per_cell=kinetic_energy_per_cell,
        shifted_interaction_energy_per_cell=
            shifted_interaction_energy_per_cell,
        internal_energy_per_cell=internal_energy_per_cell,
        grand_energy_per_cell=internal_energy_per_cell - mu * density,
    )
end

function _independent_ratio_jackknife(signs, values, block_size)
    measurement_count = length(signs)
    ranges = [
        first_index:min(first_index + block_size - 1, measurement_count)
        for first_index in 1:block_size:measurement_count
    ]
    sign_sum = sum(signs)
    signed_sum = sum(signs[index] * values[index] for index in eachindex(signs))
    leave_one_out = [
        begin
            block_sign_sum = sum(signs[index] for index in range)
            block_signed_sum = sum(
                signs[index] * values[index] for index in range
            )
            (signed_sum - block_signed_sum) / (sign_sum - block_sign_sum)
        end
        for range in ranges
    ]
    center = sum(leave_one_out) / length(leave_one_out)
    standard_error = sqrt(
        (length(leave_one_out) - 1) / length(leave_one_out) *
        sum((estimate - center)^2 for estimate in leave_one_out),
    )
    return (estimate=signed_sum / sign_sum, standard_error=standard_error)
end

function _strict_pilot_summary(;
    signs=ones(16),
    raw_total_densities=collect(range(1.1, 1.8; length=length(signs))),
    block_size=2,
    maximum_greens_recompute_error=0.0,
    acceptance_rate=0.5,
    attempted_updates=100,
    accepted_updates=50,
    rejected_updates=50,
    runtime_seconds=1.0,
)
    diagnostics = summarize_sign_feasibility(
        signs,
        raw_total_densities;
        block_size=block_size,
        maximum_greens_recompute_error=maximum_greens_recompute_error,
        acceptance_rate=acceptance_rate,
    )
    return merge(
        diagnostics,
        (
            attempted_updates=attempted_updates,
            accepted_updates=accepted_updates,
            rejected_updates=rejected_updates,
            runtime_seconds=runtime_seconds,
        ),
    )
end

function _repeated_sign_blocks(positive_count, block_length, block_count)
    block = vcat(
        ones(positive_count),
        -ones(block_length - positive_count),
    )
    return repeat(block, block_count)
end

@testset "Milestone 6A DQMC sign feasibility" begin
    default_criteria = SignFeasibilityCriteria()
    criteria = SignFeasibilityCriteria(
        average_sign_threshold=0.2,
        maximum_average_sign_standard_error=0.1,
        minimum_block_count=4,
        maximum_block_drift_zscore=2.0,
        minimum_sign_signal_to_noise=3.0,
        minimum_absolute_average_sign=1.0e-8,
        maximum_greens_recompute_error=1.0e-6,
        u0_filling_tolerance=1.0e-10,
    )

    @testset "Criteria and strict average-sign gate" begin
        @test default_criteria.average_sign_threshold ==
              DEFAULT_SIGN_FEASIBILITY_CRITERIA.average_sign_threshold == 0.2
        @test default_criteria.maximum_average_sign_standard_error ==
              DEFAULT_SIGN_FEASIBILITY_CRITERIA.maximum_average_sign_standard_error ==
              0.02
        @test default_criteria.minimum_sign_signal_to_noise ==
              DEFAULT_SIGN_FEASIBILITY_CRITERIA.minimum_sign_signal_to_noise == 5.0
        @test default_criteria.maximum_greens_recompute_error ==
              DEFAULT_SIGN_FEASIBILITY_CRITERIA.maximum_greens_recompute_error ==
              1.0e-8

        @test !average_sign_gate_passes(2 // 5, 1 // 10, 1 // 5)
        @test average_sign_gate_passes(401 // 1000, 1 // 10, 1 // 5)
        @test !average_sign_gate_passes(0.8, Inf, 0.1)
        @test average_sign_gate_passes(0.8, -0.1, 0.1)

        @test_throws ArgumentError SignFeasibilityCriteria(
            average_sign_threshold=1.0,
        )
        @test_throws ArgumentError SignFeasibilityCriteria(
            minimum_block_count=0,
        )
        @test_throws ArgumentError SignFeasibilityCriteria(
            minimum_sign_signal_to_noise=-1.0,
        )
        @test_throws ArgumentError SignFeasibilityCriteria(
            maximum_greens_recompute_error=NaN,
        )
    end

    @testset "Target-only runner configuration and dry-run" begin
        chains = SignFeasibilityRunner.FEASIBILITY_CHAINS
        u0_controls = filter(chain -> chain.phase == :u0_control, chains)
        calibration = filter(chain -> chain.phase == :mu_calibration, chains)
        fixed_mu = filter(chain -> chain.phase == :fixed_mu_diagnostics, chains)

        @test SignFeasibilityRunner.TARGET_MU == -1.897
        @test length(chains) == 14
        @test length(u0_controls) == 1
        @test length(calibration) == 9
        @test only(u0_controls).U == 0.0
        @test only(u0_controls).mu == -1.897
        @test [chain.mu for chain in calibration] ==
              [-2.55, -2.5202, -2.5202, -2.5202, -2.497, -2.197,
               -1.897, -1.597, -1.297]
        lower_bracket = only(filter(
            chain -> chain.case == "calibrate_target_mu_m2p55",
            calibration,
        ))
        @test lower_bracket.phase == :mu_calibration
        @test lower_bracket.case == "calibrate_target_mu_m2p55"
        @test lower_bracket.L == 4
        @test lower_bracket.beta == 4.0
        @test lower_bracket.dtau == 0.1
        @test round(Int, lower_bracket.beta / lower_bracket.dtau) == 40
        @test lower_bracket.U == 2.0
        @test lower_bracket.mu == -2.55
        @test lower_bracket.seed == 6116
        @test lower_bracket.burnin_sweeps == 200
        @test lower_bracket.measurement_sweeps == 800
        @test lower_bracket.sweeps_between_measurements == 2
        @test lower_bracket.block_size == 100
        @test lower_bracket.stabilization_period == 5
        @test lower_bracket.target_filling == 1.46021875

        lower_bracket_options = SignFeasibilityRunner.parse_runner_options([
            "--phase=mu-calibration",
            "--case=calibrate_target_mu_m2p55",
            "--seed=6116",
            "--require-single",
            "--strict-pilot",
        ])
        lower_bracket_selection =
            SignFeasibilityRunner.selected_chains(lower_bracket_options)
        @test length(lower_bracket_selection) == 1
        @test only(lower_bracket_selection) == lower_bracket

        candidate_chains = filter(
            chain -> chain.case == "calibrate_target_mu_m2p5202",
            calibration,
        )
        @test [chain.seed for chain in candidate_chains] == [6117, 6118, 6119]
        candidate = only(filter(chain -> chain.seed == 6117, candidate_chains))
        @test candidate.phase == :mu_calibration
        @test candidate.case == "calibrate_target_mu_m2p5202"
        @test candidate.L == 4
        @test candidate.beta == 4.0
        @test candidate.dtau == 0.1
        @test round(Int, candidate.beta / candidate.dtau) == 40
        @test candidate.U == 2.0
        @test candidate.mu == -2.5202
        @test candidate.seed == 6117
        @test candidate.burnin_sweeps == 200
        @test candidate.measurement_sweeps == 800
        @test candidate.sweeps_between_measurements == 2
        @test candidate.block_size == 100
        @test candidate.stabilization_period == 5
        @test candidate.target_filling == 1.46021875

        candidate_options = SignFeasibilityRunner.parse_runner_options([
            "--phase=mu-calibration",
            "--case=calibrate_target_mu_m2p5202",
            "--seed=6117",
            "--require-single",
            "--strict-pilot",
        ])
        candidate_selection =
            SignFeasibilityRunner.selected_chains(candidate_options)
        @test candidate_options.phases == [:mu_calibration]
        @test candidate_options.selected_case ==
              "calibrate_target_mu_m2p5202"
        @test candidate_options.selected_seed == 6117
        @test candidate_options.require_single
        @test candidate_options.strict_pilot
        @test length(candidate_selection) == 1
        @test only(candidate_selection) == candidate

        candidate_prospective_metadata = (
            result_identifier=
                "$(candidate.case)_seed$(candidate.seed)_strict_pilot",
            phase=candidate.phase,
            case=candidate.case,
            seed=candidate.seed,
            mu=candidate.mu,
            target_filling=candidate.target_filling,
            strict_pilot=candidate_options.strict_pilot,
        )
        @test candidate_prospective_metadata == (
            result_identifier=
                "calibrate_target_mu_m2p5202_seed6117_strict_pilot",
            phase=:mu_calibration,
            case="calibrate_target_mu_m2p5202",
            seed=6117,
            mu=-2.5202,
            target_filling=1.46021875,
            strict_pilot=true,
        )

        candidate_case_options = SignFeasibilityRunner.parse_runner_options([
            "--phase=mu-calibration",
            "--case=calibrate_target_mu_m2p5202",
            "--require-single",
            "--strict-pilot",
        ])
        candidate_case_selection =
            SignFeasibilityRunner.selected_chains(candidate_case_options)
        @test candidate_case_selection == collect(candidate_chains)
        @test length(candidate_case_selection) == 3
        redirect_stdout(devnull) do
            @test_throws ErrorException SignFeasibilityRunner.validate_chain_selection(
                candidate_case_options,
                candidate_case_selection,
            )
        end

        for (seed, result_identifier) in (
            (6118, "calibrate_target_mu_m2p5202_seed6118_strict_pilot"),
            (6119, "calibrate_target_mu_m2p5202_seed6119_strict_pilot"),
        )
            confirmation = only(filter(chain -> chain.seed == seed, candidate_chains))
            @test confirmation.phase == :mu_calibration
            @test confirmation.case == "calibrate_target_mu_m2p5202"
            @test confirmation.L == 4
            @test confirmation.beta == 4.0
            @test confirmation.dtau == 0.1
            @test round(Int, confirmation.beta / confirmation.dtau) == 40
            @test confirmation.U == 2.0
            @test confirmation.mu == -2.5202
            @test confirmation.seed == seed
            @test confirmation.burnin_sweeps == 200
            @test confirmation.measurement_sweeps == 800
            @test confirmation.sweeps_between_measurements == 2
            @test confirmation.block_size == 100
            @test confirmation.stabilization_period == 5
            @test confirmation.target_filling == 1.46021875

            confirmation_options = SignFeasibilityRunner.parse_runner_options([
                "--phase=mu-calibration",
                "--case=calibrate_target_mu_m2p5202",
                "--seed=$seed",
                "--require-single",
                "--strict-pilot",
            ])
            confirmation_selection =
                SignFeasibilityRunner.selected_chains(confirmation_options)
            @test confirmation_options.phases == [:mu_calibration]
            @test confirmation_options.selected_case ==
                  "calibrate_target_mu_m2p5202"
            @test confirmation_options.selected_seed == seed
            @test confirmation_options.require_single
            @test confirmation_options.strict_pilot
            @test length(confirmation_selection) == 1
            @test only(confirmation_selection) == confirmation

            confirmation_prospective_metadata = (
                result_identifier=
                    "$(confirmation.case)_seed$(confirmation.seed)_strict_pilot",
                phase=confirmation.phase,
                case=confirmation.case,
                seed=confirmation.seed,
                mu=confirmation.mu,
                target_filling=confirmation.target_filling,
                strict_pilot=confirmation_options.strict_pilot,
            )
            @test confirmation_prospective_metadata == (
                result_identifier=result_identifier,
                phase=:mu_calibration,
                case="calibrate_target_mu_m2p5202",
                seed=seed,
                mu=-2.5202,
                target_filling=1.46021875,
                strict_pilot=true,
            )
        end
        @test !isempty(fixed_mu)
        @test all(chain -> chain.mu == -1.897, fixed_mu)
        @test all(chain -> chain.mu == -1.897,
                  filter(chain -> chain.phase != :mu_calibration, chains))
        @test all(chain -> !occursin("lower_energy_vhs", chain.case), chains)
        @test all(chain -> chain.mu != -3.8369423191, chains)

        options = SignFeasibilityRunner.parse_runner_options(String[])
        @test !options.run_requested
        @test !options.require_single
        @test !options.strict_pilot
        pipe = Pipe()
        redirect_stdout(pipe) do
            @test isnothing(SignFeasibilityRunner.main(String[]))
        end
        close(pipe.in)
        output = read(pipe, String)
        @test occursin("mode=dry-run", output)
        @test occursin("dqmc_chains_constructed=0", output)
        @test occursin("outputs_written=false", output)

        one_chain_arguments = [
            "--phase=fixed-mu-diagnostics",
            "--case=fixed_target_L4_beta4",
            "--seed=6211",
            "--require-single",
            "--strict-pilot",
        ]
        one_chain_options = SignFeasibilityRunner.parse_runner_options(
            one_chain_arguments,
        )
        selected = SignFeasibilityRunner.selected_chains(one_chain_options)
        @test length(selected) == 1
        @test only(selected).seed == 6211
        one_chain_pipe = Pipe()
        redirect_stdout(one_chain_pipe) do
            @test isnothing(SignFeasibilityRunner.main(one_chain_arguments))
        end
        close(one_chain_pipe.in)
        one_chain_output = read(one_chain_pipe, String)
        @test occursin("mode=dry-run", one_chain_output)
        @test occursin("selected_chain_count=1", one_chain_output)
        @test occursin("single_chain_guard_passed=true", one_chain_output)
        @test occursin("dqmc_chains_constructed=0", one_chain_output)
        @test occursin("outputs_written=false", one_chain_output)

        run_one_options = SignFeasibilityRunner.parse_runner_options(
            [one_chain_arguments; "--run"],
        )
        redirect_stdout(devnull) do
            @test isnothing(SignFeasibilityRunner.validate_chain_selection(
                run_one_options,
                selected,
            ))
        end

        multiple_options = SignFeasibilityRunner.parse_runner_options(
            ["--run", "--require-single"],
        )
        redirect_stdout(devnull) do
            @test_throws ErrorException SignFeasibilityRunner.validate_chain_selection(
                multiple_options,
                chains,
            )
        end

        zero_options = SignFeasibilityRunner.parse_runner_options([
            "--run",
            "--require-single",
            "--case=not_a_configured_case",
        ])
        @test_throws ErrorException SignFeasibilityRunner.selected_chains(
            zero_options,
        )

        unguarded_strict_options =
            SignFeasibilityRunner.parse_runner_options([
                "--strict-pilot",
            ])
        @test_throws ErrorException SignFeasibilityRunner.validate_chain_selection(
            unguarded_strict_options,
            chains,
        )
    end

    @testset "Runner diagnostics output and overwrite refusal" begin
        chain = SignFeasibilityRunner.FEASIBILITY_CHAINS[7]
        diagnostics = _strict_pilot_summary()
        assessment = AltermagnetDQMCTrial.classify_strict_sign_pilot(
            diagnostics,
        )
        summary = merge(
            (
                target_filling=chain.target_filling,
                phase=String(chain.phase),
                case=chain.case,
                seed=chain.seed,
            ),
            diagnostics,
            assessment,
        )
        base_summary = merge(
            (
                target_filling=chain.target_filling,
                phase=String(chain.phase),
                case=chain.case,
                seed=chain.seed,
            ),
            diagnostics,
        )

        mktempdir() do directory
            options = SignFeasibilityRunner.parse_runner_options([
                "--output-dir=$directory",
                "--require-single",
                "--strict-pilot",
            ])
            summary_path = joinpath(
                directory,
                SignFeasibilityRunner.SUMMARY_FILENAME,
            )
            metadata_path = joinpath(
                directory,
                SignFeasibilityRunner.METADATA_FILENAME,
            )
            SignFeasibilityRunner.write_summary(summary_path, [summary])
            SignFeasibilityRunner.write_metadata(
                metadata_path,
                [chain],
                DEFAULT_SIGN_FEASIBILITY_CRITERIA,
                [summary],
                options,
            )

            csv_header = first(readlines(summary_path))
            metadata = read(metadata_path, String)
            for field in (
                "full_sign_denominator_usable",
                "leave_one_out_denominators_usable",
                "leave_one_out_denominator_failure_count",
                "minimum_absolute_leave_one_out_average_sign",
                "filling_jackknife_valid",
                "filling_jackknife_standard_error",
                "acceptance_rate",
                "attempted_updates",
                "accepted_updates",
                "rejected_updates",
                "block_drift_zscore",
                "maximum_greens_recompute_error",
                "runtime_seconds",
                "strict_pilot_status",
            )
                @test occursin(field, csv_header)
                @test occursin(field, metadata)
            end
            @test occursin("require_single=true", metadata)
            @test occursin("strict_pilot=true", metadata)
            @test occursin(
                "completed strict-pilot classification is recorded",
                metadata,
            )
            @test !occursin(
                "No average-sign conclusion is encoded",
                metadata,
            )

            nonpilot_summary_path = joinpath(directory, "nonpilot.csv")
            nonpilot_metadata_path = joinpath(directory, "nonpilot_metadata.txt")
            nonpilot_options = SignFeasibilityRunner.parse_runner_options([
                "--output-dir=$directory",
            ])
            SignFeasibilityRunner.write_summary(
                nonpilot_summary_path,
                [base_summary],
            )
            SignFeasibilityRunner.write_metadata(
                nonpilot_metadata_path,
                [chain],
                DEFAULT_SIGN_FEASIBILITY_CRITERIA,
                [base_summary],
                nonpilot_options,
            )
            nonpilot_csv_header = first(readlines(nonpilot_summary_path))
            nonpilot_metadata = read(nonpilot_metadata_path, String)
            @test !occursin("strict_pilot_status", nonpilot_csv_header)
            @test !occursin("strict_pilot_status", nonpilot_metadata)
            @test !occursin("strict_pilot=false", nonpilot_metadata)
            @test occursin(
                "No average-sign conclusion is encoded",
                nonpilot_metadata,
            )

            collision_options =
                SignFeasibilityRunner.parse_runner_options([
                    "--output-dir=$directory",
                ])
            @test_throws ErrorException SignFeasibilityRunner.check_output_collisions(
                collision_options,
            )
            overwrite_options =
                SignFeasibilityRunner.parse_runner_options([
                    "--output-dir=$directory",
                    "--overwrite",
                ])
            @test isnothing(SignFeasibilityRunner.check_output_collisions(
                overwrite_options,
            ))
            no_write_options =
                SignFeasibilityRunner.parse_runner_options([
                    "--output-dir=$directory",
                    "--no-write",
                ])
            @test isnothing(SignFeasibilityRunner.check_output_collisions(
                no_write_options,
            ))
            @test_throws ErrorException SignFeasibilityRunner.parse_runner_options([
                "--overwrite",
                "--no-write",
            ])

            failing_options = SignFeasibilityRunner.parse_runner_options([
                "--output-dir=$(joinpath(directory, "failing"))",
                "--require-single",
                "--strict-pilot",
            ])
            failure_pipe = Pipe()
            redirect_stdout(failure_pipe) do
                @test_throws ErrorException SignFeasibilityRunner.finish_completed_run(
                    failing_options,
                    [chain],
                    DEFAULT_SIGN_FEASIBILITY_CRITERIA,
                    [base_summary];
                    summary_writer=(path, summaries) ->
                        error("synthetic summary write failure"),
                )
            end
            close(failure_pipe.in)
            failure_output = read(failure_pipe, String)
            @test occursin("outputs_complete=false", failure_output)
            @test occursin("status=INCONCLUSIVE", failure_output)
            @test occursin("failed_gates=output", failure_output)
            @test !occursin("status=PASS", failure_output)

            no_write_strict_options =
                SignFeasibilityRunner.parse_runner_options([
                    "--no-write",
                    "--require-single",
                    "--strict-pilot",
                ])
            no_write_pipe = Pipe()
            redirect_stdout(no_write_pipe) do
                @test_throws ErrorException SignFeasibilityRunner.finish_completed_run(
                    no_write_strict_options,
                    [chain],
                    DEFAULT_SIGN_FEASIBILITY_CRITERIA,
                    [base_summary],
                )
            end
            close(no_write_pipe.in)
            no_write_output = read(no_write_pipe, String)
            @test occursin("outputs_written=false", no_write_output)
            @test occursin("status=INCONCLUSIVE", no_write_output)
            @test occursin("failed_gates=output", no_write_output)
        end
    end

    @testset "Raw density and exactly-once sign weighting" begin
        parameters = DQMCParameters(
            L=2,
            beta=0.2,
            dtau=0.1,
            U=0.0,
            mu=-0.3,
            stabilization_period=1,
            burnin_sweeps=0,
            measurement_sweeps=1,
            sweeps_between_measurements=1,
            seed=6601,
        )
        state = initialize_dqmc_state(parameters)
        positive_measurement = measure_basic_observables(state)
        @test positive_measurement.configuration_sign == 1.0

        # The determinant-sign fields affect only configuration_sign. The
        # density is the raw 1-G estimator and must not change when that sign
        # is flipped independently for this convention-focused test.
        state.signdet_green_up = -state.signdet_green_up
        negative_measurement = measure_basic_observables(state)
        @test negative_measurement.configuration_sign == -1.0
        @test negative_measurement.total_density ==
              positive_measurement.total_density

        one_measurement = summarize_sign_feasibility(
            [negative_measurement.configuration_sign],
            [negative_measurement.total_density];
            block_size=1,
            criteria=criteria,
        )
        @test one_measurement.signed_filling_sum ==
              negative_measurement.configuration_sign *
              negative_measurement.total_density

        signs = [1.0, -1.0, 1.0]
        raw_total_densities = [1.25, 2.0, 3.25]
        synthetic = summarize_sign_feasibility(
            signs,
            raw_total_densities;
            block_size=1,
            criteria=SignFeasibilityCriteria(
                average_sign_threshold=-0.9,
                maximum_average_sign_standard_error=1.0,
                minimum_block_count=3,
                maximum_block_drift_zscore=10.0,
                minimum_sign_signal_to_noise=0.0,
                minimum_absolute_average_sign=1.0e-12,
            ),
        )
        @test synthetic.configuration_sign_sum == sum(signs)
        @test synthetic.raw_filling_sum == sum(raw_total_densities)
        @test synthetic.signed_filling_sum == 1.25 - 2.0 + 3.25
        @test synthetic.total_density ==
              synthetic.signed_filling_sum / synthetic.configuration_sign_sum
    end


    @testset "General sign-reweighted thermodynamics" begin
        thermodynamic_names = (
            :density,
            :doping,
            :kinetic_energy_per_cell,
            :shifted_interaction_energy_per_cell,
            :internal_energy_per_cell,
            :grand_energy_per_cell,
        )
        mu = -1.897
        densities = [1.2, 1.7, 2.1, 2.6, 1.4, 2.8, 2.3, 1.9]
        kinetic_energies = [-2.4, -2.1, -1.8, -1.5, -2.2, -1.3, -1.6, -2.0]
        interaction_energies = [0.4, 0.2, -0.1, 0.3, 0.5, -0.2, 0.1, 0.25]
        raw_observations = [
            _synthetic_thermodynamic_observation(
                densities[index],
                kinetic_energies[index],
                interaction_energies[index],
                mu,
            )
            for index in eachindex(densities)
        ]

        @testset "Mixed signs, total ratios, and independent jackknife" begin
            signs = [1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0]
            summary = summarize_sign_reweighted_thermodynamics(
                signs,
                raw_observations;
                block_size=2,
            )
            @test summary.measurement_count == 8
            @test summary.block_count == 4
            @test summary.minimum_block_length == 2
            @test summary.maximum_block_length == 2
            @test summary.configuration_sign_sum == sum(signs) == 4.0
            @test summary.average_sign == 0.5
            @test summary.sign_denominator_usable
            @test summary.leave_one_out_denominators_usable
            @test summary.leave_one_out_denominator_failure_count == 0
            @test summary.minimum_absolute_leave_one_out_average_sign == 1 / 3
            @test summary.denominator_guard_passed
            @test summary.configuration_signs_valid
            @test summary.raw_observables_finite
            @test summary.finite_results

            for name in thermodynamic_names
                values = getproperty.(raw_observations, name)
                independent = _independent_ratio_jackknife(signs, values, 2)
                @test isapprox(
                    getproperty(summary, name),
                    independent.estimate;
                    atol=2.0e-14,
                    rtol=0.0,
                )
                error_name = Symbol(name, "_standard_error")
                @test isapprox(
                    getproperty(summary, error_name),
                    independent.standard_error;
                    atol=2.0e-14,
                    rtol=0.0,
                )
                @test getproperty(summary, Symbol(name, "_jackknife_valid"))
            end

            @test isapprox(
                summary.doping,
                (summary.density - 2) / 4;
                atol=2.0e-14,
                rtol=0.0,
            )
            @test isapprox(
                summary.internal_energy_per_cell,
                summary.kinetic_energy_per_cell +
                summary.shifted_interaction_energy_per_cell;
                atol=2.0e-14,
                rtol=0.0,
            )
            @test isapprox(
                summary.grand_energy_per_cell,
                summary.internal_energy_per_cell - mu * summary.density;
                atol=2.0e-14,
                rtol=0.0,
            )

            compatibility_criteria = SignFeasibilityCriteria(
                average_sign_threshold=-0.9,
                maximum_average_sign_standard_error=1.0,
                minimum_block_count=4,
                maximum_block_drift_zscore=10.0,
                minimum_sign_signal_to_noise=0.0,
                minimum_absolute_average_sign=0.0,
            )
            legacy = summarize_sign_feasibility(
                signs,
                densities;
                block_size=2,
                criteria=compatibility_criteria,
            )
            @test legacy.total_density == summary.density
            @test legacy.total_density_standard_error ==
                  summary.density_standard_error
            @test legacy.filling_jackknife_valid == summary.density_jackknife_valid
            @test legacy.configuration_sign_sum == summary.configuration_sign_sum
            @test legacy.raw_filling_sum == sum(densities)
            @test legacy.signed_filling_sum ==
                  sum(signs[index] * densities[index] for index in eachindex(signs))
            @test propertynames(legacy) == (
                :measurement_count,
                :requested_block_size,
                :block_count,
                :minimum_block_length,
                :maximum_block_length,
                :configuration_sign_sum,
                :raw_filling_sum,
                :signed_filling_sum,
                :average_sign,
                :average_sign_standard_error,
                :average_sign_lower_2se,
                :configured_average_sign_threshold,
                :average_sign_gate_passed,
                :block_count_sufficient,
                :block_stability_evaluable,
                :block_drift_zscore,
                :block_stable,
                :average_sign_precision_sufficient,
                :sampling_sufficient,
                :sign_denominator_usable,
                :full_sign_denominator_usable,
                :leave_one_out_denominators_usable,
                :leave_one_out_denominator_failure_count,
                :minimum_absolute_leave_one_out_average_sign,
                :denominator_guard_passed,
                :sign_signal_to_noise,
                :sign_signal_to_noise_sufficient,
                :filling_resolution_sufficient,
                :total_density,
                :total_density_standard_error,
                :filling_jackknife_standard_error,
                :filling_jackknife_valid,
                :configuration_signs_finite,
                :configuration_signs_valid,
                :raw_fillings_finite,
                :acceptance_rate,
                :acceptance_valid,
                :maximum_greens_recompute_error,
                :correction_guard_passed,
                :finite_results,
                :numerical_valid,
                :correlation_measurement_ready,
            )
        end

        @testset "All-positive signs and final short block" begin
            signs = ones(7)
            observations = raw_observations[1:7]
            summary = summarize_sign_reweighted_thermodynamics(
                signs,
                observations;
                block_size=3,
            )
            @test summary.block_count == 3
            @test summary.minimum_block_length == 1
            @test summary.maximum_block_length == 3
            @test summary.average_sign == 1.0
            @test summary.average_sign_standard_error == 0.0
            @test summary.denominator_guard_passed
            for name in thermodynamic_names
                values = getproperty.(observations, name)
                @test isapprox(
                    getproperty(summary, name),
                    sum(values) / length(values);
                    atol=2.0e-14,
                    rtol=0.0,
                )
                @test getproperty(summary, Symbol(name, "_jackknife_valid"))
            end
        end

        @testset "Full and leave-one-out denominator guards" begin
            zero_total = summarize_sign_reweighted_thermodynamics(
                [1.0, -1.0, 1.0, -1.0],
                raw_observations[1:4];
                block_size=2,
            )
            @test zero_total.configuration_sign_sum == 0.0
            @test !zero_total.sign_denominator_usable
            @test !zero_total.denominator_guard_passed
            for name in thermodynamic_names
                @test ismissing(getproperty(zero_total, name))
                @test ismissing(getproperty(
                    zero_total,
                    Symbol(name, "_standard_error"),
                ))
                @test !getproperty(zero_total, Symbol(name, "_jackknife_valid"))
            end

            zero_leave_one_out = summarize_sign_reweighted_thermodynamics(
                [1.0, 1.0, -1.0, -1.0, 1.0, 1.0],
                raw_observations[1:6];
                block_size=2,
            )
            @test zero_leave_one_out.configuration_sign_sum == 2.0
            @test zero_leave_one_out.sign_denominator_usable
            @test !zero_leave_one_out.leave_one_out_denominators_usable
            @test zero_leave_one_out.leave_one_out_denominator_failure_count == 2
            @test zero_leave_one_out.minimum_absolute_leave_one_out_average_sign ==
                  0.0
            @test !zero_leave_one_out.denominator_guard_passed
            for name in thermodynamic_names
                @test isfinite(getproperty(zero_leave_one_out, name))
                @test ismissing(getproperty(
                    zero_leave_one_out,
                    Symbol(name, "_standard_error"),
                ))
                @test !getproperty(
                    zero_leave_one_out,
                    Symbol(name, "_jackknife_valid"),
                )
            end
        end

        @testset "Input validation" begin
            @test_throws ArgumentError summarize_sign_reweighted_thermodynamics(
                Float64[],
                NamedTuple[];
                block_size=1,
            )
            @test_throws DimensionMismatch summarize_sign_reweighted_thermodynamics(
                [1.0],
                raw_observations[1:2];
                block_size=1,
            )
            @test_throws ArgumentError summarize_sign_reweighted_thermodynamics(
                [1.0],
                raw_observations[1:1];
                block_size=0,
            )
            @test_throws ArgumentError summarize_sign_reweighted_thermodynamics(
                [1.0],
                [(density=1.0,)];
                block_size=1,
            )
            @test_throws ArgumentError summarize_sign_reweighted_thermodynamics(
                [1.0],
                raw_observations[1:1];
                block_size=1,
                minimum_absolute_average_sign=1.0,
            )
        end
    end

    @testset "Independent diagnostic domains and combined readiness" begin
        ready = summarize_sign_feasibility(
            ones(8),
            fill(1.5, 8);
            block_size=2,
            criteria=criteria,
            maximum_greens_recompute_error=5.0e-7,
            acceptance_rate=0.5,
        )
        @test ready.average_sign == 1.0
        @test ready.average_sign_standard_error == 0.0
        @test ready.average_sign_lower_2se == 1.0
        @test ready.average_sign_gate_passed
        @test ready.block_count_sufficient
        @test ready.block_stable
        @test ready.average_sign_precision_sufficient
        @test ready.sampling_sufficient
        @test ready.sign_denominator_usable
        @test ready.sign_signal_to_noise == Inf
        @test ready.sign_signal_to_noise_sufficient
        @test ready.filling_resolution_sufficient
        @test ready.correction_guard_passed
        @test ready.finite_results
        @test ready.numerical_valid
        @test ready.correlation_measurement_ready

        bad_correction = summarize_sign_feasibility(
            ones(8),
            fill(1.5, 8);
            block_size=2,
            criteria=criteria,
            maximum_greens_recompute_error=2.0e-6,
            acceptance_rate=0.5,
        )
        @test bad_correction.average_sign_gate_passed
        @test bad_correction.sampling_sufficient
        @test bad_correction.filling_resolution_sufficient
        @test !bad_correction.correction_guard_passed
        @test !bad_correction.numerical_valid
        @test !bad_correction.correlation_measurement_ready

        zero_denominator = summarize_sign_feasibility(
            [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0],
            fill(2.0, 8);
            block_size=2,
            criteria=criteria,
        )
        @test zero_denominator.configuration_sign_sum == 0.0
        @test zero_denominator.signed_filling_sum == 0.0
        @test !zero_denominator.average_sign_gate_passed
        @test zero_denominator.sampling_sufficient
        @test !zero_denominator.sign_denominator_usable
        @test zero_denominator.sign_signal_to_noise == 0.0
        @test !zero_denominator.filling_resolution_sufficient
        @test ismissing(zero_denominator.total_density)
        @test zero_denominator.numerical_valid
        @test !zero_denominator.correlation_measurement_ready
    end

    @testset "Strict single-pilot classification" begin
        pass_summary = _strict_pilot_summary()
        pass_assessment =
            AltermagnetDQMCTrial.classify_strict_sign_pilot(pass_summary)
        @test pass_assessment.strict_pilot_status == :PASS
        @test pass_assessment.strict_pilot_passed
        @test pass_assessment.strict_pilot_prerequisite_gates_passed
        @test pass_assessment.strict_pilot_sign_threshold_passed
        @test pass_assessment.strict_pilot_sign_precision_passed
        @test pass_assessment.strict_pilot_sign_signal_to_noise_passed
        @test pass_assessment.strict_pilot_sign_signal_to_noise_consistent
        @test pass_assessment.strict_pilot_acceptance_consistency_passed
        @test pass_assessment.strict_pilot_runtime_passed
        @test pass_assessment.strict_pilot_output_gate_passed
        @test pass_summary.average_sign == 1.0
        @test pass_summary.average_sign_standard_error == 0.0
        @test pass_summary.sign_signal_to_noise == Inf
        @test pass_summary.full_sign_denominator_usable
        @test pass_summary.leave_one_out_denominators_usable
        @test pass_summary.leave_one_out_denominator_failure_count == 0
        @test pass_summary.minimum_absolute_leave_one_out_average_sign == 1.0
        @test pass_summary.denominator_guard_passed
        @test pass_summary.filling_jackknife_valid
        @test isfinite(pass_summary.filling_jackknife_standard_error)

        definitive_signs = _repeated_sign_blocks(6, 10, 8)
        definitive_summary = _strict_pilot_summary(
            signs=definitive_signs,
            block_size=10,
        )
        definitive_assessment =
            AltermagnetDQMCTrial.classify_strict_sign_pilot(
                definitive_summary,
        )
        @test definitive_summary.average_sign == 0.2
        # Identical 0.2 block means can differ at the final Float64 summation
        # bit; 2e-16 is a roundoff-only tolerance at this scale.
        @test isapprox(
            definitive_summary.average_sign_standard_error,
            0.0;
            atol=2.0e-16,
            rtol=0.0,
        )
        @test definitive_assessment.strict_pilot_prerequisite_gates_passed
        @test !definitive_assessment.strict_pilot_sign_threshold_passed
        @test definitive_assessment.strict_pilot_status ==
              :DEFINITIVE_SIGN_FAIL
        @test definitive_assessment.strict_pilot_failed_gates ==
              "sign_threshold"

        precision_failure = merge(
            pass_summary,
            (
                average_sign_standard_error=0.03,
                average_sign_lower_2se=0.94,
                sign_signal_to_noise=1.0 / 0.03,
            ),
        )
        loo_failure = merge(
            pass_summary,
            (
                leave_one_out_denominators_usable=false,
                leave_one_out_denominator_failure_count=1,
                minimum_absolute_leave_one_out_average_sign=0.0,
                denominator_guard_passed=false,
            ),
        )
        green_failure = merge(
            pass_summary,
            (
                maximum_greens_recompute_error=nextfloat(1.0e-8),
                correction_guard_passed=false,
            ),
        )
        drift_failure = merge(
            pass_summary,
            (
                block_drift_zscore=Inf,
                block_stable=false,
                sampling_sufficient=false,
            ),
        )
        jackknife_failure = merge(
            pass_summary,
            (
                total_density_standard_error=missing,
                filling_jackknife_standard_error=missing,
                filling_jackknife_valid=false,
            ),
        )
        update_failure = merge(
            pass_summary,
            (
                attempted_updates=100,
                accepted_updates=0,
                rejected_updates=100,
            ),
        )
        missing_acceptance = merge(
            pass_summary,
            (acceptance_rate=missing,),
        )
        nonfinite_acceptance = merge(
            pass_summary,
            (acceptance_rate=NaN,),
        )
        inconsistent_acceptance = merge(
            pass_summary,
            (acceptance_rate=0.9,),
        )
        forged_infinite_sign_snr = merge(
            pass_summary,
            (
                average_sign_standard_error=0.01,
                average_sign_lower_2se=0.98,
                sign_signal_to_noise=Inf,
            ),
        )
        raw_filling_failure = merge(
            pass_summary,
            (
                raw_filling_sum=NaN,
                raw_fillings_finite=false,
            ),
        )
        missing_filling = merge(
            pass_summary,
            (total_density=missing,),
        )
        nonfinite_filling = merge(
            pass_summary,
            (total_density=Inf,),
        )
        block_count_failure = merge(
            pass_summary,
            (
                block_count=7,
                block_count_sufficient=false,
                sampling_sufficient=false,
            ),
        )
        runtime_failure = merge(
            pass_summary,
            (runtime_seconds=Inf,),
        )

        inconclusive_cases = (
            (precision_failure, "sign_precision"),
            (loo_failure, "denominators"),
            (green_failure, "greens_correction"),
            (drift_failure, "block_drift"),
            (jackknife_failure, "filling_jackknife"),
            (update_failure, "update_counts"),
            (missing_acceptance, "acceptance"),
            (nonfinite_acceptance, "acceptance"),
            (inconsistent_acceptance, "acceptance_consistency"),
            (
                forged_infinite_sign_snr,
                "sign_signal_to_noise_consistency",
            ),
            (raw_filling_failure, "raw_filling"),
            (missing_filling, "reweighted_filling"),
            (nonfinite_filling, "reweighted_filling"),
            (block_count_failure, "block_count"),
            (runtime_failure, "runtime"),
        )
        for (summary, failed_gate) in inconclusive_cases
            assessment = AltermagnetDQMCTrial.classify_strict_sign_pilot(
                summary,
            )
            @test assessment.strict_pilot_status == :INCONCLUSIVE
            @test !assessment.strict_pilot_passed
            @test occursin(failed_gate, assessment.strict_pilot_failed_gates)
        end

        procedural_failure =
            AltermagnetDQMCTrial.classify_strict_sign_pilot(
                pass_summary;
                procedural_gate_passed=false,
            )
        @test procedural_failure.strict_pilot_status == :INCONCLUSIVE
        @test occursin(
            "procedural",
            procedural_failure.strict_pilot_failed_gates,
        )

        output_failure = AltermagnetDQMCTrial.classify_strict_sign_pilot(
            pass_summary;
            output_gate_passed=false,
        )
        @test output_failure.strict_pilot_status == :INCONCLUSIVE
        @test !output_failure.strict_pilot_output_gate_passed
        @test occursin("output", output_failure.strict_pilot_failed_gates)

        loo_signs = vcat(
            [1.0, 1.0],
            [-1.0, -1.0],
            [1.0, 1.0],
            [-1.0, -1.0],
            [1.0, 1.0],
            [-1.0, -1.0],
            [1.0, 1.0],
            [1.0, -1.0],
        )
        loo_summary = _strict_pilot_summary(
            signs=loo_signs,
            block_size=2,
        )
        @test loo_summary.full_sign_denominator_usable
        @test !loo_summary.leave_one_out_denominators_usable
        @test loo_summary.leave_one_out_denominator_failure_count > 0
        @test loo_summary.minimum_absolute_leave_one_out_average_sign == 0.0
        @test AltermagnetDQMCTrial.classify_strict_sign_pilot(
            loo_summary,
        ).strict_pilot_status == :INCONCLUSIVE
    end

    @testset "Statistical and numerical edge cases" begin
        one_block = summarize_sign_feasibility(
            ones(2),
            fill(1.0, 2);
            block_size=2,
            criteria=criteria,
        )
        @test ismissing(one_block.average_sign_standard_error)
        @test !one_block.average_sign_gate_passed
        @test !one_block.block_stability_evaluable
        @test !one_block.sampling_sufficient
        @test one_block.numerical_valid

        drifting = summarize_sign_feasibility(
            [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0],
            fill(1.0, 8);
            block_size=2,
            criteria=criteria,
        )
        @test drifting.block_stability_evaluable
        @test drifting.block_drift_zscore == Inf
        @test !drifting.block_stable
        @test !drifting.sampling_sufficient

        imprecise = summarize_sign_feasibility(
            [1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0],
            fill(1.0, 8);
            block_size=2,
            criteria=criteria,
        )
        @test imprecise.block_stable
        @test !imprecise.average_sign_precision_sufficient
        @test !imprecise.sampling_sufficient

        nonfinite = summarize_sign_feasibility(
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, NaN, 1.0];
            block_size=1,
            criteria=criteria,
        )
        @test !nonfinite.raw_fillings_finite
        @test !nonfinite.finite_results
        @test !nonfinite.numerical_valid

        invalid_sign = summarize_sign_feasibility(
            [1.0, 0.0, 1.0, 1.0],
            fill(1.0, 4);
            block_size=1,
            criteria=criteria,
        )
        @test !invalid_sign.configuration_signs_valid
        @test !invalid_sign.numerical_valid

        @test_throws ArgumentError summarize_sign_feasibility(
            Float64[],
            Float64[];
            block_size=1,
        )
        @test_throws DimensionMismatch summarize_sign_feasibility(
            [1.0],
            [1.0, 2.0];
            block_size=1,
        )
        @test_throws ArgumentError summarize_sign_feasibility(
            [1.0],
            [1.0];
            block_size=0,
        )
    end
end

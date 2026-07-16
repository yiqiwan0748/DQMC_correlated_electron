using AltermagnetDQMCTrial
using Test

@testset "Milestone 6A DQMC sign feasibility" begin
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
        @test !average_sign_gate_passes(1 // 2, 1 // 10, 3 // 10)
        @test average_sign_gate_passes(1 // 2, 1 // 10, 299 // 1000)
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

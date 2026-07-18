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
        @test length(chains) == 10
        @test length(u0_controls) == 1
        @test only(u0_controls).U == 0.0
        @test only(u0_controls).mu == -1.897
        @test [chain.mu for chain in calibration] ==
              [-2.497, -2.197, -1.897, -1.597, -1.297]
        @test !isempty(fixed_mu)
        @test all(chain -> chain.mu == -1.897, fixed_mu)
        @test all(chain -> chain.mu == -1.897,
                  filter(chain -> chain.phase != :mu_calibration, chains))
        @test all(chain -> !occursin("lower_energy_vhs", chain.case), chains)
        @test all(chain -> chain.mu != -3.8369423191, chains)

        options = SignFeasibilityRunner.parse_runner_options(String[])
        @test !options.run_requested
        pipe = Pipe()
        redirect_stdout(pipe) do
            @test isnothing(SignFeasibilityRunner.main(String[]))
        end
        close(pipe.in)
        output = read(pipe, String)
        @test occursin("mode=dry-run", output)
        @test occursin("dqmc_chains_constructed=0", output)
        @test occursin("outputs_written=false", output)
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
                :sign_signal_to_noise,
                :sign_signal_to_noise_sufficient,
                :filling_resolution_sufficient,
                :total_density,
                :total_density_standard_error,
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

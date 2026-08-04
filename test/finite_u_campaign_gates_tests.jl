using Test

module FiniteUCampaignGates
include(joinpath(@__DIR__, "..", "src", "finite_u_campaign_gates.jl"))
end

function calibration_point(
    mu,
    density;
    standard_error=0.001,
    jackknife_valid=true,
    completed=true,
    burnin_sweeps=200,
    measurement_sweeps=800,
    sweeps_between_measurements=2,
    block_size=100,
    stabilization_period=5,
    measurement_count=800,
    full_denominator_usable=true,
    leave_one_out_denominators_usable=true,
)
    return (
        mu=mu,
        total_density=density,
        filling_jackknife_standard_error=standard_error,
        filling_jackknife_valid=jackknife_valid,
        completed=completed,
        burnin_sweeps=burnin_sweeps,
        measurement_sweeps=measurement_sweeps,
        sweeps_between_measurements=sweeps_between_measurements,
        block_size=block_size,
        stabilization_period=stabilization_period,
        measurement_count=measurement_count,
        full_sign_denominator_usable=full_denominator_usable,
        leave_one_out_denominators_usable=
            leave_one_out_denominators_usable,
    )
end

function chain_result(
    seed;
    completed=true,
    density=FiniteUCampaignGates.FIXED_DENSITY_TARGET_FILLING,
    density_standard_error=0.001,
    jackknife_valid=true,
    full_denominator_usable=true,
    leave_one_out_denominators_usable=true,
    average_sign=0.5,
    average_sign_standard_error=0.01,
    strict_pilot_status=:PASS,
    burnin_sweeps=400,
    measurement_sweeps=1600,
    sweeps_between_measurements=2,
    block_size=100,
    stabilization_period=5,
    measurement_count=1600,
)
    return (
        seed=seed,
        completed=completed,
        total_density=density,
        filling_jackknife_standard_error=density_standard_error,
        filling_jackknife_valid=jackknife_valid,
        full_sign_denominator_usable=full_denominator_usable,
        leave_one_out_denominators_usable=
            leave_one_out_denominators_usable,
        average_sign=average_sign,
        average_sign_standard_error=average_sign_standard_error,
        strict_pilot_status=strict_pilot_status,
        burnin_sweeps=burnin_sweeps,
        measurement_sweeps=measurement_sweeps,
        sweeps_between_measurements=sweeps_between_measurements,
        block_size=block_size,
        stabilization_period=stabilization_period,
        measurement_count=measurement_count,
    )
end

function resource_record(;
    projected_two_chain_seconds=86400.0,
    peak_rss_bytes=FiniteUCampaignGates.RESOURCE_PEAK_RSS_BYTES_LIMIT,
    swap_bytes=0,
    oom_occurred=false,
    numerical_failure=false,
    maximum_greens_recompute_error=1.0e-8,
)
    return (
        projected_two_chain_seconds=projected_two_chain_seconds,
        peak_rss_bytes=peak_rss_bytes,
        swap_bytes=swap_bytes,
        oom_occurred=oom_occurred,
        numerical_failure=numerical_failure,
        maximum_greens_recompute_error=
            maximum_greens_recompute_error,
    )
end

function without_field(record::NamedTuple, field::Symbol)
    names = Tuple(name for name in propertynames(record) if name != field)
    values = Tuple(getproperty(record, name) for name in names)
    return NamedTuple{names}(values)
end

@testset "Finite-U campaign gates" begin
    @testset "Three-point fixed-density calibration" begin
        target = FiniteUCampaignGates.FIXED_DENSITY_TARGET_FILLING
        @test target == 1.46021875
        @test FiniteUCampaignGates.FIXED_DENSITY_CALIBRATION_WORKLOAD == (
            burnin_sweeps=200,
            measurement_sweeps=800,
            sweeps_between_measurements=2,
            block_size=100,
            stabilization_period=5,
        )

        compatible_points = [
            calibration_point(-2.4602, target + 0.03),
            calibration_point(-2.5802, target - 0.03),
            calibration_point(-2.5202, target),
        ]
        compatible =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration(
                compatible_points,
            )
        @test compatible.status == :TARGET_COMPATIBLE
        @test compatible.advance
        @test compatible.selected_mu == -2.5202
        @test compatible.compatible_point_mu == -2.5202
        @test isnothing(compatible.interpolated_mu_unrounded)
        @test isnothing(compatible.bracket_mus)
        @test compatible.sorted_mus == (-2.5802, -2.5202, -2.4602)
        @test isempty(compatible.reasons)
        @test !compatible.automatic_extension_authorized

        boundary_points = [
            calibration_point(-0.1, 7 // 10; standard_error=1 // 100),
            calibration_point(0.0, 6 // 5; standard_error=1 // 10),
            calibration_point(0.1, 7 // 5; standard_error=1 // 100),
        ]
        boundary =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration(
                boundary_points;
                target_filling=1 // 1,
            )
        @test boundary.status == :TARGET_COMPATIBLE
        @test boundary.selected_mu == 0.0

        tied_compatible_points = [
            calibration_point(-0.1, 9 // 10; standard_error=1 // 10),
            calibration_point(0.0, 11 // 10; standard_error=1 // 10),
            calibration_point(0.1, 13 // 10; standard_error=1 // 100),
        ]
        tied =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration(
                tied_compatible_points;
                target_filling=1 // 1,
            )
        @test tied.status == :TARGET_COMPATIBLE
        @test tied.selected_mu == -0.1

        lower_density = 1.4501234567
        upper_density = 1.4709876543
        interpolating_points = [
            calibration_point(-2.4602, 1.5001234567;
                              standard_error=1.0e-7),
            calibration_point(-2.5802, lower_density;
                              standard_error=1.0e-7),
            calibration_point(-2.5202, upper_density;
                              standard_error=1.0e-7),
        ]
        interpolated =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration(
                interpolating_points,
            )
        expected_unrounded =
            -2.5802 +
            (target - lower_density) *
            (-2.5202 - -2.5802) /
            (upper_density - lower_density)
        @test interpolated.status == :INTERPOLATED
        @test interpolated.advance
        @test interpolated.interpolated_mu_unrounded == expected_unrounded
        @test interpolated.selected_mu ==
              round(expected_unrounded; digits=4)
        @test interpolated.bracket_mus == (-2.5802, -2.5202)
        @test isnothing(interpolated.compatible_point_mu)
        @test isempty(interpolated.reasons)
        @test !interpolated.automatic_extension_authorized

        wrong_count =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration(
                compatible_points[1:2],
            )
        @test wrong_count.status == :STOP_INVALID
        @test !wrong_count.advance
        @test wrong_count.reasons == (:requires_exactly_three_points,)
        @test isnothing(wrong_count.selected_mu)
        @test !wrong_count.automatic_extension_authorized

        malformed = merge(
            calibration_point(-2.5802, target - 0.02),
            (
                mu=NaN,
                total_density=Inf,
                filling_jackknife_standard_error=0.0,
                filling_jackknife_valid=false,
                completed=false,
                burnin_sweeps=201,
                measurement_sweeps=799,
                sweeps_between_measurements=1,
                block_size=99,
                stabilization_period=4,
                measurement_count=799,
                full_sign_denominator_usable=false,
                leave_one_out_denominators_usable=false,
            ),
        )
        malformed_decision =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration([
                malformed,
                compatible_points[2],
                compatible_points[3],
            ])
        @test malformed_decision.status == :STOP_INVALID
        @test !malformed_decision.advance
        for reason in (
            :point_1_mu_not_finite,
            :point_1_density_not_finite,
            :point_1_density_standard_error_not_positive_finite,
            :point_1_density_jackknife_not_valid,
            :point_1_not_completed,
            :point_1_full_denominator_not_usable,
            :point_1_leave_one_out_denominators_not_usable,
            :point_1_burnin_sweeps_mismatch,
            :point_1_measurement_sweeps_mismatch,
            :point_1_sweeps_between_measurements_mismatch,
            :point_1_block_size_mismatch,
            :point_1_stabilization_period_mismatch,
            :point_1_measurement_count_mismatch,
        )
            @test reason in malformed_decision.reasons
        end
        @test !malformed_decision.automatic_extension_authorized

        missing_density = without_field(
            compatible_points[1],
            :total_density,
        )
        missing =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration([
                missing_density,
                compatible_points[2],
                compatible_points[3],
            ])
        @test missing.status == :STOP_INVALID
        @test missing.reasons == (:point_1_missing_total_density,)

        duplicate_mu =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration([
                calibration_point(-2.5, target - 0.03),
                calibration_point(-2.5, target),
                calibration_point(-2.4, target + 0.03),
            ])
        @test duplicate_mu.status == :STOP_INVALID
        @test duplicate_mu.reasons ==
              (:mu_values_not_strictly_increasing,)

        nonmonotone =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration([
                calibration_point(-2.6, target - 0.03),
                calibration_point(-2.5, target + 0.03),
                calibration_point(-2.4, target + 0.02),
            ])
        @test nonmonotone.status == :STOP_INVALID
        @test nonmonotone.reasons ==
              (:density_not_monotone_nondecreasing,)

        unbracketed =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration([
                calibration_point(-2.6, target + 0.02;
                                  standard_error=0.001),
                calibration_point(-2.5, target + 0.03;
                                  standard_error=0.001),
                calibration_point(-2.4, target + 0.04;
                                  standard_error=0.001),
            ])
        @test unbracketed.status == :STOP_UNBRACKETED
        @test !unbracketed.advance
        @test unbracketed.reasons ==
              (:target_not_compatible_or_adjacent_mean_bracketed,)
        @test isnothing(unbracketed.selected_mu)
        @test !unbracketed.automatic_extension_authorized

        invalid_target =
            FiniteUCampaignGates.assess_three_point_fixed_density_calibration(
                compatible_points;
                target_filling=NaN,
            )
        @test invalid_target.status == :STOP_INVALID
        @test invalid_target.reasons == (:target_filling_not_finite,)

        @test_throws ArgumentError FiniteUCampaignGates.assess_three_point_fixed_density_calibration(
            compatible_points;
            expected_workload=(burnin_sweeps=200,),
        )
        @test_throws ArgumentError FiniteUCampaignGates.assess_three_point_fixed_density_calibration(
            compatible_points;
            expected_workload=(
                burnin_sweeps=200,
                measurement_sweeps=0,
                sweeps_between_measurements=2,
                block_size=100,
                stabilization_period=5,
            ),
        )
    end

    @testset "Registered two-chain gate" begin
        target = FiniteUCampaignGates.FIXED_DENSITY_TARGET_FILLING
        @test FiniteUCampaignGates.FIXED_DENSITY_PRODUCTION_WORKLOAD == (
            burnin_sweeps=400,
            measurement_sweeps=1600,
            sweeps_between_measurements=2,
            block_size=100,
            stabilization_period=5,
        )
        first_pass = chain_result(
            6311;
            average_sign=0.52,
            average_sign_standard_error=0.01,
            strict_pilot_status="pass",
        )
        second_pass = chain_result(
            6312;
            average_sign=0.44,
            average_sign_standard_error=0.02,
        )
        advance = FiniteUCampaignGates.assess_two_chain_gate(
            [second_pass, first_pass];
            expected_seeds=(6311, 6312),
        )
        @test advance.status == :ADVANCE
        @test advance.advance
        @test isempty(advance.reasons)
        @test advance.expected_seeds == (6311, 6312)
        @test advance.observed_seeds == (6312, 6311)
        @test advance.density_compatible_seeds == (6312, 6311)
        @test advance.conservative_average_sign_lower_2se == 0.4

        boundary_pass = FiniteUCampaignGates.assess_two_chain_gate(
            [
                chain_result(
                    1;
                    density=6 // 5,
                    density_standard_error=1 // 10,
                ),
                chain_result(
                    2;
                    density=4 // 5,
                    density_standard_error=1 // 10,
                ),
            ];
            expected_seeds=(1, 2),
            target_filling=1 // 1,
        )
        @test boundary_pass.status == :ADVANCE
        @test boundary_pass.density_compatible_seeds == (1, 2)

        definitive = FiniteUCampaignGates.assess_two_chain_gate(
            [
                merge(first_pass, (
                    strict_pilot_status=:DEFINITIVE_SIGN_FAIL,
                )),
                second_pass,
            ];
            expected_seeds=(6311, 6312),
        )
        @test definitive.status == :DEFINITIVE_SIGN_FAIL
        @test !definitive.advance
        @test definitive.reasons ==
              (:strict_pilot_definitive_sign_fail,)
        @test definitive.conservative_average_sign_lower_2se == 0.4

        definitive_with_inconclusive =
            FiniteUCampaignGates.assess_two_chain_gate(
                [
                    merge(first_pass, (
                        strict_pilot_status=:DEFINITIVE_SIGN_FAIL,
                    )),
                    merge(second_pass, (
                        strict_pilot_status=:INCONCLUSIVE,
                    )),
                ];
                expected_seeds=(6311, 6312),
            )
        @test definitive_with_inconclusive.status == :NOT_ADJUDICABLE
        @test definitive_with_inconclusive.reasons ==
              (:strict_pilot_inconclusive,)

        inconclusive = FiniteUCampaignGates.assess_two_chain_gate(
            [
                first_pass,
                merge(second_pass, (
                    strict_pilot_status=:INCONCLUSIVE,
                )),
            ];
            expected_seeds=(6311, 6312),
        )
        @test inconclusive.status == :NOT_ADJUDICABLE
        @test !inconclusive.advance
        @test inconclusive.reasons == (:strict_pilot_inconclusive,)

        missing_result = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass];
            expected_seeds=(6311, 6312),
        )
        @test missing_result.status == :NOT_ADJUDICABLE
        @test missing_result.reasons ==
              (:requires_exactly_two_chain_results,)

        bad_expected_count = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass, second_pass];
            expected_seeds=(6311,),
        )
        @test bad_expected_count.reasons ==
              (:requires_exactly_two_expected_seeds,)
        duplicate_expected = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass, second_pass];
            expected_seeds=(6311, 6311),
        )
        @test duplicate_expected.reasons ==
              (:expected_seeds_must_be_distinct,)
        noninteger_expected = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass, second_pass];
            expected_seeds=(6311, 6312.0),
        )
        @test noninteger_expected.reasons ==
              (:expected_seeds_must_be_integers,)

        duplicate_results = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass, first_pass];
            expected_seeds=(6311, 6312),
        )
        @test duplicate_results.status == :NOT_ADJUDICABLE
        @test :result_seeds_must_be_distinct in duplicate_results.reasons
        @test :missing_expected_seed in duplicate_results.reasons

        unexpected = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass, chain_result(6399)];
            expected_seeds=(6311, 6312),
        )
        @test unexpected.status == :NOT_ADJUDICABLE
        @test :unexpected_result_seed in unexpected.reasons
        @test :missing_expected_seed in unexpected.reasons

        missing_field = FiniteUCampaignGates.assess_two_chain_gate(
            [
                without_field(first_pass, :strict_pilot_status),
                second_pass,
            ];
            expected_seeds=(6311, 6312),
        )
        @test missing_field.status == :NOT_ADJUDICABLE
        @test missing_field.reasons ==
              (:chain_1_missing_strict_pilot_status,)

        malformed_chain = merge(
            first_pass,
            (
                completed=false,
                total_density=NaN,
                filling_jackknife_standard_error=0.0,
                filling_jackknife_valid=false,
                full_sign_denominator_usable=false,
                leave_one_out_denominators_usable=false,
                average_sign=Inf,
                average_sign_standard_error=-0.1,
                strict_pilot_status=:UNKNOWN,
            ),
        )
        density_miss = merge(
            second_pass,
            (total_density=target + 0.1,),
        )
        malformed_pair = FiniteUCampaignGates.assess_two_chain_gate(
            [malformed_chain, density_miss];
            expected_seeds=(6311, 6312),
        )
        @test malformed_pair.status == :NOT_ADJUDICABLE
        @test !malformed_pair.advance
        @test isnothing(
            malformed_pair.conservative_average_sign_lower_2se,
        )
        for reason in (
            :chain_1_not_completed,
            :chain_1_density_not_finite,
            :chain_1_density_standard_error_not_positive_finite,
            :chain_1_density_jackknife_not_valid,
            :chain_1_full_denominator_not_usable,
            :chain_1_leave_one_out_denominators_not_usable,
            :chain_1_average_sign_not_finite,
            :chain_1_average_sign_standard_error_not_nonnegative_finite,
            :chain_1_strict_pilot_status_unrecognized,
            :chain_2_density_not_target_compatible,
        )
            @test reason in malformed_pair.reasons
        end

        finite_density_miss = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass, density_miss];
            expected_seeds=(6311, 6312),
        )
        @test finite_density_miss.status == :NOT_ADJUDICABLE
        @test finite_density_miss.conservative_average_sign_lower_2se ==
              0.4
        @test finite_density_miss.reasons ==
              (:chain_2_density_not_target_compatible,)

        workload_miss = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass, merge(second_pass, (
                measurement_sweeps=1599,
                measurement_count=1599,
            ))];
            expected_seeds=(6311, 6312),
        )
        @test workload_miss.status == :NOT_ADJUDICABLE
        @test workload_miss.reasons == (
            :chain_2_measurement_sweeps_mismatch,
            :chain_2_measurement_count_mismatch,
        )

        invalid_target = FiniteUCampaignGates.assess_two_chain_gate(
            [first_pass, second_pass];
            expected_seeds=(6311, 6312),
            target_filling=Inf,
        )
        @test invalid_target.status == :NOT_ADJUDICABLE
        @test invalid_target.reasons == (:target_filling_not_finite,)
    end

    @testset "Resource gate" begin
        @test FiniteUCampaignGates.RESOURCE_PROJECTED_TWO_CHAIN_SECONDS_LIMIT ==
              86400.0
        @test FiniteUCampaignGates.RESOURCE_PEAK_RSS_BYTES_LIMIT ==
              5_905_580_032
        @test FiniteUCampaignGates.RESOURCE_GREENS_RECOMPUTE_ERROR_LIMIT ==
              1.0e-8

        boundary = FiniteUCampaignGates.assess_resource_gate(
            resource_record(),
        )
        @test boundary.status == :RESOURCE_FEASIBLE
        @test boundary.feasible
        @test isempty(boundary.reasons)
        @test boundary.limits.projected_two_chain_seconds == 86400.0
        @test boundary.limits.peak_rss_bytes == 5_905_580_032
        @test boundary.limits.maximum_greens_recompute_error == 1.0e-8

        zero_metrics = FiniteUCampaignGates.assess_resource_gate(
            resource_record(
                projected_two_chain_seconds=0.0,
                peak_rss_bytes=0,
                maximum_greens_recompute_error=0.0,
            ),
        )
        @test zero_metrics.status == :RESOURCE_FEASIBLE

        failed = FiniteUCampaignGates.assess_resource_gate(
            resource_record(
                projected_two_chain_seconds=nextfloat(86400.0),
                peak_rss_bytes=5_905_580_033,
                swap_bytes=1,
                oom_occurred=true,
                numerical_failure=true,
                maximum_greens_recompute_error=nextfloat(1.0e-8),
            ),
        )
        @test failed.status == :RESOURCE_INFEASIBLE
        @test !failed.feasible
        @test failed.reasons == (
            :projected_two_chain_time_exceeds_limit,
            :peak_rss_exceeds_limit,
            :swap_used,
            :oom_occurred,
            :numerical_failure,
            :greens_recompute_error_exceeds_limit,
        )

        missing = FiniteUCampaignGates.assess_resource_gate(
            without_field(resource_record(), :peak_rss_bytes),
        )
        @test missing.status == :NOT_ADJUDICABLE
        @test !missing.feasible
        @test missing.reasons == (:missing_peak_rss_bytes,)
        @test isnothing(missing.observed)

        malformed = FiniteUCampaignGates.assess_resource_gate(
            resource_record(
                projected_two_chain_seconds=-1.0,
                peak_rss_bytes=Inf,
                swap_bytes=-1,
                oom_occurred=0,
                numerical_failure="false",
                maximum_greens_recompute_error=NaN,
            ),
        )
        @test malformed.status == :NOT_ADJUDICABLE
        @test !malformed.feasible
        @test malformed.reasons == (
            :projected_two_chain_seconds_not_nonnegative_finite,
            :peak_rss_bytes_not_nonnegative_finite,
            :swap_bytes_not_nonnegative_finite,
            :oom_occurred_not_boolean,
            :numerical_failure_not_boolean,
            :maximum_greens_recompute_error_not_nonnegative_finite,
        )
    end
end

using LinearAlgebra
using Test

const ALM_OBSERVABLE_ATOL = 2.0e-12

function _alm_test_measure(greens_up, greens_down, L; configuration_sign=1.0)
    return AltermagnetDQMCTrial._altermagnetic_observables_from_greens(
        greens_up,
        greens_down,
        L,
        configuration_sign,
    )
end

function _alm_test_green_from_occupations(occupations)
    return Matrix(Diagonal(one(eltype(occupations)) .- occupations))
end

function _alm_test_orbital_swap(matrix, L)
    dimension = 2 * L^2
    permutation = collect(1:dimension)
    for site in 0:(L^2 - 1)
        index_x = 2 * site + 1
        index_y = index_x + 1
        permutation[index_x] = index_y
        permutation[index_y] = index_x
    end
    return matrix[permutation, permutation]
end

function _alm_test_density_pair(
    density_up,
    density_down,
    spin_a,
    index_a,
    spin_b,
    index_b,
)
    density_a = spin_a === :up ? density_up : density_down
    density_b = spin_b === :up ? density_up : density_down
    if spin_a === spin_b
        return density_a[index_a, index_a] * density_a[index_b, index_b] +
               (index_a == index_b ? density_a[index_a, index_a] : 0.0) -
               density_a[index_a, index_b] * density_a[index_b, index_a]
    end
    return density_a[index_a, index_a] * density_b[index_b, index_b]
end

"""Independent direct Wick sum over all four spin-orbital flavors."""
function _alm_test_wick_oracle(density_up, density_down, L)
    number_of_sites = L^2
    dimension = 2 * number_of_sites
    flavors = Tuple{Symbol,Int,Float64}[]
    for index in 1:dimension
        orbital_weight = isodd(index) ? 0.5 : -0.5
        push!(flavors, (:up, index, orbital_weight))
        push!(flavors, (:down, index, -orbital_weight))
    end

    second_moment = 0.0
    for (spin_a, index_a, weight_a) in flavors
        for (spin_b, index_b, weight_b) in flavors
            second_moment += weight_a * weight_b * _alm_test_density_pair(
                density_up,
                density_down,
                spin_a,
                index_a,
                spin_b,
                index_b,
            )
        end
    end

    mx = sum(
        density_up[index, index] - density_down[index, index]
        for index in 1:2:dimension
    ) / number_of_sites
    my = sum(
        density_up[index, index] - density_down[index, index]
        for index in 2:2:dimension
    ) / number_of_sites
    return (
        mx=mx,
        my=my,
        m_fm=(mx + my) / 2,
        delta_alm=(mx - my) / 2,
        S_alm_q0=second_moment / number_of_sites,
        m_alm_squared=second_moment / number_of_sites^2,
    )
end

@testset "Altermagnetic equal-time observables" begin
    @testset "Eq. 4 product states" begin
        greens_up = _alm_test_green_from_occupations([1.0, 0.0])
        greens_down = _alm_test_green_from_occupations([0.0, 1.0])
        ordered = _alm_test_measure(greens_up, greens_down, 1)
        @test propertynames(ordered) == (
            :mx,
            :my,
            :m_fm,
            :delta_alm,
            :S_alm_q0,
            :m_alm_squared,
            :configuration_sign,
        )
        @test ordered.mx == 1.0
        @test ordered.my == -1.0
        @test ordered.m_fm == 0.0
        @test ordered.delta_alm == 1.0
        @test ordered.S_alm_q0 == 1.0
        @test ordered.m_alm_squared == 1.0
        @test ordered.configuration_sign == 1.0

        ferromagnetic_up = _alm_test_green_from_occupations([1.0, 1.0])
        ferromagnetic_down = _alm_test_green_from_occupations([0.0, 0.0])
        ferromagnetic = _alm_test_measure(
            ferromagnetic_up,
            ferromagnetic_down,
            1,
        )
        @test ferromagnetic.mx == 1.0
        @test ferromagnetic.my == 1.0
        @test ferromagnetic.m_fm == 1.0
        @test ferromagnetic.delta_alm == 0.0
        @test ferromagnetic.S_alm_q0 == 0.0
        @test ferromagnetic.m_alm_squared == 0.0
    end

    @testset "Same-mode contact terms" begin
        # Four independent half-filled modes have <Delta_alm> = 0 but
        # <Delta_alm^2> = 1/4.  Squaring the one-point estimator would
        # incorrectly return zero and therefore fail this test.
        half_filled_green = Matrix(0.5I, 2, 2)
        measured = _alm_test_measure(
            half_filled_green,
            half_filled_green,
            1,
        )
        @test measured.delta_alm == 0.0
        @test measured.S_alm_q0 == 0.25
        @test measured.m_alm_squared == 0.25
    end

    @testset "Independent Wick oracle" begin
        density_up = [0.67 0.19; 0.19 0.28]
        density_down = [0.21 -0.11; -0.11 0.72]
        identity_matrix = Matrix{Float64}(I, 2, 2)
        greens_up = identity_matrix - density_up
        greens_down = identity_matrix - density_down
        measured = _alm_test_measure(
            greens_up,
            greens_down,
            1;
            configuration_sign=-1.0,
        )
        expected = _alm_test_wick_oracle(density_up, density_down, 1)
        for name in propertynames(expected)
            @test isapprox(
                getproperty(measured, name),
                getproperty(expected, name);
                atol=ALM_OBSERVABLE_ATOL,
                rtol=0.0,
            )
        end
        @test measured.configuration_sign == -1.0
    end

    @testset "Spin and orbital symmetries" begin
        density_up = [0.67 0.19; 0.19 0.28]
        density_down = [0.21 -0.11; -0.11 0.72]
        identity_matrix = Matrix{Float64}(I, 2, 2)
        greens_up = identity_matrix - density_up
        greens_down = identity_matrix - density_down
        original = _alm_test_measure(greens_up, greens_down, 1)
        spin_flipped = _alm_test_measure(greens_down, greens_up, 1)
        orbital_flipped = _alm_test_measure(
            _alm_test_orbital_swap(greens_up, 1),
            _alm_test_orbital_swap(greens_down, 1),
            1,
        )

        @test isapprox(spin_flipped.mx, -original.mx; atol=ALM_OBSERVABLE_ATOL)
        @test isapprox(spin_flipped.my, -original.my; atol=ALM_OBSERVABLE_ATOL)
        @test isapprox(
            spin_flipped.m_fm,
            -original.m_fm;
            atol=ALM_OBSERVABLE_ATOL,
        )
        @test isapprox(
            spin_flipped.delta_alm,
            -original.delta_alm;
            atol=ALM_OBSERVABLE_ATOL,
        )
        @test isapprox(orbital_flipped.mx, original.my; atol=ALM_OBSERVABLE_ATOL)
        @test isapprox(orbital_flipped.my, original.mx; atol=ALM_OBSERVABLE_ATOL)
        @test isapprox(
            orbital_flipped.m_fm,
            original.m_fm;
            atol=ALM_OBSERVABLE_ATOL,
        )
        @test isapprox(
            orbital_flipped.delta_alm,
            -original.delta_alm;
            atol=ALM_OBSERVABLE_ATOL,
        )
        for transformed in (spin_flipped, orbital_flipped)
            @test isapprox(
                transformed.delta_alm^2,
                original.delta_alm^2;
                atol=ALM_OBSERVABLE_ATOL,
            )
            @test isapprox(
                transformed.S_alm_q0,
                original.S_alm_q0;
                atol=ALM_OBSERVABLE_ATOL,
            )
            @test isapprox(
                transformed.m_alm_squared,
                original.m_alm_squared;
                atol=ALM_OBSERVABLE_ATOL,
            )
        end
    end

    @testset "Lattice normalization and live state" begin
        L = 2
        dimension = 2 * L^2
        occupations_up = [isodd(index) ? 1.0 : 0.0 for index in 1:dimension]
        occupations_down = [isodd(index) ? 0.0 : 1.0 for index in 1:dimension]
        ordered = _alm_test_measure(
            _alm_test_green_from_occupations(occupations_up),
            _alm_test_green_from_occupations(occupations_down),
            L,
        )
        @test ordered.delta_alm == 1.0
        @test ordered.S_alm_q0 == L^2
        @test ordered.m_alm_squared == 1.0

        parameters = DQMCParameters(
            L=2,
            beta=0.2,
            dtau=0.1,
            U=0.0,
            stabilization_period=1,
            burnin_sweeps=0,
            measurement_sweeps=1,
            seed=8461,
        )
        state = initialize_dqmc_state(parameters)
        live = AltermagnetDQMCTrial.measure_altermagnetic_observables(state)
        @test live.configuration_sign ==
              measure_basic_observables(state).configuration_sign
        @test live.m_alm_squared == live.S_alm_q0 / L^2
        @test all(isfinite, values(live))
    end

    @testset "Dimension validation" begin
        @test_throws ArgumentError _alm_test_measure(zeros(2, 2), zeros(2, 2), 0)
        @test_throws DimensionMismatch _alm_test_measure(
            zeros(2, 3),
            zeros(2, 2),
            1,
        )
        @test_throws DimensionMismatch _alm_test_measure(
            zeros(2, 2),
            zeros(4, 4),
            1,
        )
    end
end

function _alm_test_synthetic_observation(value, configuration_sign)
    mx = value
    my = 0.25 - 0.5 * value
    S_alm_q0 = value^2 + 0.125
    return (
        mx=mx,
        my=my,
        m_fm=(mx + my) / 2,
        delta_alm=(mx - my) / 2,
        S_alm_q0=S_alm_q0,
        m_alm_squared=S_alm_q0 / 4,
        configuration_sign=configuration_sign,
    )
end

function _alm_test_independent_ratio_jackknife(signs, values, block_size)
    measurement_count = length(signs)
    block_ranges = [
        first_index:min(first_index + block_size - 1, measurement_count)
        for first_index in 1:block_size:measurement_count
    ]
    sign_sum = sum(signs)
    signed_sum = sum(
        signs[index] * values[index] for index in eachindex(signs)
    )
    leave_one_out = [
        begin
            block_sign_sum = sum(signs[index] for index in block_range)
            block_signed_sum = sum(
                signs[index] * values[index] for index in block_range
            )
            (signed_sum - block_signed_sum) / (sign_sum - block_sign_sum)
        end
        for block_range in block_ranges
    ]
    center = sum(leave_one_out) / length(leave_one_out)
    standard_error = sqrt(
        (length(leave_one_out) - 1) / length(leave_one_out) *
        sum((estimate - center)^2 for estimate in leave_one_out),
    )
    return (estimate=signed_sum / sign_sum, standard_error=standard_error)
end

@testset "Sign-reweighted altermagnetic summaries" begin
    observable_names = (
        :mx,
        :my,
        :m_fm,
        :delta_alm,
        :S_alm_q0,
        :m_alm_squared,
    )
    signs = [1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0]
    raw_observations = [
        _alm_test_synthetic_observation(0.2 * index - 0.7, signs[index])
        for index in eachindex(signs)
    ]

    @testset "Exactly-once weighting and common jackknife" begin
        summary = AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
            signs,
            raw_observations;
            block_size=2,
        )
        @test summary.measurement_count == 8
        @test summary.block_count == 4
        @test summary.minimum_block_length == 2
        @test summary.maximum_block_length == 2
        @test summary.configuration_sign_sum == 4.0
        @test summary.average_sign == 0.5
        @test summary.sign_denominator_usable
        @test summary.leave_one_out_denominators_usable
        @test summary.leave_one_out_denominator_failure_count == 0
        @test summary.minimum_absolute_leave_one_out_average_sign == 1 / 3
        @test summary.denominator_guard_passed
        @test summary.configuration_signs_finite
        @test summary.configuration_signs_valid
        @test summary.raw_observables_finite
        @test summary.finite_results

        for name in observable_names
            values = getproperty.(raw_observations, name)
            independent = _alm_test_independent_ratio_jackknife(signs, values, 2)
            @test isapprox(
                getproperty(summary, name),
                independent.estimate;
                atol=2.0e-14,
                rtol=0.0,
            )
            @test isapprox(
                getproperty(summary, Symbol(name, "_standard_error")),
                independent.standard_error;
                atol=2.0e-14,
                rtol=0.0,
            )
            @test getproperty(summary, Symbol(name, "_jackknife_valid"))
        end

        # The embedded signs equal the separate signs. Multiplying by both
        # would cancel the signs and produce this unequal unweighted mean.
        unweighted_mx = sum(getproperty.(raw_observations, :mx)) /
                        length(raw_observations)
        @test !isapprox(summary.mx, unweighted_mx; atol=2.0e-14, rtol=0.0)
    end

    @testset "Common failed leave-one-out denominator" begin
        guarded_signs = [1.0, 1.0, -1.0, -1.0, 1.0, 1.0]
        guarded_observations = raw_observations[1:6]
        summary = AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
            guarded_signs,
            guarded_observations;
            block_size=2,
        )
        @test summary.configuration_sign_sum == 2.0
        @test summary.sign_denominator_usable
        @test !summary.leave_one_out_denominators_usable
        @test summary.leave_one_out_denominator_failure_count == 2
        @test summary.minimum_absolute_leave_one_out_average_sign == 0.0
        @test !summary.denominator_guard_passed
        for name in observable_names
            @test isfinite(getproperty(summary, name))
            @test ismissing(getproperty(summary, Symbol(name, "_standard_error")))
            @test !getproperty(summary, Symbol(name, "_jackknife_valid"))
        end
    end

    @testset "Final short block" begin
        positive_signs = ones(7)
        observations = raw_observations[1:7]
        summary = AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
            positive_signs,
            observations;
            block_size=3,
        )
        @test summary.block_count == 3
        @test summary.minimum_block_length == 1
        @test summary.maximum_block_length == 3
        @test summary.average_sign == 1.0
        @test summary.average_sign_standard_error == 0.0
        @test summary.denominator_guard_passed
        for name in observable_names
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

    @testset "Finiteness and input validation" begin
        nonfinite_raw = copy(raw_observations[1:4])
        nonfinite_raw[2] = merge(nonfinite_raw[2], (S_alm_q0=NaN,))
        raw_summary =
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                ones(4),
                nonfinite_raw;
                block_size=2,
            )
        @test raw_summary.configuration_signs_finite
        @test raw_summary.configuration_signs_valid
        @test !raw_summary.raw_observables_finite
        @test !raw_summary.finite_results
        @test ismissing(raw_summary.S_alm_q0)
        @test !raw_summary.S_alm_q0_jackknife_valid
        @test isfinite(raw_summary.mx)

        sign_summary =
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                [1.0, NaN, 1.0, -1.0],
                raw_observations[1:4];
                block_size=2,
            )
        @test !sign_summary.configuration_signs_finite
        @test !sign_summary.configuration_signs_valid
        @test sign_summary.raw_observables_finite
        @test !sign_summary.finite_results
        @test !sign_summary.sign_denominator_usable
        for name in observable_names
            @test ismissing(getproperty(sign_summary, name))
        end

        invalid_sign_summary =
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                [1.0, 0.0, 1.0, 1.0],
                raw_observations[1:4];
                block_size=2,
            )
        @test invalid_sign_summary.configuration_signs_finite
        @test !invalid_sign_summary.configuration_signs_valid

        @test_throws ArgumentError (
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                Float64[],
                NamedTuple[];
                block_size=1,
            )
        )
        @test_throws DimensionMismatch (
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                [1.0],
                raw_observations[1:2];
                block_size=1,
            )
        )
        @test_throws ArgumentError (
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                [1.0],
                raw_observations[1:1];
                block_size=0,
            )
        )
        @test_throws ArgumentError (
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                [1.0],
                [(mx=1.0,)];
                block_size=1,
            )
        )
        nonreal = merge(raw_observations[1], (delta_alm="not real",))
        @test_throws ArgumentError (
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                [1.0],
                [nonreal];
                block_size=1,
            )
        )
        @test_throws ArgumentError (
            AltermagnetDQMCTrial.summarize_sign_reweighted_altermagnetism(
                [1.0],
                raw_observations[1:1];
                block_size=1,
                minimum_absolute_average_sign=1.0,
            )
        )
    end
end

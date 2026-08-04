using LinearAlgebra
using Serialization
using Test

const UNEQUAL_TIME_ATOL = 5.0e-11

function _serialized_bytes(value)
    buffer = IOBuffer()
    serialize(buffer, value)
    return take!(buffer)
end

function _exact_free_momentum_green(kx, ky, tau, beta, mu)
    kinetic_mu = Matrix(bloch_hamiltonian(kx, ky))
    for orbital in 1:2
        kinetic_mu[orbital, orbital] -= mu
    end
    decomposition = eigen(Symmetric(kinetic_mu))
    spectral_green = [
        exp(-tau * energy) / (1 + exp(-beta * energy))
        for energy in decomposition.values
    ]
    return decomposition.vectors * Diagonal(spectral_green) *
           transpose(decomposition.vectors)
end

@testset "Unequal-time Green measurement boundary" begin
    L = 3
    beta = 0.4
    dtau = 0.1
    mu = -0.37
    parameters = DQMCParameters(
        L=L,
        beta=beta,
        dtau=dtau,
        U=0.0,
        mu=mu,
        stabilization_period=2,
        burnin_sweeps=0,
        measurement_sweeps=1,
        seed=9101,
    )
    dimensions = validate_dqmc_parameters(parameters)
    fixed_fields = fill(Int8(1), dimensions.Nsp, dimensions.Ltau)
    state = initialize_dqmc_state(parameters; hs_fields=fixed_fields)
    state_before = _serialized_bytes(state)

    snapshot = AltermagnetDQMCTrial.measure_unequal_time_greens(state)

    @test snapshot isa AltermagnetDQMCTrial.UnequalTimeGreensSnapshot{Float64}
    @test snapshot.tau == [time_slice * dtau for time_slice in 0:dimensions.Ltau]
    @test snapshot.momenta == allowed_momenta(L)
    @test size(snapshot.real_space) == (2, 2, 2, L, L, dimensions.Ltau + 1)
    @test size(snapshot.momentum_space) ==
          (2, 2, 2, L, L, dimensions.Ltau + 1)
    @test snapshot.configuration_sign == 1.0
    @test all(isfinite, snapshot.real_space)
    @test all(isfinite, snapshot.momentum_space)
    @test maximum(snapshot.maximum_stabilization_error) < 1.0e-10

    # The complete sampling state, including RNG and mutable framework
    # calculators/workspaces, must be bit-for-bit unchanged by measurement.
    @test _serialized_bytes(state) == state_before

    @testset "Exact finite-lattice U=0 oracle and spin symmetry" begin
        for (time_index, tau) in pairs(snapshot.tau),
            (ky_index, ky) in pairs(snapshot.momenta),
            (kx_index, kx) in pairs(snapshot.momenta)
            expected = _exact_free_momentum_green(
                kx,
                ky,
                tau,
                beta,
                mu,
            )
            up = Matrix(@view snapshot.momentum_space[
                1, :, :, kx_index, ky_index, time_index,
            ])
            down = Matrix(@view snapshot.momentum_space[
                2, :, :, kx_index, ky_index, time_index,
            ])
            @test isapprox(up, expected; atol=UNEQUAL_TIME_ATOL, rtol=0.0)
            @test isapprox(down, expected; atol=UNEQUAL_TIME_ATOL, rtol=0.0)
            @test isapprox(up, down; atol=UNEQUAL_TIME_ATOL, rtol=0.0)
        end
    end

    @testset "Endpoints and fermionic boundary discontinuity" begin
        identity_orbitals = Matrix{ComplexF64}(I, 2, 2)
        for spin in 1:2, ky_index in 1:L, kx_index in 1:L
            green_zero = Matrix(@view snapshot.momentum_space[
                spin, :, :, kx_index, ky_index, 1,
            ])
            green_beta = Matrix(@view snapshot.momentum_space[
                spin, :, :, kx_index, ky_index, end,
            ])
            # With the documented no-leading-minus convention,
            # G(beta,0) = I-G(0,0).  This is the endpoint form of the
            # antiperiodic fermionic boundary condition across the tau=0 jump.
            @test isapprox(
                green_zero + green_beta,
                identity_orbitals;
                atol=UNEQUAL_TIME_ATOL,
                rtol=0.0,
            )
        end

        for spin in 1:2, destination_orbital in 1:2,
            source_orbital in 1:2, displacement_y in 0:(L - 1),
            displacement_x in 0:(L - 1)
            endpoint_sum = snapshot.real_space[
                spin,
                destination_orbital,
                source_orbital,
                displacement_x + 1,
                displacement_y + 1,
                1,
            ] + snapshot.real_space[
                spin,
                destination_orbital,
                source_orbital,
                displacement_x + 1,
                displacement_y + 1,
                end,
            ]
            expected = destination_orbital == source_orbital &&
                       displacement_x == 0 && displacement_y == 0 ? 1.0 : 0.0
            @test isapprox(
                endpoint_sum,
                expected;
                atol=UNEQUAL_TIME_ATOL,
                rtol=0.0,
            )
        end
    end

    @testset "Orbital C4 covariance" begin
        orbital_exchange = ComplexF64[0 1; 1 0]
        for spin in 1:2, time_index in eachindex(snapshot.tau),
            ky_index in 1:L, kx_index in 1:L
            green_kx_ky = Matrix(@view snapshot.momentum_space[
                spin, :, :, kx_index, ky_index, time_index,
            ])
            green_ky_kx = Matrix(@view snapshot.momentum_space[
                spin, :, :, ky_index, kx_index, time_index,
            ])
            @test isapprox(
                green_kx_ky,
                orbital_exchange * green_ky_kx * orbital_exchange;
                atol=UNEQUAL_TIME_ATOL,
                rtol=0.0,
            )
        end
    end

    @testset "Finite Fourier transform normalization" begin
        transformed =
            AltermagnetDQMCTrial.fourier_transform_unequal_time_greens(
                snapshot.real_space,
            )
        @test transformed.momenta == snapshot.momenta
        @test transformed.values == snapshot.momentum_space

        for spin in 1:2, time_index in eachindex(snapshot.tau),
            source_orbital in 1:2, destination_orbital in 1:2,
            displacement_y in 0:(L - 1), displacement_x in 0:(L - 1)
            reconstructed = 0.0im
            for (ky_index, ky) in pairs(snapshot.momenta),
                (kx_index, kx) in pairs(snapshot.momenta)
                reconstructed += cis(
                    kx * displacement_x + ky * displacement_y,
                ) * snapshot.momentum_space[
                    spin,
                    destination_orbital,
                    source_orbital,
                    kx_index,
                    ky_index,
                    time_index,
                ]
            end
            reconstructed /= L^2
            @test isapprox(
                reconstructed,
                snapshot.real_space[
                    spin,
                    destination_orbital,
                    source_orbital,
                    displacement_x + 1,
                    displacement_y + 1,
                    time_index,
                ];
                atol=UNEQUAL_TIME_ATOL,
                rtol=0.0,
            )
        end
    end

    @testset "Boundary and transform validation" begin
        post_sweep = initialize_dqmc_state(
            parameters;
            hs_fields=fixed_fields,
        )
        sweep_all_time_slices!(post_sweep)
        @test !post_sweep.greens_calculator_up.forward
        @test post_sweep.greens_calculator_up.l == dimensions.Ltau
        post_sweep_before = _serialized_bytes(post_sweep)
        post_sweep_snapshot =
            AltermagnetDQMCTrial.measure_unequal_time_greens(post_sweep)
        @test all(isfinite, post_sweep_snapshot.real_space)
        @test _serialized_bytes(post_sweep) == post_sweep_before

        mid_sweep = initialize_dqmc_state(parameters; hs_fields=fixed_fields)
        mid_sweep.greens_calculator_up.l = 2
        mid_sweep.greens_calculator_down.l = 2
        @test_throws ArgumentError begin
            AltermagnetDQMCTrial.measure_unequal_time_greens(mid_sweep)
        end

        unsynchronized = initialize_dqmc_state(
            parameters;
            hs_fields=fixed_fields,
        )
        unsynchronized.greens_calculator_down.forward = false
        @test_throws ArgumentError begin
            AltermagnetDQMCTrial.measure_unequal_time_greens(unsynchronized)
        end

        @test_throws DimensionMismatch begin
            AltermagnetDQMCTrial.fourier_transform_unequal_time_greens(
                zeros(1, 2, 2, L, L, 1),
            )
        end
        @test_throws DimensionMismatch begin
            AltermagnetDQMCTrial.fourier_transform_unequal_time_greens(
                zeros(2, 2, 2, L, L + 1, 1),
            )
        end
        nonfinite = zeros(2, 2, 2, L, L, 1)
        nonfinite[1] = NaN
        @test_throws ArgumentError begin
            AltermagnetDQMCTrial.fourier_transform_unequal_time_greens(
                nonfinite,
            )
        end
    end
end

const UNEQUAL_TIME_SUMMARY_TAU = [0.0, 0.5, 1.0]
const UNEQUAL_TIME_SUMMARY_MOMENTA = [0.0, pi]

function _unequal_time_summary_snapshot(
    configuration_sign,
    channel;
    tau=UNEQUAL_TIME_SUMMARY_TAU,
    momenta=UNEQUAL_TIME_SUMMARY_MOMENTA,
)
    floating_tau = Float64.(tau)
    floating_momenta = Float64.(momenta)
    L = length(floating_momenta)
    number_of_times = length(floating_tau)
    shape = (2, 2, 2, L, L, number_of_times)
    real_space = zeros(Float64, shape)
    momentum_space = zeros(ComplexF64, shape)
    momentum_space[2, 1, 2, 2, 1, :] .= channel
    return AltermagnetDQMCTrial.UnequalTimeGreensSnapshot(
        floating_tau,
        floating_momenta,
        real_space,
        momentum_space,
        Float64(configuration_sign),
        (0.0, 0.0),
    )
end

function _summarize_test_unequal_time_channel(
    snapshots;
    block_size=2,
    minimum_absolute_average_sign=0.0,
)
    return AltermagnetDQMCTrial.
           summarize_sign_reweighted_unequal_time_channel(
        snapshots;
        spin=2,
        destination_orbital=1,
        source_orbital=2,
        kx_index=2,
        ky_index=1,
        block_size=block_size,
        minimum_absolute_average_sign=minimum_absolute_average_sign,
    )
end

function _independent_unequal_time_channel_statistics(
    signs,
    channels,
    block_size,
)
    measurement_count = length(signs)
    block_ranges = [
        first_index:min(first_index + block_size - 1, measurement_count)
        for first_index in 1:block_size:measurement_count
    ]
    number_of_times = length(first(channels))
    block_signs = [sum(signs[block_range]) for block_range in block_ranges]
    block_numerators = zeros(ComplexF64, length(block_ranges), number_of_times)
    for (block, block_range) in enumerate(block_ranges)
        for measurement in block_range, time_index in 1:number_of_times
            block_numerators[block, time_index] +=
                signs[measurement] * channels[measurement][time_index]
        end
    end
    full_sign = sum(block_signs)
    full_numerator = vec(sum(block_numerators; dims=1))
    leave_one_out = Matrix{ComplexF64}(
        undef,
        length(block_ranges),
        number_of_times,
    )
    for block in eachindex(block_ranges)
        leave_one_out[block, :] .= (
            full_numerator .- view(block_numerators, block, :)
        ) ./ (full_sign - block_signs[block])
    end
    concatenated = hcat(real.(leave_one_out), imag.(leave_one_out))
    center = vec(sum(concatenated; dims=1)) ./ length(block_ranges)
    covariance = zeros(Float64, 2 * number_of_times, 2 * number_of_times)
    for block in eachindex(block_ranges)
        deviation = view(concatenated, block, :) .- center
        covariance .+= deviation * transpose(deviation)
    end
    covariance .*= (length(block_ranges) - 1) / length(block_ranges)
    return (
        block_ranges=block_ranges,
        block_signs=block_signs,
        block_numerators=block_numerators,
        full_sign=full_sign,
        full_numerator=full_numerator,
        green=full_numerator ./ full_sign,
        covariance=covariance,
    )
end

@testset "Sign-reweighted unequal-time channel summary" begin
    signs = [1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0]
    channels = [
        ComplexF64[
            0.15 * measurement - 0.4 +
                im * (0.03 * measurement^2 - 0.2),
            -0.25 + 0.08 * measurement +
                im * (0.35 - 0.04 * measurement),
            0.6 - 0.05 * measurement +
                im * (0.02 * measurement^2 - 0.1),
        ] for measurement in eachindex(signs)
    ]
    snapshots = [
        _unequal_time_summary_snapshot(signs[index], channels[index])
        for index in eachindex(signs)
    ]

    @testset "Exactly-once weighting and complex covariance" begin
        snapshots_before = _serialized_bytes(snapshots)
        summary = _summarize_test_unequal_time_channel(snapshots)
        independent = _independent_unequal_time_channel_statistics(
            signs,
            channels,
            2,
        )
        @test _serialized_bytes(snapshots) == snapshots_before
        @test summary.measurement_count == 8
        @test summary.requested_block_size == 2
        @test summary.block_count == 4
        @test summary.block_lengths == [2, 2, 2, 2]
        @test summary.minimum_block_length == 2
        @test summary.maximum_block_length == 2
        @test summary.spin == 2
        @test summary.destination_orbital == 1
        @test summary.source_orbital == 2
        @test summary.kx_index == 2
        @test summary.ky_index == 1
        @test summary.kx == Float64(pi)
        @test summary.ky == 0.0
        @test summary.tau == UNEQUAL_TIME_SUMMARY_TAU
        @test summary.block_sign_denominators == independent.block_signs
        @test summary.block_average_signs == [1.0, 0.0, 1.0, 0.0]
        @test summary.block_signed_numerators ≈ independent.block_numerators
        @test summary.full_sign_denominator == independent.full_sign == 4.0
        @test summary.full_average_sign == 0.5
        @test summary.full_signed_numerator ≈ independent.full_numerator
        @test summary.full_denominator_usable
        @test summary.leave_one_out_sign_denominators == [2.0, 4.0, 2.0, 4.0]
        @test summary.leave_one_out_average_signs ==
              Union{Missing,Float64}[1 / 3, 2 / 3, 1 / 3, 2 / 3]
        @test all(summary.leave_one_out_denominator_usable)
        @test summary.leave_one_out_denominators_usable
        @test summary.leave_one_out_denominator_failure_count == 0
        @test summary.minimum_absolute_leave_one_out_average_sign == 1 / 3
        @test summary.denominator_guard_passed
        @test summary.channel_values_finite
        @test summary.central_estimate_valid
        @test summary.jackknife_valid
        @test summary.covariance_ready
        @test summary.green ≈ independent.green
        @test summary.covariance ≈ independent.covariance
        @test summary.covariance ≈ transpose(summary.covariance)
        expected_standard_errors = sqrt.(diag(independent.covariance))
        @test summary.real_standard_errors ≈ expected_standard_errors[1:3]
        @test summary.imag_standard_errors ≈ expected_standard_errors[4:6]

        # Multiplying by the embedded sign twice would produce this different
        # numerator, because s^2=1.
        double_weighted = vec(sum(reduce(hcat, channels); dims=2)) ./
                          sum(signs)
        @test !isapprox(summary.green, double_weighted; atol=1.0e-12)
    end

    @testset "Common denominator failures" begin
        zero_signs = [1.0, -1.0, 1.0, -1.0]
        zero_snapshots = [
            _unequal_time_summary_snapshot(
                zero_signs[index],
                channels[index],
            ) for index in eachindex(zero_signs)
        ]
        zero_summary = _summarize_test_unequal_time_channel(zero_snapshots)
        @test zero_summary.full_sign_denominator == 0.0
        @test !zero_summary.full_denominator_usable
        @test !zero_summary.central_estimate_valid
        @test !zero_summary.jackknife_valid
        @test !zero_summary.covariance_ready
        @test ismissing(zero_summary.green)
        @test ismissing(zero_summary.real_standard_errors)
        @test ismissing(zero_summary.imag_standard_errors)
        @test ismissing(zero_summary.covariance)

        guarded_signs = [1.0, 1.0, -1.0, -1.0, 1.0, 1.0]
        guarded_snapshots = [
            _unequal_time_summary_snapshot(
                guarded_signs[index],
                channels[index],
            ) for index in eachindex(guarded_signs)
        ]
        guarded = _summarize_test_unequal_time_channel(guarded_snapshots)
        @test guarded.full_sign_denominator == 2.0
        @test guarded.full_denominator_usable
        @test guarded.central_estimate_valid
        @test !guarded.leave_one_out_denominators_usable
        @test guarded.leave_one_out_denominator_failure_count == 2
        @test guarded.minimum_absolute_leave_one_out_average_sign == 0.0
        @test !guarded.denominator_guard_passed
        @test !guarded.jackknife_valid
        @test !guarded.covariance_ready
        @test !ismissing(guarded.green)
        @test ismissing(guarded.real_standard_errors)
        @test ismissing(guarded.imag_standard_errors)
        @test ismissing(guarded.covariance)

        thresholded = _summarize_test_unequal_time_channel(
            snapshots;
            minimum_absolute_average_sign=1 / 3,
        )
        @test thresholded.full_denominator_usable
        @test !thresholded.leave_one_out_denominators_usable
        @test thresholded.leave_one_out_denominator_failure_count == 2
        @test !thresholded.jackknife_valid
        @test ismissing(thresholded.covariance)
    end

    @testset "Final short block" begin
        short_snapshots = [
            _unequal_time_summary_snapshot(1.0, channels[index])
            for index in 1:7
        ]
        summary = _summarize_test_unequal_time_channel(
            short_snapshots;
            block_size=3,
        )
        @test summary.block_count == 3
        @test summary.block_lengths == [3, 3, 1]
        @test summary.minimum_block_length == 1
        @test summary.maximum_block_length == 3
        @test summary.block_sign_denominators == [3.0, 3.0, 1.0]
        @test summary.block_signed_numerators[3, :] == channels[7]
        @test summary.green ≈ vec(sum(reduce(hcat, channels[1:7]); dims=2)) ./ 7
        @test summary.jackknife_valid
        @test summary.covariance_ready
    end

    @testset "Grid, channel, and sign validation" begin
        mismatched_tau = copy(snapshots)
        mismatched_tau[2] = _unequal_time_summary_snapshot(
            signs[2],
            channels[2];
            tau=[0.0, 0.4, 1.0],
        )
        @test_throws ArgumentError begin
            _summarize_test_unequal_time_channel(mismatched_tau)
        end

        mismatched_momenta = copy(snapshots)
        mismatched_momenta[2] = _unequal_time_summary_snapshot(
            signs[2],
            channels[2];
            momenta=[0.0, pi / 2],
        )
        @test_throws ArgumentError begin
            _summarize_test_unequal_time_channel(mismatched_momenta)
        end

        bad_shape = AltermagnetDQMCTrial.UnequalTimeGreensSnapshot(
            copy(UNEQUAL_TIME_SUMMARY_TAU),
            copy(UNEQUAL_TIME_SUMMARY_MOMENTA),
            zeros(Float64, 2, 2, 2, 2, 2, 3),
            zeros(ComplexF64, 2, 2, 2, 2, 1, 3),
            1.0,
            (0.0, 0.0),
        )
        @test_throws DimensionMismatch begin
            _summarize_test_unequal_time_channel([bad_shape])
        end

        for invalid_sign in (0.0, 0.5, NaN, Inf)
            invalid_snapshot = _unequal_time_summary_snapshot(
                invalid_sign,
                channels[1],
            )
            @test_throws ArgumentError begin
                _summarize_test_unequal_time_channel([invalid_snapshot])
            end
        end

        @test_throws ArgumentError begin
            _summarize_test_unequal_time_channel(
                AltermagnetDQMCTrial.UnequalTimeGreensSnapshot[],
            )
        end
        @test_throws ArgumentError begin
            _summarize_test_unequal_time_channel(snapshots; block_size=0)
        end
        @test_throws ArgumentError begin
            _summarize_test_unequal_time_channel(
                snapshots;
                minimum_absolute_average_sign=1.0,
            )
        end

        invalid_channel_keywords = (
            (spin=0,),
            (destination_orbital=3,),
            (source_orbital=0,),
            (kx_index=3,),
            (ky_index=0,),
        )
        for keywords in invalid_channel_keywords
            @test_throws ArgumentError begin
                AltermagnetDQMCTrial.
                summarize_sign_reweighted_unequal_time_channel(
                    snapshots;
                    spin=get(keywords, :spin, 2),
                    destination_orbital=get(
                        keywords,
                        :destination_orbital,
                        1,
                    ),
                    source_orbital=get(keywords, :source_orbital, 2),
                    kx_index=get(keywords, :kx_index, 2),
                    ky_index=get(keywords, :ky_index, 1),
                    block_size=2,
                )
            end
        end
    end
end

const FREE_GREEN_IDENTITY_ATOL = 2.0e-12
const FREE_PROPAGATOR_ATOL = 2.0e-11
const FREE_DENSITY_ATOL = 1.0e-11

# The Green tolerance covers symmetric-eigensolver reconstruction and dense
# multiplication roundoff. The propagator tolerance is looser because powers
# of dense slice matrices accumulate roundoff. The density tolerance enforces
# the requested 1e-11 real-space/momentum-space agreement independently.

free_maximum_abs(matrix) = maximum(abs, matrix)

function stable_test_fermi(scaled_energy)
    if scaled_energy >= 0
        tail = exp(-scaled_energy)
        return tail / (1 + tail)
    else
        return 1 / (1 + exp(scaled_energy))
    end
end

function momentum_space_filling(L, beta, mu; parameters=DEFAULT_PARAMETERS)
    occupation = 0.0
    for ky in allowed_momenta(L), kx in allowed_momenta(L)
        for energy in band_energies(kx, ky, parameters)
            occupation += stable_test_fermi(beta * (energy - mu))
        end
    end
    return 2 * occupation / L^2
end

function sampled_finite_temperature_filling(samples, beta, mu)
    occupation = sum(
        stable_test_fermi(beta * (energy - mu)) for energy in samples.energies
    )
    return 2 * occupation / samples.Nk^2
end

@testset "Milestone 4 finite-temperature free fermions" begin
    @testset "API checks, type behavior, and chemical-potential shift" begin
        kinetic = [1.0 0.25; 0.25 -2.0]
        original = copy(kinetic)
        shifted = shifted_kinetic_matrix(kinetic, 0.75)
        @test shifted == [0.25 0.25; 0.25 -2.75]
        @test kinetic == original
        @test issymmetric(shifted)

        rational_kinetic = [1 // 1 1 // 4; 1 // 4 -2 // 1]
        rational_shifted = shifted_kinetic_matrix(rational_kinetic, 1 // 2)
        @test eltype(rational_shifted) === Rational{Int}
        @test rational_shifted == [1 // 2 1 // 4; 1 // 4 -5 // 2]

        float32_kinetic = Float32[1 0.25; 0.25 -2]
        float32_shifted = shifted_kinetic_matrix(float32_kinetic, Float32(0.5))
        @test eltype(float32_shifted) === Float32
        @test eltype(free_equal_time_greens(float32_shifted, Float32(2))) ===
              Float32
        @test eltype(free_density_matrix(float32_shifted, Float32(2))) ===
              Float32

        @test_throws DimensionMismatch shifted_kinetic_matrix(zeros(2, 3), 0.0)
        @test_throws ArgumentError shifted_kinetic_matrix(
            [0.0 1.0; 0.0 0.0],
            0.0,
        )
        @test_throws ArgumentError shifted_kinetic_matrix([NaN 0.0; 0.0 0.0], 0.0)
        @test_throws ArgumentError shifted_kinetic_matrix(zeros(2, 2), Inf)
        @test_throws ArgumentError free_equal_time_greens(zeros(2, 2), -1.0)
        @test_throws ArgumentError free_density_matrix(zeros(2, 2), Inf)
        @test_throws ArgumentError dense_slice_propagators(zeros(2, 2), -0.1)
        @test_throws ArgumentError free_grand_log_partition(
            zeros(2, 2),
            1.0;
            spin_degeneracy=0,
        )
        @test_throws OverflowError dense_slice_propagators(
            [1000.0 0.0; 0.0 -1000.0],
            1.0,
        )
        @test_throws DimensionMismatch orbital_resolved_density(zeros(2, 2), 2)
        @test_throws DimensionMismatch total_filling_per_physical_site(
            zeros(2, 2),
            2,
        )

        identity_2 = Matrix{Float64}(I, 2, 2)
        B_zero, B_zero_inv = dense_slice_propagators(kinetic, 0.0)
        @test B_zero == identity_2
        @test B_zero_inv == identity_2
    end

    @testset "Green and density matrix conventions" begin
        for L in (2, 3, 4)
            kinetic = realspace_kinetic_matrix(L)
            dimension = 2 * L^2
            identity_matrix = Matrix{Float64}(I, dimension, dimension)
            for beta in (0.0, 0.5, 2.0, 8.0),
                mu in (-5.5, -1.897, 0.0, 9.5)
                kinetic_mu = shifted_kinetic_matrix(kinetic, mu)
                greens = free_equal_time_greens(kinetic_mu, beta)
                density = free_density_matrix(kinetic_mu, beta)
                greens_eigenvalues = eigvals(Symmetric(greens))
                density_eigenvalues = eigvals(Symmetric(density))

                @test free_maximum_abs(greens - transpose(greens)) <=
                      FREE_GREEN_IDENTITY_ATOL
                @test free_maximum_abs(density - transpose(density)) <=
                      FREE_GREEN_IDENTITY_ATOL
                @test free_maximum_abs(greens + density - identity_matrix) <=
                      FREE_GREEN_IDENTITY_ATOL
                @test minimum(greens_eigenvalues) >= -FREE_GREEN_IDENTITY_ATOL
                @test maximum(greens_eigenvalues) <= 1 + FREE_GREEN_IDENTITY_ATOL
                @test minimum(density_eigenvalues) >= -FREE_GREEN_IDENTITY_ATOL
                @test maximum(density_eigenvalues) <= 1 + FREE_GREEN_IDENTITY_ATOL
            end
        end
    end

    @testset "Dense slice propagators and direct-safe Green validation" begin
        cases = (
            (0.1, 4, 0.0),
            (0.2, 5, -1.897),
            (0.1, 20, -3.837),
        )
        for L in 2:4
            kinetic = realspace_kinetic_matrix(L)
            dimension = 2 * L^2
            identity_matrix = Matrix{Float64}(I, dimension, dimension)
            for (dtau, Ltau, mu) in cases
                beta = Ltau * dtau
                kinetic_mu = shifted_kinetic_matrix(kinetic, mu)
                B, B_inv = dense_slice_propagators(kinetic_mu, dtau)
                product = B^Ltau
                reference_product = exp(-beta * kinetic_mu)
                product_scale = max(1.0, free_maximum_abs(reference_product))

                @test free_maximum_abs(B * B_inv - identity_matrix) <=
                      FREE_PROPAGATOR_ATOL
                @test free_maximum_abs(B_inv * B - identity_matrix) <=
                      FREE_PROPAGATOR_ATOL
                @test free_maximum_abs(product - reference_product) /
                      product_scale <= FREE_PROPAGATOR_ATOL
                @test isapprox(
                    B,
                    exp(-dtau * kinetic_mu);
                    atol=FREE_PROPAGATOR_ATOL,
                    rtol=FREE_PROPAGATOR_ATOL,
                )
                @test isapprox(
                    B_inv,
                    exp(dtau * kinetic_mu);
                    atol=FREE_PROPAGATOR_ATOL,
                    rtol=FREE_PROPAGATOR_ATOL,
                )

                # These small systems and beta values are well conditioned
                # enough for an independent definition-level inverse check.
                direct_greens = inv(identity_matrix + product)
                greens = free_equal_time_greens(kinetic_mu, beta)
                @test free_maximum_abs(greens - direct_greens) <=
                      FREE_PROPAGATOR_ATOL
            end
        end
    end

    @testset "Momentum filling and C4-related orbital densities" begin
        for L in 2:6
            kinetic = realspace_kinetic_matrix(L)
            for beta in (0.0, 1.0, 4.0, 8.0),
                mu in (-5.5, -3.837, -1.897, 0.0, 9.5)
                density = free_density_matrix(
                    shifted_kinetic_matrix(kinetic, mu),
                    beta,
                )
                realspace_filling = total_filling_per_physical_site(density, L)
                momentum_filling = momentum_space_filling(L, beta, mu)
                orbital_filling = orbital_resolved_density(density, L)

                @test 0 - FREE_DENSITY_ATOL <= realspace_filling <=
                      4 + FREE_DENSITY_ATOL
                @test abs(realspace_filling - momentum_filling) <
                      FREE_DENSITY_ATOL
                @test abs(orbital_filling.n_x - orbital_filling.n_y) <
                      FREE_DENSITY_ATOL
                @test abs(
                    realspace_filling -
                    (orbital_filling.n_x + orbital_filling.n_y),
                ) < FREE_DENSITY_ATOL
            end
        end
    end

    @testset "Filling limits, monotonicity, and van Hove crossover" begin
        L = 4
        kinetic = realspace_kinetic_matrix(L)

        beta_zero_density = free_density_matrix(
            shifted_kinetic_matrix(kinetic, 3.7),
            0.0,
        )
        @test total_filling_per_physical_site(beta_zero_density, L) == 2.0

        empty_density = free_density_matrix(
            shifted_kinetic_matrix(kinetic, -100.0),
            2.0,
        )
        full_density = free_density_matrix(
            shifted_kinetic_matrix(kinetic, 100.0),
            2.0,
        )
        @test total_filling_per_physical_site(empty_density, L) < 1.0e-30
        @test abs(4 - total_filling_per_physical_site(full_density, L)) <
              FREE_GREEN_IDENTITY_ATOL

        chemical_potentials = collect(range(-20.0, 20.0; length=81))
        fillings = [
            total_filling_per_physical_site(
                free_density_matrix(
                    shifted_kinetic_matrix(kinetic, mu),
                    3.0,
                ),
                L,
            ) for mu in chemical_potentials
        ]
        @test all(diff(fillings) .>= -FREE_DENSITY_ATOL)
        @test first(fillings) < 1.0e-15
        @test last(fillings) > 4 - 1.0e-14

        samples = sample_band_energies(200)
        spectrum = sorted_band_energies(samples)
        vhs_cases = ((-3.837, 0.423), (-1.897, 1.460))
        for (mu, milestone_2_filling) in vhs_cases
            zero_temperature = zero_temperature_filling(
                samples,
                mu;
                sorted_energies=spectrum,
            )
            thermal_fillings = [
                sampled_finite_temperature_filling(samples, beta, mu) for
                beta in (2.0, 8.0, 32.0, 128.0)
            ]
            @test abs(last(thermal_fillings) - zero_temperature) <
                  abs(first(thermal_fillings) - zero_temperature)
            @test abs(last(thermal_fillings) - zero_temperature) < 5.0e-3
            @test abs(zero_temperature - milestone_2_filling) < 3.0e-3
        end
    end

    @testset "Grand partition function and positive free weight" begin
        for L in 2:4
            kinetic = realspace_kinetic_matrix(L)
            dimension = 2 * L^2
            for beta in (0.0, 1.0, 4.0), mu in (-3.837, 0.0, 9.5)
                kinetic_mu = shifted_kinetic_matrix(kinetic, mu)
                log_z_one_spin = free_grand_log_partition(kinetic_mu, beta)
                log_z_both_spins = free_grand_log_partition(
                    kinetic_mu,
                    beta;
                    spin_degeneracy=2,
                )
                @test isfinite(log_z_one_spin)
                @test isapprox(
                    log_z_both_spins,
                    2 * log_z_one_spin;
                    atol=FREE_GREEN_IDENTITY_ATOL,
                    rtol=FREE_GREEN_IDENTITY_ATOL,
                )
                if iszero(beta)
                    @test isapprox(
                        log_z_one_spin,
                        dimension * log(2);
                        atol=FREE_GREEN_IDENTITY_ATOL,
                        rtol=0.0,
                    )
                end
            end
        end

        # A direct determinant is used only for this numerically safe case.
        L = 2
        beta = 1.0
        dtau = 0.1
        Ltau = 10
        kinetic_mu = shifted_kinetic_matrix(realspace_kinetic_matrix(L), 0.0)
        B, _ = dense_slice_propagators(kinetic_mu, dtau)
        identity_matrix = Matrix{Float64}(I, 2 * L^2, 2 * L^2)
        determinant_logabs, determinant_sign = logabsdet(
            identity_matrix + B^Ltau,
        )
        log_z_one_spin = free_grand_log_partition(kinetic_mu, beta)
        log_z_both_spins = free_grand_log_partition(
            kinetic_mu,
            beta;
            spin_degeneracy=2,
        )
        @test determinant_sign > 0
        @test isapprox(
            determinant_logabs,
            log_z_one_spin;
            atol=FREE_PROPAGATOR_ATOL,
            rtol=FREE_PROPAGATOR_ATOL,
        )
        @test exp(log_z_one_spin) > 0
        @test exp(log_z_both_spins) > 0
    end
end

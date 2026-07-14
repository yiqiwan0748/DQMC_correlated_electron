using AltermagnetDQMCTrial
using LinearAlgebra
using Test

# These tolerances cover Float64 trigonometric and eigensolver roundoff while
# staying strict relative to the order-one energy scale of the model.
const TEST_ATOL = 1.0e-12
const TEST_RTOL = 1.0e-12

@testset "Noninteracting two-orbital Bloch model" begin
    parameters = ModelParameters()

    @testset "Parameters and numeric types" begin
        @test parameters isa ModelParameters{Float64}
        @test (parameters.t1, parameters.t2, parameters.t3, parameters.t4) ==
              (-1.0, -1.75, -0.85, -0.65)
        @test (DEFAULT_PARAMETERS.t1, DEFAULT_PARAMETERS.t2,
               DEFAULT_PARAMETERS.t3, DEFAULT_PARAMETERS.t4) ==
              (-1.0, -1.75, -0.85, -0.65)

        parameters32 = ModelParameters(
            t1=Float32(-1.0),
            t2=Float32(-1.75),
            t3=Float32(-0.85),
            t4=Float32(-0.65),
        )
        @test eltype(bloch_hamiltonian(Float32(0.2), Float32(-0.4), parameters32)) ===
              Float32
        @test eltype(band_energies(Float32(0.2), Float32(-0.4), parameters32)) ===
              Float32
    end

    @testset "Matrix structure" begin
        hamiltonian = bloch_hamiltonian(0.37, -1.11, parameters)
        @test size(hamiltonian) == (2, 2)
        @test eltype(hamiltonian) <: Real
        @test issymmetric(hamiltonian)
        @test ishermitian(hamiltonian)
    end

    @testset "High-symmetry points" begin
        cases = (
            ("Gamma", 0.0, 0.0, [8.9 0.0; 0.0 8.9], [8.9, 8.9]),
            ("X", pi, 0.0, [-1.9 0.0; 0.0 -4.9], [-4.9, -1.9]),
            ("Y", 0.0, pi, [-4.9 0.0; 0.0 -1.9], [-4.9, -1.9]),
            ("M", pi, pi, [-2.1 0.0; 0.0 -2.1], [-2.1, -2.1]),
        )

        for (name, kx, ky, expected_hamiltonian, expected_energies) in cases
            @testset "$name" begin
                @test isapprox(
                    bloch_hamiltonian(kx, ky, parameters),
                    expected_hamiltonian;
                    atol=TEST_ATOL,
                    rtol=TEST_RTOL,
                )
                @test isapprox(
                    band_energies(kx, ky, parameters),
                    expected_energies;
                    atol=TEST_ATOL,
                    rtol=TEST_RTOL,
                )
            end
        end
    end

    @testset "Generic momenta" begin
        cases = (
            (
                0.37,
                -1.11,
                [4.830508311237669 -0.842136014544862;
                 -0.842136014544862 5.562007054534160],
            ),
            (
                -2.2,
                0.73,
                [-0.059914118625784 -1.401820424464128;
                 -1.401820424464128 -2.060427398026108],
            ),
            (
                pi / 3,
                pi / 7,
                [5.685038113092580 0.976961284948442;
                 0.976961284948442 5.083584811238952],
            ),
        )

        for (kx, ky, expected_hamiltonian) in cases
            @test isapprox(
                bloch_hamiltonian(kx, ky, parameters),
                expected_hamiltonian;
                atol=TEST_ATOL,
                rtol=TEST_RTOL,
            )
        end
    end

    momenta = (
        (0.37, -1.11),
        (-2.2, 0.73),
        (pi / 3, pi / 7),
        (0.0, 0.61),
        (pi, 0.0),
        (pi, pi),
    )

    @testset "Orbital-exchange symmetry" begin
        for (kx, ky) in momenta
            h_kxky = bloch_hamiltonian(kx, ky, parameters)
            h_kykx = bloch_hamiltonian(ky, kx, parameters)

            # epsilon_x(kx, ky) = epsilon_y(ky, kx)
            @test isapprox(
                h_kxky[1, 1],
                h_kykx[2, 2];
                atol=TEST_ATOL,
                rtol=TEST_RTOL,
            )
            @test isapprox(
                band_energies(kx, ky, parameters),
                band_energies(ky, kx, parameters);
                atol=TEST_ATOL,
                rtol=TEST_RTOL,
            )
        end
    end

    @testset "Analytic eigenvalues" begin
        for (kx, ky) in momenta
            hamiltonian = bloch_hamiltonian(kx, ky, parameters)
            numerical = eigvals(Symmetric(hamiltonian))
            analytic = band_energies(kx, ky, parameters)
            @test isapprox(analytic, numerical; atol=TEST_ATOL, rtol=TEST_RTOL)
        end
    end
end

@testset "Milestone 2 noninteracting analysis" begin
    @testset "Band path and high-symmetry ticks" begin
        path = band_path(points_per_segment=8)
        @test path.tick_labels == ["Γ", "X", "M", "Γ"]
        @test path.tick_indices == [1, 9, 17, 25]
        @test isapprox(
            path.tick_positions,
            [0.0, pi, 2pi, 2pi + sqrt(2) * pi];
            atol=TEST_ATOL,
            rtol=TEST_RTOL,
        )
        @test (path.kx[1], path.ky[1]) == (0.0, 0.0)
        @test isapprox(
            [path.kx[9], path.ky[9]],
            [pi, 0.0];
            atol=TEST_ATOL,
            rtol=TEST_RTOL,
        )
        @test isapprox(
            [path.kx[17], path.ky[17]],
            [pi, pi];
            atol=TEST_ATOL,
            rtol=TEST_RTOL,
        )
        @test (path.kx[end], path.ky[end]) == (0.0, 0.0)
        @test all(diff(path.coordinate) .> 0)

        expected_endpoint_energies = (
            [8.9, 8.9],
            [-4.9, -1.9],
            [-2.1, -2.1],
            [8.9, 8.9],
        )
        for (index, expected) in zip(path.tick_indices, expected_endpoint_energies)
            actual = [path.lower_energy[index], path.upper_energy[index]]
            @test isapprox(actual, expected; atol=TEST_ATOL, rtol=TEST_RTOL)
        end
    end

    @testset "Grid DOS and filling" begin
        samples = sample_band_energies(64)
        spectrum = sorted_band_energies(samples)
        dos = density_of_states(samples; bins=400, smoothing_width=0.08)
        @test samples.Nk == 64
        @test size(samples.energies) == (2, 64^2)
        @test length(spectrum) == 2 * 64^2
        @test issorted(spectrum)

        # Normalization is explicitly restored after the truncated Gaussian
        # convolution; 2e-12 covers only Float64 summation roundoff.
        @test isapprox(
            sum(dos.density) * dos.bin_width,
            2.0;
            atol=2.0e-12,
            rtol=0.0,
        )

        energy_min, energy_max = extrema(samples.energies)
        chemical_potentials = collect(range(
            energy_min - 1,
            energy_max + 1;
            length=401,
        ))
        fillings = zero_temperature_filling(
            samples,
            chemical_potentials;
            sorted_energies=spectrum,
        )
        @test all(diff(fillings) .>= 0)
        @test first(fillings) == 0.0
        @test last(fillings) == 4.0
        @test zero_temperature_filling(
            samples,
            energy_min - 1;
            sorted_energies=spectrum,
        ) == 0.0
        @test zero_temperature_filling(
            samples,
            energy_max + 1;
            sorted_energies=spectrum,
        ) == 4.0

        peaks = estimate_van_hove_peaks(
            dos,
            samples;
            sorted_energies=spectrum,
        )
        @test length(peaks) == 2
        @test peaks[1].energy < peaks[2].energy
        @test all(0 <= peak.filling <= 4 for peak in peaks)
        @test all(peak.delta_c == (peak.filling - 2) / 4 for peak in peaks)
        @test infer_doping_convention(peaks).convention == :delta_c
    end

    @testset "Deterministic fixed-input output" begin
        first_samples = sample_band_energies(24)
        second_samples = sample_band_energies(24)
        @test first_samples.k_values == second_samples.k_values
        @test first_samples.energies == second_samples.energies
        first_dos = density_of_states(first_samples; bins=120, smoothing_width=0.05)
        second_dos = density_of_states(second_samples; bins=120, smoothing_width=0.05)
        @test first_dos.energy == second_dos.energy
        @test first_dos.counts == second_dos.counts
        @test first_dos.density == second_dos.density
    end
end

include("realspace_kinetic_tests.jl")
include("free_fermion_tests.jl")
include("two_orbital_dqmc_tests.jl")

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

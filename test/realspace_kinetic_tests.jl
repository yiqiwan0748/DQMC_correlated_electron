const REALSPACE_STRUCTURAL_ATOL = 10 * eps(Float64)
const REALSPACE_FOURIER_ATOL = 1.0e-10
const REALSPACE_SPECTRUM_ATOL = 1.0e-10

# The structural tolerance is at machine-precision scale. The Fourier and
# spectral tolerances allow accumulated Float64 phase and eigensolver roundoff
# while remaining ten orders of magnitude below the hopping-energy scale.

function orbital_block(kinetic, source_x, source_y, destination_x, destination_y, L)
    return [
        kinetic[
            single_particle_index(source_x, source_y, source_orbital, L),
            single_particle_index(
                destination_x,
                destination_y,
                destination_orbital,
                L,
            ),
        ] for source_orbital in 1:2, destination_orbital in 1:2
    ]
end

@testset "Milestone 3 finite-size real-space kinetic matrix" begin
    parameters = ModelParameters()

    @testset "Checked orbital-fast indexing and wrapping" begin
        for L in (1, 2, 3, 5)
            expected_index = 1
            for y in 1:L, x in 1:L, orbital in 1:2
                p = single_particle_index(x, y, orbital, L)
                @test p == expected_index
                @test single_particle_coordinates(p, L) == (x, y, orbital)
                expected_index += 1
            end
            @test periodic_coordinate(1, L) == 1
            @test periodic_coordinate(L, L) == L
            @test periodic_coordinate(L + 1, L) == 1
            @test periodic_coordinate(0, L) == L
            @test periodic_coordinate(-1, L) == periodic_coordinate(L - 1, L)
        end

        @test_throws ArgumentError single_particle_index(0, 1, 1, 3)
        @test_throws ArgumentError single_particle_index(4, 1, 1, 3)
        @test_throws ArgumentError single_particle_index(1, 0, 1, 3)
        @test_throws ArgumentError single_particle_index(1, 4, 1, 3)
        @test_throws ArgumentError single_particle_index(1, 1, 0, 3)
        @test_throws ArgumentError single_particle_index(1, 1, 3, 3)
        @test_throws ArgumentError single_particle_coordinates(0, 3)
        @test_throws ArgumentError single_particle_coordinates(19, 3)
        @test_throws ArgumentError periodic_coordinate(1, 0)
        @test_throws ArgumentError periodic_coordinate(1, -2)
        @test_throws ArgumentError single_particle_index(1, 1, 1, 0)
        @test_throws ArgumentError single_particle_coordinates(1, 0)
        @test_throws ArgumentError realspace_kinetic_matrix(0)
        @test_throws ArgumentError realspace_kinetic_matrix(-1)
        @test_throws ArgumentError allowed_momenta(0)
        @test_throws DimensionMismatch fourier_block(zeros(2, 2), 2, 0.0, 0.0)
        @test_throws DimensionMismatch realspace_band_energies(zeros(2, 3))
        @test_throws ArgumentError realspace_band_energies([0.0 1.0; 0.0 0.0])
    end

    @testset "Dimensions, numeric type, and Hermiticity" begin
        for L in 1:6
            kinetic = realspace_kinetic_matrix(L, parameters)
            @test kinetic isa Matrix{Float64}
            @test size(kinetic) == (2 * L^2, 2 * L^2)
            @test issymmetric(kinetic)
            @test ishermitian(kinetic)
            @test maximum(abs, kinetic - transpose(kinetic)) <=
                  REALSPACE_STRUCTURAL_ATOL
        end

        rational_parameters = ModelParameters(1 // 1, 2 // 1, 3 // 1, 4 // 1)
        @test eltype(realspace_kinetic_matrix(3, rational_parameters)) ===
              Rational{Int}
    end

    @testset "Translation invariance" begin
        L = 4
        kinetic = realspace_kinetic_matrix(L, parameters)
        @test all(
            orbital_block(
                kinetic,
                x,
                y,
                periodic_coordinate(x + dx, L),
                periodic_coordinate(y + dy, L),
                L,
            ) == orbital_block(
                kinetic,
                1,
                1,
                periodic_coordinate(1 + dx, L),
                periodic_coordinate(1 + dy, L),
                L,
            ) for y in 1:L for x in 1:L for dy in 0:(L - 1) for
            dx in 0:(L - 1)
        )
    end

    @testset "Isolated hopping limits" begin
        L = 5
        x, y = 2, 3
        p1 = single_particle_index(x, y, 1, L)
        p2 = single_particle_index(x, y, 2, L)

        only_t1 = ModelParameters(1.25, 0.0, 0.0, 0.0)
        kinetic_t1 = realspace_kinetic_matrix(L, only_t1)
        @test kinetic_t1[p1, single_particle_index(x + 1, y, 1, L)] == -1.25
        @test kinetic_t1[p1, single_particle_index(x, y + 1, 1, L)] == 0.0
        @test kinetic_t1[p2, single_particle_index(x + 1, y, 2, L)] == 0.0
        @test kinetic_t1[p2, single_particle_index(x, y + 1, 2, L)] == -1.25
        @test count(value -> !iszero(value), kinetic_t1) == 4 * L^2

        only_t2 = ModelParameters(0.0, 1.5, 0.0, 0.0)
        kinetic_t2 = realspace_kinetic_matrix(L, only_t2)
        @test kinetic_t2[p1, single_particle_index(x + 1, y, 1, L)] == 0.0
        @test kinetic_t2[p1, single_particle_index(x, y + 1, 1, L)] == -1.5
        @test kinetic_t2[p2, single_particle_index(x + 1, y, 2, L)] == -1.5
        @test kinetic_t2[p2, single_particle_index(x, y + 1, 2, L)] == 0.0
        @test count(value -> !iszero(value), kinetic_t2) == 4 * L^2

        only_t3 = ModelParameters(0.0, 0.0, 0.75, 0.0)
        kinetic_t3 = realspace_kinetic_matrix(L, only_t3)
        @test kinetic_t3[p1, single_particle_index(x + 1, y + 1, 1, L)] == -0.75
        @test kinetic_t3[p1, single_particle_index(x + 1, y - 1, 1, L)] == -0.75
        @test kinetic_t3[p2, single_particle_index(x + 1, y + 1, 2, L)] == -0.75
        @test kinetic_t3[p1, single_particle_index(x + 1, y + 1, 2, L)] == 0.0
        @test count(value -> !iszero(value), kinetic_t3) == 8 * L^2

        only_t4 = ModelParameters(0.0, 0.0, 0.0, 0.625)
        kinetic_t4 = realspace_kinetic_matrix(L, only_t4)
        @test kinetic_t4[p1, single_particle_index(x + 1, y + 1, 2, L)] == 0.625
        @test kinetic_t4[p1, single_particle_index(x + 1, y - 1, 2, L)] == -0.625
        @test kinetic_t4[p2, single_particle_index(x + 1, y + 1, 1, L)] == 0.625
        @test kinetic_t4[p1, single_particle_index(x + 1, y + 1, 1, L)] == 0.0
        @test count(value -> !iszero(value), kinetic_t4) == 8 * L^2

        isolated_parameters = (only_t1, only_t2, only_t3, only_t4)
        kx, ky = 2pi / L, 4pi / L
        for isolated in isolated_parameters
            isolated_kinetic = realspace_kinetic_matrix(L, isolated)
            mismatch = maximum(abs, fourier_block(
                isolated_kinetic,
                L,
                kx,
                ky,
            ) - bloch_hamiltonian(kx, ky, isolated))
            @test mismatch < REALSPACE_FOURIER_ATOL
        end
    end

    @testset "t4 diagonal sign and L=2 duplicate-bond accumulation" begin
        L = 4
        t4 = 0.7
        kinetic_t4 = realspace_kinetic_matrix(
            L,
            ModelParameters(0.0, 0.0, 0.0, t4),
        )
        source = single_particle_index(2, 2, 1, L)
        plus_destination = single_particle_index(3, 3, 2, L)
        minus_destination = single_particle_index(3, 1, 2, L)
        @test kinetic_t4[source, plus_destination] == t4
        @test kinetic_t4[source, minus_destination] == -t4
        @test kinetic_t4[source, plus_destination] ==
              -kinetic_t4[source, minus_destination]

        L = 2
        duplicate_parameters = ModelParameters(1.0, 2.0, 3.0, 4.0)
        kinetic = realspace_kinetic_matrix(L, duplicate_parameters)
        source_1 = single_particle_index(1, 1, 1, L)
        source_2 = single_particle_index(1, 1, 2, L)
        @test kinetic[source_1, single_particle_index(2, 1, 1, L)] == -2.0
        @test kinetic[source_2, single_particle_index(2, 1, 2, L)] == -4.0
        @test kinetic[source_1, single_particle_index(1, 2, 1, L)] == -4.0
        @test kinetic[source_2, single_particle_index(1, 2, 2, L)] == -2.0
        @test kinetic[source_1, single_particle_index(2, 2, 1, L)] == -12.0
        @test kinetic[source_1, single_particle_index(2, 2, 2, L)] == 0.0
        @test iszero(realspace_kinetic_matrix(
            L,
            ModelParameters(0.0, 0.0, 0.0, 1.0),
        ))
    end

    @testset "Fourier blocks reproduce the Bloch Hamiltonian" begin
        for L in 2:6
            kinetic = realspace_kinetic_matrix(L, parameters)
            maximum_mismatch = 0.0
            for ky in allowed_momenta(L), kx in allowed_momenta(L)
                realspace_block = fourier_block(kinetic, L, kx, ky)
                bloch_block = bloch_hamiltonian(kx, ky, parameters)
                mismatch = maximum(abs, realspace_block - bloch_block)
                maximum_mismatch = max(maximum_mismatch, mismatch)
                @test all(isapprox.(
                    realspace_block,
                    bloch_block;
                    atol=REALSPACE_FOURIER_ATOL,
                    rtol=0.0,
                ))
            end
            @test maximum_mismatch < REALSPACE_FOURIER_ATOL
        end
    end

    @testset "Complete spectra and validated global range" begin
        for L in 2:6
            kinetic = realspace_kinetic_matrix(L, parameters)
            realspace_spectrum = realspace_band_energies(kinetic)
            momentum_spectrum = sort!([
                energy for ky in allowed_momenta(L) for
                kx in allowed_momenta(L) for
                energy in band_energies(kx, ky, parameters)
            ])
            @test length(realspace_spectrum) == 2 * L^2
            @test issorted(realspace_spectrum)
            @test maximum(abs, realspace_spectrum - momentum_spectrum) <
                  REALSPACE_SPECTRUM_ATOL

            minimum_energy, maximum_energy = extrema(realspace_spectrum)
            @test minimum_energy >= -4.9 - REALSPACE_SPECTRUM_ATOL
            @test maximum_energy <= 8.9 + REALSPACE_SPECTRUM_ATOL
            @test isapprox(
                maximum_energy,
                8.9;
                atol=REALSPACE_SPECTRUM_ATOL,
                rtol=0.0,
            )
            if iseven(L)
                @test isapprox(
                    minimum_energy,
                    -4.9;
                    atol=REALSPACE_SPECTRUM_ATOL,
                    rtol=0.0,
                )
            end
        end
    end
end

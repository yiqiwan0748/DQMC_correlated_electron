using AltermagnetDQMCTrial
using LinearAlgebra
using Printf

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIRECTORY = joinpath(
    REPOSITORY_ROOT,
    "results",
    "realspace_kinetic",
)
const LATTICE_SIZES = 2:6
const VALIDATION_TOLERANCE = 1.0e-10
const SYMMETRY_TOLERANCE = 10 * eps(Float64)

function momentum_spectrum(L, parameters)
    return sort!([
        energy for ky in allowed_momenta(L) for
        kx in allowed_momenta(L) for
        energy in band_energies(kx, ky, parameters)
    ])
end

function validate_lattice_size(L, parameters)
    kinetic = realspace_kinetic_matrix(L, parameters)
    symmetry_error = maximum(abs, kinetic - transpose(kinetic))

    fourier_error = 0.0
    for ky in allowed_momenta(L), kx in allowed_momenta(L)
        mismatch = maximum(abs, fourier_block(kinetic, L, kx, ky) -
                                bloch_hamiltonian(kx, ky, parameters))
        fourier_error = max(fourier_error, mismatch)
    end

    realspace_spectrum = realspace_band_energies(kinetic)
    expected_spectrum = momentum_spectrum(L, parameters)
    spectrum_error = maximum(abs, realspace_spectrum - expected_spectrum)
    minimum_eigenvalue, maximum_eigenvalue = extrema(realspace_spectrum)

    symmetry_error <= SYMMETRY_TOLERANCE || error(
        "L=$L symmetry error $symmetry_error exceeds $SYMMETRY_TOLERANCE",
    )
    fourier_error < VALIDATION_TOLERANCE || error(
        "L=$L Fourier error $fourier_error exceeds $VALIDATION_TOLERANCE",
    )
    spectrum_error < VALIDATION_TOLERANCE || error(
        "L=$L spectrum error $spectrum_error exceeds $VALIDATION_TOLERANCE",
    )
    minimum_eigenvalue >= -4.9 - VALIDATION_TOLERANCE || error(
        "L=$L minimum eigenvalue lies below the validated global range",
    )
    maximum_eigenvalue <= 8.9 + VALIDATION_TOLERANCE || error(
        "L=$L maximum eigenvalue lies above the validated global range",
    )

    return (
        L=L,
        matrix_dimension=size(kinetic, 1),
        max_symmetry_error=symmetry_error,
        max_fourier_block_error=fourier_error,
        max_spectrum_error=spectrum_error,
        minimum_eigenvalue=minimum_eigenvalue,
        maximum_eigenvalue=maximum_eigenvalue,
    )
end

parameters = ModelParameters()
rows = [validate_lattice_size(L, parameters) for L in LATTICE_SIZES]
mkpath(OUTPUT_DIRECTORY)

open(joinpath(OUTPUT_DIRECTORY, "validation_summary.csv"), "w") do io
    println(
        io,
        "L,matrix_dimension,max_symmetry_error,max_fourier_block_error," *
        "max_spectrum_error,minimum_eigenvalue,maximum_eigenvalue",
    )
    for row in rows
        @printf(
            io,
            "%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g\n",
            row.L,
            row.matrix_dimension,
            row.max_symmetry_error,
            row.max_fourier_block_error,
            row.max_spectrum_error,
            row.minimum_eigenvalue,
            row.maximum_eigenvalue,
        )
    end
end

open(joinpath(OUTPUT_DIRECTORY, "metadata.txt"), "w") do io
    println(io, "milestone=3")
    println(io, "model=finite-size real-space kinetic matrix for the validated two-orbital noninteracting model")
    println(io, "coordinate_convention=x,y in 1,...,L; phase coordinates are (x-1,y-1); lattice spacing=1")
    println(io, "orbital_convention=orbital 1:d_xz; orbital 2:d_yz")
    println(io, "linear_index_convention=site=x+L*(y-1); p=2*(site-1)+orbital (orbital-fast, 1-based)")
    println(io, "spin_sector_dimension=2*L^2")
    println(io, "hoppings=t1:-1.0,t2:-1.75,t3:-0.85,t4:-0.65")
    println(io, "positive_displacement_block_(+1,0)=[-t1 0; 0 -t2]")
    println(io, "positive_displacement_block_(0,+1)=[-t2 0; 0 -t1]")
    println(io, "positive_displacement_block_(+1,+1)=[-t3 +t4; +t4 -t3]")
    println(io, "positive_displacement_block_(+1,-1)=[-t3 -t4; -t4 -t3]")
    println(io, "periodic_boundary_convention=destination coordinates wrap independently into 1,...,L; every positive-displacement block and its transpose are accumulated with +=")
    println(io, "duplicate_bond_convention=coincident periodic bonds, including L=2 aliases, retain every contribution")
    println(io, "allowed_momenta=kx,ky=2*pi*n/L, n=0,...,L-1")
    println(io, "fourier_transform_convention=H_ab(k)=L^-2 sum_(r,r') exp(-i*k.r) K_(r,a),(r',b) exp(+i*k.r')")
    println(io, "chemical_potential=absent; no diagonal chemical-potential shift is included")
    println(io, "validation_lattice_sizes=$(join(LATTICE_SIZES, ','))")
    println(io, "maximum_allowed_fourier_and_spectrum_error=$VALIDATION_TOLERANCE")
end

for row in rows
    @printf(
        "L=%d dimension=%d symmetry_error=%.16g Fourier_error=%.16g spectrum_error=%.16g range=[%.16g,%.16g]\n",
        row.L,
        row.matrix_dimension,
        row.max_symmetry_error,
        row.max_fourier_block_error,
        row.max_spectrum_error,
        row.minimum_eigenvalue,
        row.maximum_eigenvalue,
    )
end
println("output_directory=$OUTPUT_DIRECTORY")

using AltermagnetDQMCTrial
using LinearAlgebra
using Printf

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIRECTORY = joinpath(REPOSITORY_ROOT, "results", "free_fermion")
const LATTICE_SIZES = (2, 4, 6)
const BETAS = (1.0, 2.0, 4.0, 8.0)
const CHEMICAL_POTENTIALS = (-5.5, -3.837, -1.897, 0.0, 9.5)
const DTAU = 0.1

const GREEN_IDENTITY_TOLERANCE = 2.0e-12
const PROPAGATOR_INVERSE_TOLERANCE = 2.0e-11
const PROPAGATOR_PRODUCT_TOLERANCE = 5.0e-11
const DENSITY_COMPARISON_TOLERANCE = 1.0e-11

matrix_maximum_abs(matrix) = maximum(abs, matrix)

function stable_fermi(scaled_energy)
    if scaled_energy >= 0
        tail = exp(-scaled_energy)
        return tail / (1 + tail)
    else
        return 1 / (1 + exp(scaled_energy))
    end
end

function momentum_filling(L, beta, mu, parameters)
    one_spin_occupation = 0.0
    for ky in allowed_momenta(L), kx in allowed_momenta(L)
        for energy in band_energies(kx, ky, parameters)
            one_spin_occupation += stable_fermi(beta * (energy - mu))
        end
    end
    return 2 * one_spin_occupation / L^2
end

function relative_matrix_error(actual, reference)
    scale = max(1.0, matrix_maximum_abs(reference))
    return matrix_maximum_abs(actual - reference) / scale
end

function validate_point(L, beta, mu, kinetic, parameters)
    kinetic_mu = shifted_kinetic_matrix(kinetic, mu)
    dimension = size(kinetic_mu, 1)
    identity_matrix = Matrix{Float64}(I, dimension, dimension)

    greens = free_equal_time_greens(kinetic_mu, beta)
    density = free_density_matrix(kinetic_mu, beta)
    orbital_density = orbital_resolved_density(density, L)
    filling_realspace = total_filling_per_physical_site(density, L)
    filling_momentum = momentum_filling(L, beta, mu, parameters)

    greens_symmetry_error = matrix_maximum_abs(greens - transpose(greens))
    density_symmetry_error = matrix_maximum_abs(density - transpose(density))
    identity_error = matrix_maximum_abs(greens + density - identity_matrix)
    greens_eigenvalues = eigvals(Symmetric(greens))
    density_eigenvalues = eigvals(Symmetric(density))

    Ltau = round(Int, beta / DTAU)
    isapprox(Ltau * DTAU, beta; atol=10 * eps(Float64), rtol=0.0) ||
        error("beta/dtau must be integral")
    B, B_inv = dense_slice_propagators(kinetic_mu, DTAU)
    propagator_inverse_error = max(
        matrix_maximum_abs(B * B_inv - identity_matrix),
        matrix_maximum_abs(B_inv * B - identity_matrix),
    )
    accumulated_propagator = B^Ltau
    beta_propagator = exp(-beta * kinetic_mu)
    propagator_product_error = relative_matrix_error(
        accumulated_propagator,
        beta_propagator,
    )
    filling_mismatch = abs(filling_realspace - filling_momentum)

    log_z_one_spin = free_grand_log_partition(kinetic_mu, beta)
    log_z_both_spins = free_grand_log_partition(
        kinetic_mu,
        beta;
        spin_degeneracy=2,
    )

    direct_greens_error = 0.0
    determinant_sign = 1.0
    determinant_log_error = 0.0
    if beta == 1.0 && mu == 0.0
        direct_greens = inv(identity_matrix + accumulated_propagator)
        direct_greens_error = matrix_maximum_abs(greens - direct_greens)
        determinant_logabs, determinant_sign = logabsdet(
            identity_matrix + accumulated_propagator,
        )
        determinant_log_error = abs(determinant_logabs - log_z_one_spin)
    end

    greens_symmetry_error <= GREEN_IDENTITY_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu Green symmetry validation failed",
    )
    density_symmetry_error <= GREEN_IDENTITY_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu density symmetry validation failed",
    )
    identity_error <= GREEN_IDENTITY_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu G+C identity validation failed",
    )
    minimum(greens_eigenvalues) >= -GREEN_IDENTITY_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu has a Green eigenvalue below zero",
    )
    maximum(greens_eigenvalues) <= 1 + GREEN_IDENTITY_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu has a Green eigenvalue above one",
    )
    minimum(density_eigenvalues) >= -GREEN_IDENTITY_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu has a density eigenvalue below zero",
    )
    maximum(density_eigenvalues) <= 1 + GREEN_IDENTITY_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu has a density eigenvalue above one",
    )
    propagator_inverse_error <= PROPAGATOR_INVERSE_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu slice inverse validation failed",
    )
    propagator_product_error <= PROPAGATOR_PRODUCT_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu accumulated propagator validation failed",
    )
    filling_mismatch < DENSITY_COMPARISON_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu filling validation failed",
    )
    -DENSITY_COMPARISON_TOLERANCE <= filling_realspace <=
        4 + DENSITY_COMPARISON_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu filling lies outside [0,4]",
    )
    abs(orbital_density.n_x - orbital_density.n_y) <
        DENSITY_COMPARISON_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu orbital densities differ",
    )
    abs(
        filling_realspace - (orbital_density.n_x + orbital_density.n_y),
    ) < DENSITY_COMPARISON_TOLERANCE || error(
        "L=$L beta=$beta mu=$mu orbital densities do not sum to total",
    )
    isfinite(log_z_one_spin) || error("one-spin log Z is not finite")
    isapprox(
        log_z_both_spins,
        2 * log_z_one_spin;
        atol=GREEN_IDENTITY_TOLERANCE,
        rtol=GREEN_IDENTITY_TOLERANCE,
    ) || error("two-spin log Z is not twice the one-spin value")
    direct_greens_error <= PROPAGATOR_PRODUCT_TOLERANCE || error(
        "definition-level direct Green validation failed",
    )
    determinant_sign > 0 || error("one-spin determinant factor is not positive")
    determinant_log_error <= PROPAGATOR_PRODUCT_TOLERANCE || error(
        "one-spin determinant log does not match the spectral log Z",
    )

    return (
        L=L,
        beta=beta,
        mu=mu,
        matrix_dimension=dimension,
        total_filling_realspace=filling_realspace,
        total_filling_momentum=filling_momentum,
        orbital_x_filling=orbital_density.n_x,
        orbital_y_filling=orbital_density.n_y,
        greens_symmetry_error=greens_symmetry_error,
        identity_error_G_plus_C=identity_error,
        propagator_inverse_error=propagator_inverse_error,
        propagator_product_error=propagator_product_error,
        filling_mismatch=filling_mismatch,
        minimum_G_eigenvalue=minimum(greens_eigenvalues),
        maximum_G_eigenvalue=maximum(greens_eigenvalues),
        direct_greens_error=direct_greens_error,
        determinant_log_error=determinant_log_error,
        log_z_one_spin=log_z_one_spin,
        log_z_both_spins=log_z_both_spins,
    )
end

parameters = ModelParameters()
rows = NamedTuple[]
for L in LATTICE_SIZES
    kinetic = realspace_kinetic_matrix(L, parameters)
    for beta in BETAS, mu in CHEMICAL_POTENTIALS
        push!(rows, validate_point(L, beta, mu, kinetic, parameters))
    end
end

mkpath(OUTPUT_DIRECTORY)
open(joinpath(OUTPUT_DIRECTORY, "validation_summary.csv"), "w") do io
    println(
        io,
        "L,beta,mu,matrix_dimension,total_filling_realspace," *
        "total_filling_momentum,orbital_x_filling,orbital_y_filling," *
        "greens_symmetry_error,identity_error_G_plus_C," *
        "propagator_inverse_error,propagator_product_error," *
        "filling_mismatch,minimum_G_eigenvalue,maximum_G_eigenvalue",
    )
    for row in rows
        @printf(
            io,
            "%d,%.16g,%.16g,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n",
            row.L,
            row.beta,
            row.mu,
            row.matrix_dimension,
            row.total_filling_realspace,
            row.total_filling_momentum,
            row.orbital_x_filling,
            row.orbital_y_filling,
            row.greens_symmetry_error,
            row.identity_error_G_plus_C,
            row.propagator_inverse_error,
            row.propagator_product_error,
            row.filling_mismatch,
            row.minimum_G_eigenvalue,
            row.maximum_G_eigenvalue,
        )
    end
end

maximum_greens_symmetry_error = maximum(row.greens_symmetry_error for row in rows)
maximum_identity_error = maximum(row.identity_error_G_plus_C for row in rows)
maximum_propagator_inverse_error = maximum(
    row.propagator_inverse_error for row in rows
)
maximum_propagator_product_error = maximum(
    row.propagator_product_error for row in rows
)
maximum_filling_mismatch = maximum(row.filling_mismatch for row in rows)
maximum_orbital_mismatch = maximum(
    abs(row.orbital_x_filling - row.orbital_y_filling) for row in rows
)
maximum_direct_greens_error = maximum(row.direct_greens_error for row in rows)
maximum_determinant_log_error = maximum(
    row.determinant_log_error for row in rows
)

representative_lower_vhs = only(filter(
    row -> row.L == 6 && row.beta == 8.0 && row.mu == -3.837,
    rows,
))
representative_upper_vhs = only(filter(
    row -> row.L == 6 && row.beta == 8.0 && row.mu == -1.897,
    rows,
))
representative_partition = only(filter(
    row -> row.L == 4 && row.beta == 4.0 && row.mu == -1.897,
    rows,
))

open(joinpath(OUTPUT_DIRECTORY, "metadata.txt"), "w") do io
    println(io, "milestone=4")
    println(io, "model=finite-temperature free fermions for the validated two-orbital real-space Hamiltonian")
    println(io, "kinetic_convention=K_mu=K-mu*I; positive mu lowers every one-particle energy")
    println(io, "green_convention=G_ij=<c_i c_j^dagger>=(I+exp(-beta*K_mu))^-1")
    println(io, "density_convention=C_ij=<c_j^dagger c_i>; for this real-symmetric model C=I-G numerically")
    println(io, "index_convention=site=x+L*(y-1); p=2*(site-1)+orbital; orbital 1:d_xz; orbital 2:d_yz")
    println(io, "matrix_spin_convention=K_mu,G,C,B,and B_inv describe one spin sector")
    println(io, "spin_degeneracy=2 for every reported orbital and total filling")
    println(io, "orbital_normalization=n_x=2*sum_site C_(site,x),(site,x)/L^2 and likewise for n_y")
    println(io, "total_filling_normalization=n=2*tr(C)/L^2=n_x+n_y; physical range 0<=n<=4")
    println(io, "slice_propagators=B=exp(-dtau*K_mu); B_inv=exp(+dtau*K_mu)")
    println(io, "implementation=symmetric eigendecomposition; no direct inverse is used to construct G")
    println(io, "fermi_stability=logistic branches use exp(-abs(beta*(epsilon-mu))) to avoid overflow")
    println(io, "partition_convention=log(Z_one_spin)=sum_a log(1+exp(-beta*lambda_a)); log(Z_two_spin)=2*log(Z_one_spin)")
    println(io, "partition_stability=log1p softplus form; the determinant and Z are not explicitly formed in the primary implementation")
    println(io, "free_weight_sign=each one-spin factor 1+exp(-beta*lambda_a) is positive, hence the two-spin free-fermion weight is positive")
    println(io, "interacting_sign_scope=this free result does not imply that the interacting doped model is sign-problem free")
    println(io, "validation_lattice_sizes=$(join(LATTICE_SIZES, ','))")
    println(io, "validation_betas=$(join(BETAS, ','))")
    println(io, "validation_chemical_potentials=$(join(CHEMICAL_POTENTIALS, ','))")
    println(io, "validation_dtau=$DTAU; Ltau=beta/dtau")
    println(io, "green_identity_tolerance=$GREEN_IDENTITY_TOLERANCE")
    println(io, "propagator_inverse_tolerance=$PROPAGATOR_INVERSE_TOLERANCE")
    println(io, "propagator_product_tolerance=$PROPAGATOR_PRODUCT_TOLERANCE")
    println(io, "density_comparison_tolerance=$DENSITY_COMPARISON_TOLERANCE")
    println(io, "propagator_inverse_error=max(maxabs(B*B_inv-I),maxabs(B_inv*B-I))")
    println(io, "propagator_product_error=maxabs(B^Ltau-exp(-beta*K_mu))/max(1,maxabs(exp(-beta*K_mu)))")
    println(io, "direct_inverse_scope=G is compared with inv(I+B^Ltau) only for beta=1,mu=0 safe validation points")
    println(io, "maximum_direct_safe_G_error=$maximum_direct_greens_error")
    println(io, "maximum_direct_safe_log_determinant_error=$maximum_determinant_log_error")
    println(io, "representative_log_Z_one_spin_L4_beta4_mu-1.897=$(representative_partition.log_z_one_spin)")
    println(io, "representative_log_Z_two_spin_L4_beta4_mu-1.897=$(representative_partition.log_z_both_spins)")
    println(io, "representative_lower_vHs_filling_L6_beta8_mu-3.837=$(representative_lower_vhs.total_filling_realspace)")
    println(io, "representative_upper_vHs_filling_L6_beta8_mu-1.897=$(representative_upper_vhs.total_filling_realspace)")
    println(io, "maximum_G_symmetry_error=$maximum_greens_symmetry_error")
    println(io, "maximum_G_plus_C_identity_error=$maximum_identity_error")
    println(io, "maximum_propagator_inverse_error=$maximum_propagator_inverse_error")
    println(io, "maximum_propagator_product_error=$maximum_propagator_product_error")
    println(io, "maximum_filling_mismatch=$maximum_filling_mismatch")
    println(io, "maximum_orbital_filling_mismatch=$maximum_orbital_mismatch")
end

println("validated_points=$(length(rows))")
@printf("maximum_G_symmetry_error=%.16g\n", maximum_greens_symmetry_error)
@printf("maximum_G_plus_C_identity_error=%.16g\n", maximum_identity_error)
@printf("maximum_propagator_inverse_error=%.16g\n", maximum_propagator_inverse_error)
@printf("maximum_propagator_product_error=%.16g\n", maximum_propagator_product_error)
@printf("maximum_filling_mismatch=%.16g\n", maximum_filling_mismatch)
@printf("maximum_orbital_filling_mismatch=%.16g\n", maximum_orbital_mismatch)
@printf(
    "representative_lower_vHs_L6_beta8_mu-3.837=%.16g\n",
    representative_lower_vhs.total_filling_realspace,
)
@printf(
    "representative_upper_vHs_L6_beta8_mu-1.897=%.16g\n",
    representative_upper_vhs.total_filling_realspace,
)
println("one_spin_free_determinant_positive=true")
println("two_spin_free_weight_positive=true")
println("output_directory=$OUTPUT_DIRECTORY")

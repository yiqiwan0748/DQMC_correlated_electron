using LinearAlgebra: I, Symmetric, eigen, issymmetric

function _checked_real_symmetric_matrix(matrix::AbstractMatrix{<:Real}, name)
    size(matrix, 1) == size(matrix, 2) ||
        throw(DimensionMismatch("$name must be square"))
    issymmetric(matrix) || throw(ArgumentError("$name must be real symmetric"))
    all(isfinite, matrix) ||
        throw(ArgumentError("$name must contain only finite entries"))
    return size(matrix, 1)
end

function _checked_nonnegative_finite(value::Real, name)
    isfinite(value) || throw(ArgumentError("$name must be finite"))
    value >= zero(value) || throw(ArgumentError("$name must be nonnegative"))
    return value
end

function _floating_symmetric_matrix(
    matrix::AbstractMatrix{<:Real},
    scale::Real,
    matrix_name,
    scale_name,
)
    dimension = _checked_real_symmetric_matrix(matrix, matrix_name)
    _checked_nonnegative_finite(scale, scale_name)
    T = promote_type(float(eltype(matrix)), typeof(float(scale)))
    return Matrix{T}(matrix), convert(T, scale), dimension
end

function _spectral_reconstruction(eigenvectors, spectral_values)
    scaled_vectors = eigenvectors .* reshape(spectral_values, 1, :)
    reconstructed = scaled_vectors * transpose(eigenvectors)
    # Averaging removes last-bit asymmetry from dense multiplication while
    # leaving the spectral result unchanged to working precision.
    return (reconstructed + transpose(reconstructed)) / 2
end

function _fermi_green_and_density(x::T) where {T<:AbstractFloat}
    if x >= zero(T)
        tail = exp(-x)
        denominator = one(T) + tail
        return one(T) / denominator, tail / denominator
    else
        tail = exp(x)
        denominator = one(T) + tail
        return tail / denominator, one(T) / denominator
    end
end

function _free_correlation_matrix(
    kinetic_mu::AbstractMatrix{<:Real},
    beta::Real,
    return_density::Bool,
)
    floating_kinetic, floating_beta, dimension = _floating_symmetric_matrix(
        kinetic_mu,
        beta,
        "K_mu",
        "beta",
    )
    T = eltype(floating_kinetic)
    if iszero(floating_beta)
        half_identity = Matrix{T}(I, dimension, dimension)
        half_identity ./= convert(T, 2)
        return half_identity
    end

    decomposition = eigen(Symmetric(floating_kinetic))
    spectral_values = similar(decomposition.values)
    for index in eachindex(decomposition.values)
        scaled_energy = floating_beta * decomposition.values[index]
        green_value, density_value = _fermi_green_and_density(scaled_energy)
        spectral_values[index] = return_density ? density_value : green_value
    end
    return _spectral_reconstruction(decomposition.vectors, spectral_values)
end

"""
    shifted_kinetic_matrix(K, mu)

Return the dense real-symmetric one-particle matrix `K_mu = K - mu*I`.
`K` is the validated hopping matrix for one spin sector and is not mutated.
"""
function shifted_kinetic_matrix(
    kinetic::AbstractMatrix{<:Real},
    chemical_potential::Real,
)
    dimension = _checked_real_symmetric_matrix(kinetic, "K")
    isfinite(chemical_potential) || throw(ArgumentError("mu must be finite"))
    T = promote_type(eltype(kinetic), typeof(chemical_potential))
    shifted = Matrix{T}(kinetic)
    mu = convert(T, chemical_potential)
    for index in 1:dimension
        shifted[index, index] -= mu
    end
    return shifted
end

"""
    dense_slice_propagators(K_mu, dtau)

Return `(B, B_inv)` with `B = exp(-dtau*K_mu)` and
`B_inv = exp(+dtau*K_mu)`. Both matrices are reconstructed from one symmetric
eigendecomposition. An `OverflowError` is raised if either exponential is not
representable in the floating-point type; ordinary simulation parameters stay
well inside this limit.
"""
function dense_slice_propagators(
    kinetic_mu::AbstractMatrix{<:Real},
    dtau::Real,
)
    floating_kinetic, floating_dtau, dimension = _floating_symmetric_matrix(
        kinetic_mu,
        dtau,
        "K_mu",
        "dtau",
    )
    T = eltype(floating_kinetic)
    if iszero(floating_dtau)
        identity_matrix = Matrix{T}(I, dimension, dimension)
        return identity_matrix, copy(identity_matrix)
    end

    decomposition = eigen(Symmetric(floating_kinetic))
    slice_exponents = floating_dtau .* decomposition.values
    exponent_limit = log(floatmax(T))
    if any(
        exponent -> !isfinite(exponent) || abs(exponent) > exponent_limit,
        slice_exponents,
    )
        throw(OverflowError(
            "dtau*K_mu has eigenvalues whose exponentials are not representable",
        ))
    end
    forward_values = exp.(-slice_exponents)
    inverse_values = exp.(slice_exponents)
    forward = _spectral_reconstruction(
        decomposition.vectors,
        forward_values,
    )
    inverse = _spectral_reconstruction(
        decomposition.vectors,
        inverse_values,
    )
    return forward, inverse
end

"""
    free_equal_time_greens(K_mu, beta)

Return the one-spin DQMC equal-time Green matrix
`G[i,j] = <c_i c_j^dagger> = (I + exp(-beta*K_mu))^(-1)`.
The implementation applies a branch-stable logistic function to the
eigenvalues of real-symmetric `K_mu`; it does not form or invert the
potentially ill-conditioned matrix in the displayed definition.
"""
free_equal_time_greens(
    kinetic_mu::AbstractMatrix{<:Real},
    beta::Real,
) = _free_correlation_matrix(kinetic_mu, beta, false)

"""
    free_density_matrix(K_mu, beta)

Return the one-spin density matrix `C[i,j] = <c_j^dagger c_i>`. For this
real-symmetric model, `C = I - G`; its spectral occupations are evaluated
directly with a branch-stable Fermi function so that very small occupations do
not suffer subtractive cancellation.
"""
free_density_matrix(
    kinetic_mu::AbstractMatrix{<:Real},
    beta::Real,
) = _free_correlation_matrix(kinetic_mu, beta, true)

function _checked_density_matrix_dimension(
    density::AbstractMatrix{<:Real},
    L::Integer,
)
    _, number_of_sites, single_particle_dimension = _lattice_dimensions(L)
    size(density) == (single_particle_dimension, single_particle_dimension) ||
        throw(DimensionMismatch("C must have size (2L^2, 2L^2)"))
    all(isfinite, density) ||
        throw(ArgumentError("C must contain only finite entries"))
    return number_of_sites, single_particle_dimension
end

"""
    orbital_resolved_density(C, L)

Return `(n_x, n_y)`, the `d_xz` and `d_yz` fillings per physical site,
including the factor of two for spin degeneracy. `C` is a one-spin density
matrix in the established orbital-fast real-space indexing.
"""
function orbital_resolved_density(
    density::AbstractMatrix{<:Real},
    L::Integer,
)
    number_of_sites, _ = _checked_density_matrix_dimension(density, L)
    occupation_x = zero(eltype(density))
    occupation_y = zero(eltype(density))
    for site in 0:(number_of_sites - 1)
        occupation_x += density[2 * site + 1, 2 * site + 1]
        occupation_y += density[2 * site + 2, 2 * site + 2]
    end
    spin_degeneracy = 2
    return (
        n_x=spin_degeneracy * occupation_x / number_of_sites,
        n_y=spin_degeneracy * occupation_y / number_of_sites,
    )
end

"""
    total_filling_per_physical_site(C, L)

Return the total filling per physical site, including spin degeneracy. The
normalization is `2*tr(C)/L^2`, so the physical range is zero to four.
"""
function total_filling_per_physical_site(
    density::AbstractMatrix{<:Real},
    L::Integer,
)
    number_of_sites, single_particle_dimension =
        _checked_density_matrix_dimension(density, L)
    one_spin_particles = zero(eltype(density))
    for index in 1:single_particle_dimension
        one_spin_particles += density[index, index]
    end
    return 2 * one_spin_particles / number_of_sites
end

function _softplus(x::T) where {T<:AbstractFloat}
    return max(x, zero(T)) + log1p(exp(-abs(x)))
end

"""
    free_grand_log_partition(K_mu, beta; spin_degeneracy=1)

Return `log(Z)` for the free grand-canonical model. For one spin,
`log(Z) = sum_a log(1 + exp(-beta*lambda_a))`, where `lambda_a` are the
eigenvalues of `K_mu`. Set `spin_degeneracy=2` for the two-spin result. The
logarithms use a stable softplus expression and are summed without forming the
possibly overflowing determinant or partition function.
"""
function free_grand_log_partition(
    kinetic_mu::AbstractMatrix{<:Real},
    beta::Real;
    spin_degeneracy::Integer=1,
)
    spin_degeneracy > 0 ||
        throw(ArgumentError("spin_degeneracy must be positive"))
    floating_kinetic, floating_beta, dimension = _floating_symmetric_matrix(
        kinetic_mu,
        beta,
        "K_mu",
        "beta",
    )
    T = eltype(floating_kinetic)
    if iszero(floating_beta)
        return convert(T, spin_degeneracy * dimension) * log(convert(T, 2))
    end

    decomposition = eigen(Symmetric(floating_kinetic))
    log_partition_one_spin = zero(T)
    for energy in decomposition.values
        log_partition_one_spin += _softplus(-floating_beta * energy)
    end
    return spin_degeneracy * log_partition_one_spin
end

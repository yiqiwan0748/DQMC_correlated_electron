"""
    ModelParameters{T<:Real}

Immutable hopping parameters for the two-orbital noninteracting model.
"""
struct ModelParameters{T<:Real}
    t1::T
    t2::T
    t3::T
    t4::T
end

function ModelParameters(t1::Real, t2::Real, t3::Real, t4::Real)
    promoted = promote(t1, t2, t3, t4)
    T = typeof(promoted[1])
    return ModelParameters{T}(promoted...)
end

"""
    ModelParameters(; t1=-1.0, t2=-1.75, t3=-0.85, t4=-0.65)

Construct hopping parameters, using the published `Float64` values by default.
"""
ModelParameters(; t1::Real=-1.0, t2::Real=-1.75, t3::Real=-0.85, t4::Real=-0.65) =
    ModelParameters(t1, t2, t3, t4)

const DEFAULT_PARAMETERS = ModelParameters()

"""
    bloch_hamiltonian(kx, ky, parameters=DEFAULT_PARAMETERS)

Return the `2 x 2` real symmetric Bloch Hamiltonian in the orbital basis
`(d_xz, d_yz)`.
"""
function bloch_hamiltonian(
    kx::Real,
    ky::Real,
    parameters::ModelParameters=DEFAULT_PARAMETERS,
)
    (; t1, t2, t3, t4) = parameters
    cos_kx = cos(kx)
    cos_ky = cos(ky)
    sin_kx = sin(kx)
    sin_ky = sin(ky)

    epsilon_x = -2 * t1 * cos_kx - 2 * t2 * cos_ky - 4 * t3 * cos_kx * cos_ky
    epsilon_y = -2 * t2 * cos_kx - 2 * t1 * cos_ky - 4 * t3 * cos_kx * cos_ky
    epsilon_xy = -4 * t4 * sin_kx * sin_ky

    return [epsilon_x epsilon_xy; epsilon_xy epsilon_y]
end

"""
    band_energies(kx, ky, parameters=DEFAULT_PARAMETERS)

Return the two eigenvalues of the Bloch Hamiltonian in ascending order.
"""
function band_energies(
    kx::Real,
    ky::Real,
    parameters::ModelParameters=DEFAULT_PARAMETERS,
)
    hamiltonian = bloch_hamiltonian(kx, ky, parameters)
    epsilon_x = hamiltonian[1, 1]
    epsilon_y = hamiltonian[2, 2]
    epsilon_xy = hamiltonian[1, 2]

    center = epsilon_x / 2 + epsilon_y / 2
    radius = hypot(epsilon_x / 2 - epsilon_y / 2, epsilon_xy)
    return [center - radius, center + radius]
end

using LinearAlgebra: Hermitian, eigvals, ishermitian

function _checked_lattice_size(L::Integer)
    L > 0 || throw(ArgumentError("L must be a positive integer"))
    lattice_size = try
        Int(L)
    catch
        throw(ArgumentError("L must fit in an Int"))
    end
    return lattice_size
end

function _lattice_dimensions(L::Integer)
    lattice_size = _checked_lattice_size(L)
    number_of_sites = try
        Base.checked_mul(lattice_size, lattice_size)
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("L^2 must fit in an Int"))
    end
    single_particle_dimension = try
        Base.checked_mul(2, number_of_sites)
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("2*L^2 must fit in an Int"))
    end
    return lattice_size, number_of_sites, single_particle_dimension
end

"""
    periodic_coordinate(coordinate, L)

Wrap an integer coordinate periodically into `1:L`. The lattice size must be
positive.
"""
function periodic_coordinate(coordinate::Integer, L::Integer)
    lattice_size = _checked_lattice_size(L)
    remainder = mod(coordinate, lattice_size)
    return iszero(remainder) ? lattice_size : Int(remainder)
end

"""
    single_particle_index(x, y, orbital, L)

Return the orbital-fast one-particle index for the checked coordinates
`x, y in 1:L` and `orbital in 1:2`:

```text
site = x + L*(y - 1)
p = 2*(site - 1) + orbital
```
"""
function single_particle_index(
    x::Integer,
    y::Integer,
    orbital::Integer,
    L::Integer,
)
    lattice_size, _, _ = _lattice_dimensions(L)
    1 <= x <= lattice_size ||
        throw(ArgumentError("x must lie in 1:L"))
    1 <= y <= lattice_size ||
        throw(ArgumentError("y must lie in 1:L"))
    1 <= orbital <= 2 ||
        throw(ArgumentError("orbital must be 1 (d_xz) or 2 (d_yz)"))

    site = Int(x) + lattice_size * (Int(y) - 1)
    return 2 * (site - 1) + Int(orbital)
end

"""
    single_particle_coordinates(p, L)

Invert [`single_particle_index`](@ref), returning `(x, y, orbital)` after
checking that `p` lies in `1:2L^2`.
"""
function single_particle_coordinates(p::Integer, L::Integer)
    lattice_size, _, single_particle_dimension = _lattice_dimensions(L)
    1 <= p <= single_particle_dimension ||
        throw(ArgumentError("p must lie in 1:2L^2"))

    zero_based_index = Int(p) - 1
    orbital = mod(zero_based_index, 2) + 1
    zero_based_site = div(zero_based_index, 2)
    x = mod(zero_based_site, lattice_size) + 1
    y = div(zero_based_site, lattice_size) + 1
    return (x, y, orbital)
end

function _positive_displacement_blocks(parameters::ModelParameters{T}) where {T<:Real}
    (; t1, t2, t3, t4) = parameters
    zero_hopping = zero(T)
    T_x = T[-t1 zero_hopping; zero_hopping -t2]
    T_y = T[-t2 zero_hopping; zero_hopping -t1]
    T_plus = T[-t3 t4; t4 -t3]
    T_minus = T[-t3 -t4; -t4 -t3]
    return (
        (1, 0, T_x),
        (0, 1, T_y),
        (1, 1, T_plus),
        (1, -1, T_minus),
    )
end

"""
    realspace_kinetic_matrix(L, parameters=DEFAULT_PARAMETERS)

Construct the dense `2L^2 x 2L^2` real-space kinetic matrix for one spin
sector on an `L x L` periodic lattice. For each source cell, the four positive
displacement blocks are accumulated into `K`, followed by their transposes for
the reverse hoppings. Accumulation is intentional: on small periodic lattices,
distinct directed bonds can map to the same matrix entry.

No chemical-potential shift is included.
"""
function realspace_kinetic_matrix(
    L::Integer,
    parameters::ModelParameters=DEFAULT_PARAMETERS,
)
    lattice_size, _, single_particle_dimension = _lattice_dimensions(L)
    T = typeof(parameters.t1)
    kinetic = zeros(T, single_particle_dimension, single_particle_dimension)

    for y in 1:lattice_size, x in 1:lattice_size
        for (dx, dy, block) in _positive_displacement_blocks(parameters)
            destination_x = periodic_coordinate(x + dx, lattice_size)
            destination_y = periodic_coordinate(y + dy, lattice_size)
            for source_orbital in 1:2, destination_orbital in 1:2
                source = single_particle_index(
                    x,
                    y,
                    source_orbital,
                    lattice_size,
                )
                destination = single_particle_index(
                    destination_x,
                    destination_y,
                    destination_orbital,
                    lattice_size,
                )
                hopping = block[source_orbital, destination_orbital]
                kinetic[source, destination] += hopping
                kinetic[destination, source] += hopping
            end
        end
    end

    return kinetic
end

"""
    allowed_momenta(L)

Return the periodic-lattice momenta `2*pi*n/L` for `n = 0, ..., L-1`.
"""
function allowed_momenta(L::Integer)
    lattice_size = _checked_lattice_size(L)
    momentum_step = 2 * pi / lattice_size
    return [n * momentum_step for n in 0:(lattice_size - 1)]
end

raw"""
    fourier_block(K, L, kx, ky)

Project a real-space matrix onto one allowed momentum using

```math
H_{\alpha\beta}(k) = \frac{1}{L^2}\sum_{r,r'}
e^{-i k\cdot r} K_{r\alpha,r'\beta} e^{i k\cdot r'}.
```

Cell coordinates in the phases are `(x-1, y-1)`; using `(x, y)` would give
the same block because the common origin shift cancels.
"""
function fourier_block(
    kinetic::AbstractMatrix{<:Number},
    L::Integer,
    kx::Real,
    ky::Real,
)
    lattice_size, _, single_particle_dimension = _lattice_dimensions(L)
    size(kinetic) == (single_particle_dimension, single_particle_dimension) ||
        throw(DimensionMismatch(
            "K must have size (2L^2, 2L^2) for the supplied L",
        ))

    first_phase = cis(kx * 0 + ky * 0)
    coefficient_type = promote_type(eltype(kinetic), typeof(first_phase))
    plane_waves = zeros(coefficient_type, single_particle_dimension, 2)
    for y in 1:lattice_size, x in 1:lattice_size
        phase = cis(kx * (x - 1) + ky * (y - 1)) / lattice_size
        for orbital in 1:2
            p = single_particle_index(x, y, orbital, lattice_size)
            plane_waves[p, orbital] = phase
        end
    end

    return adjoint(plane_waves) * kinetic * plane_waves
end

"""
    realspace_band_energies(K)

Return the complete eigenvalue spectrum of a Hermitian real-space kinetic
matrix in ascending order.
"""
function realspace_band_energies(kinetic::AbstractMatrix{<:Number})
    size(kinetic, 1) == size(kinetic, 2) ||
        throw(DimensionMismatch("K must be square"))
    ishermitian(kinetic) ||
        throw(ArgumentError("K must be Hermitian"))
    return sort!(eigvals(Hermitian(kinetic)))
end

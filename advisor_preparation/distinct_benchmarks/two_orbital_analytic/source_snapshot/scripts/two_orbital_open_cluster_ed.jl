#!/usr/bin/env julia

module TwoOrbitalOpenClusterED

using LinearAlgebra
import LinearAlgebra: mul!

# Keep the validated model definitions namespaced.  In particular, the open
# cluster constructor below consumes the same positive-displacement hopping
# blocks as the package's periodic real-space constructor.
module ValidatedKinetic
include(joinpath(@__DIR__, "..", "src", "noninteracting.jl"))
include(joinpath(@__DIR__, "..", "src", "realspace_kinetic.jl"))
end

const ModelParameters = ValidatedKinetic.ModelParameters
const DEFAULT_PARAMETERS = ValidatedKinetic.DEFAULT_PARAMETERS

export ModelParameters,
       DEFAULT_PARAMETERS,
       SectorBasis,
       SectorHamiltonian,
       open_cluster_kinetic_matrix,
       sector_basis,
       sector_hamiltonian,
       dense_sector_hamiltonian,
       ground_state_observables,
       canonical_spin_sector_energies,
       spin_flip_error,
       scan_u

"""Return an orbital-fast index for an `Lx`-by-`Ly` open cluster."""
function _open_index(x::Integer, y::Integer, orbital::Integer, Lx::Integer)
    return 2 * ((Int(x) - 1) + Int(Lx) * (Int(y) - 1)) + Int(orbital)
end

function _checked_extent(value::Integer, name::AbstractString)
    value isa Bool && throw(ArgumentError("$name must be a positive integer"))
    value > 0 || throw(ArgumentError("$name must be a positive integer"))
    return try
        Int(value)
    catch
        throw(ArgumentError("$name must fit in an Int"))
    end
end

"""
    open_cluster_kinetic_matrix(Lx=2, Ly=Lx, parameters=DEFAULT_PARAMETERS)

Construct the one-spin kinetic matrix for an `Lx`-by-`Ly` cluster with open
boundaries.  Spatial modes use the orbital-fast order `(d_xz, d_yz)` in every
cell.  Each physical bond is added once together with its transpose.

The bond blocks are exactly those used by `realspace_kinetic_matrix`: nearest
neighbors contain `t1,t2`, same-orbital diagonals contain `t3`, and the two
diagonal orientations carry opposite `t4` signs.  Unlike the periodic 2-by-2
cluster, open boundaries therefore retain nonzero interorbital hopping.
"""
function open_cluster_kinetic_matrix(
    Lx::Integer=2,
    Ly::Integer=Lx,
    parameters::ModelParameters=DEFAULT_PARAMETERS,
)
    nx = _checked_extent(Lx, "Lx")
    ny = _checked_extent(Ly, "Ly")
    number_of_cells = try
        Base.checked_mul(nx, ny)
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("Lx*Ly must fit in an Int"))
    end
    number_of_modes = try
        Base.checked_mul(2, number_of_cells)
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("2*Lx*Ly must fit in an Int"))
    end

    T = typeof(parameters.t1)
    kinetic = zeros(T, number_of_modes, number_of_modes)
    for y in 1:ny, x in 1:nx
        for (dx, dy, block) in ValidatedKinetic._positive_displacement_blocks(parameters)
            destination_x = x + dx
            destination_y = y + dy
            1 <= destination_x <= nx || continue
            1 <= destination_y <= ny || continue
            for source_orbital in 1:2, destination_orbital in 1:2
                source = _open_index(x, y, source_orbital, nx)
                destination = _open_index(
                    destination_x,
                    destination_y,
                    destination_orbital,
                    nx,
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
Basis for one fixed `(Nup,Ndown)` sector.

`up_states` and `down_states` are bit patterns on `nmodes` spatial orbitals.
The full basis order is lexicographic in those two arrays, with the down-spin
state varying fastest.
"""
struct SectorBasis
    nmodes::Int
    Nup::Int
    Ndown::Int
    up_states::Vector{UInt64}
    down_states::Vector{UInt64}
    up_indices::Dict{UInt64,Int}
    down_indices::Dict{UInt64,Int}
end

Base.length(basis::SectorBasis) = length(basis.up_states) * length(basis.down_states)
Base.size(basis::SectorBasis) = (length(basis),)

function _fixed_weight_states(nmodes::Int, particles::Int)
    states = UInt64[]
    function append_states!(first_mode::Int, remaining::Int, state::UInt64)
        if iszero(remaining)
            push!(states, state)
            return
        end
        last_mode = nmodes - remaining + 1
        for mode in first_mode:last_mode
            append_states!(
                mode + 1,
                remaining - 1,
                state | (UInt64(1) << (mode - 1)),
            )
        end
    end
    append_states!(1, particles, UInt64(0))
    return states
end

"""
    sector_basis(nmodes, Nup, Ndown)
    sector_basis(kinetic, Nup, Ndown)

Build only the requested fixed-spin-number basis.  Bit storage limits this
small-cluster implementation to at most 64 spatial modes per spin.
"""
function sector_basis(nmodes::Integer, Nup::Integer, Ndown::Integer)
    modes = _checked_extent(nmodes, "nmodes")
    modes <= 64 || throw(ArgumentError("nmodes must not exceed 64"))
    (Nup isa Bool || Ndown isa Bool) &&
        throw(ArgumentError("Nup and Ndown must be particle counts"))
    0 <= Nup <= modes || throw(ArgumentError("Nup must lie in 0:nmodes"))
    0 <= Ndown <= modes || throw(ArgumentError("Ndown must lie in 0:nmodes"))
    nup = Int(Nup)
    ndown = Int(Ndown)
    up_states = _fixed_weight_states(modes, nup)
    down_states = _fixed_weight_states(modes, ndown)
    return SectorBasis(
        modes,
        nup,
        ndown,
        up_states,
        down_states,
        Dict(state => index for (index, state) in enumerate(up_states)),
        Dict(state => index for (index, state) in enumerate(down_states)),
    )
end

function sector_basis(
    kinetic::AbstractMatrix,
    Nup::Integer,
    Ndown::Integer,
)
    size(kinetic, 1) == size(kinetic, 2) ||
        throw(DimensionMismatch("kinetic matrix must be square"))
    return sector_basis(size(kinetic, 1), Nup, Ndown)
end

@inline function _basis_state(basis::SectorBasis, index::Int)
    number_of_down_states = length(basis.down_states)
    up_index = div(index - 1, number_of_down_states) + 1
    down_index = mod(index - 1, number_of_down_states) + 1
    return basis.up_states[up_index], basis.down_states[down_index]
end

@inline function _basis_index(
    basis::SectorBasis,
    up_state::UInt64,
    down_state::UInt64,
)
    return (basis.up_indices[up_state] - 1) * length(basis.down_states) +
           basis.down_indices[down_state]
end

@inline _occupied(state::UInt64, mode::Int) =
    !iszero(state & (UInt64(1) << (mode - 1)))

@inline function _hop(state::UInt64, destination::Int, source::Int)
    source_bit = UInt64(1) << (source - 1)
    iszero(state & source_bit) && return nothing
    intermediate = state ⊻ source_bit
    destination_bit = UInt64(1) << (destination - 1)
    !iszero(intermediate & destination_bit) && return nothing
    parity = count_ones(state & (source_bit - UInt64(1))) +
             count_ones(intermediate & (destination_bit - UInt64(1)))
    return intermediate | destination_bit, (iseven(parity) ? 1 : -1)
end

"""
Matrix-free Hamiltonian in one fixed `(Nup,Ndown)` sector.

The diagonal arrays are also retained for inexpensive ground-state
observables.  Off-diagonal hopping is stored in compressed-column form; this
does not allocate the dense 4900-by-4900 half-filled 2-by-2 block.
"""
struct SectorHamiltonian{T<:AbstractFloat}
    basis::SectorBasis
    kinetic::Matrix{T}
    U::T
    h::T
    diagonal::Vector{T}
    oalm_values::Vector{T}
    double_occupancies::Vector{T}
    column_pointers::Vector{Int}
    row_indices::Vector{Int}
    hopping_values::Vector{T}
end

Base.size(hamiltonian::SectorHamiltonian) =
    (length(hamiltonian.basis), length(hamiltonian.basis))
Base.size(hamiltonian::SectorHamiltonian, dimension::Integer) = size(hamiltonian)[dimension]
Base.eltype(::Type{SectorHamiltonian{T}}) where {T} = T
Base.eltype(::SectorHamiltonian{T}) where {T} = T

function _checked_real(value::Real, name::AbstractString)
    isfinite(value) || throw(ArgumentError("$name must be finite"))
    return value
end

"""
    sector_hamiltonian(kinetic, U, Nup, Ndown; h=0)
    sector_hamiltonian(kinetic, Nup, Ndown; U=0, h=0)

Build the requested canonical block of

`K + U sum_a (n_a_up-1/2)(n_a_down-1/2) - h*Oalm`,

where

`Oalm = 1/2 sum_cell [(n_x_up-n_x_down)-(n_y_up-n_y_down)]`.

The returned object supports `size`, multiplication by a vector, and explicit
conversion with [`dense_sector_hamiltonian`](@ref).
"""
function sector_hamiltonian(
    kinetic::AbstractMatrix{<:Real},
    U::Real,
    Nup::Integer,
    Ndown::Integer;
    h::Real=0,
)
    size(kinetic, 1) == size(kinetic, 2) ||
        throw(DimensionMismatch("kinetic matrix must be square"))
    nmodes = size(kinetic, 1)
    iseven(nmodes) || throw(ArgumentError(
        "the two-orbital kinetic dimension must be even",
    ))
    nmodes <= 64 || throw(ArgumentError("at most 64 spatial modes are supported"))
    issymmetric(kinetic) || throw(ArgumentError("kinetic matrix must be real symmetric"))
    interaction_input = _checked_real(U, "U")
    interaction_input >= 0 || throw(ArgumentError("U must be nonnegative"))
    _checked_real(h, "h")

    T = promote_type(
        typeof(float(zero(eltype(kinetic)))),
        typeof(float(interaction_input)),
        typeof(float(h)),
    )
    converted_kinetic = Matrix{T}(kinetic)
    interaction = convert(T, interaction_input)
    field = convert(T, h)
    basis = sector_basis(nmodes, Nup, Ndown)
    dimension = length(basis)
    diagonal = zeros(T, dimension)
    oalm_values = zeros(T, dimension)
    double_occupancies = zeros(T, dimension)
    column_pointers = Vector{Int}(undef, dimension + 1)
    row_indices = Int[]
    hopping_values = T[]
    capacity_estimate = dimension * min(4nmodes, nmodes^2)
    sizehint!(row_indices, capacity_estimate)
    sizehint!(hopping_values, capacity_estimate)

    half = inv(convert(T, 2))
    quarter = inv(convert(T, 4))
    for column in 1:dimension
        column_pointers[column] = length(row_indices) + 1
        up_state, down_state = _basis_state(basis, column)
        double_count = count_ones(up_state & down_state)
        double_occupancies[column] = convert(T, double_count)

        oalm_integer = 0
        one_body_diagonal = zero(T)
        for mode in 1:nmodes
            up_occupied = _occupied(up_state, mode)
            down_occupied = _occupied(down_state, mode)
            orbital_sign = isodd(mode) ? 1 : -1
            oalm_integer += orbital_sign * (up_occupied - down_occupied)
            if up_occupied
                one_body_diagonal += converted_kinetic[mode, mode]
            end
            if down_occupied
                one_body_diagonal += converted_kinetic[mode, mode]
            end
        end
        oalm = convert(T, oalm_integer) * half
        oalm_values[column] = oalm
        centered = interaction * (
            convert(T, double_count) -
            convert(T, basis.Nup + basis.Ndown) * half +
            convert(T, nmodes) * quarter
        )
        diagonal[column] = one_body_diagonal + centered - field * oalm

        for spin in 1:2
            state = spin == 1 ? up_state : down_state
            for source in 1:nmodes, destination in 1:nmodes
                source == destination && continue
                hopping = converted_kinetic[destination, source]
                iszero(hopping) && continue
                result = _hop(state, destination, source)
                result === nothing && continue
                final_spin_state, sign = result
                final_up = spin == 1 ? final_spin_state : up_state
                final_down = spin == 2 ? final_spin_state : down_state
                push!(row_indices, _basis_index(basis, final_up, final_down))
                push!(hopping_values, sign * hopping)
            end
        end
    end
    column_pointers[end] = length(row_indices) + 1

    return SectorHamiltonian(
        basis,
        converted_kinetic,
        interaction,
        field,
        diagonal,
        oalm_values,
        double_occupancies,
        column_pointers,
        row_indices,
        hopping_values,
    )
end

function sector_hamiltonian(
    kinetic::AbstractMatrix{<:Real},
    Nup::Integer,
    Ndown::Integer;
    U::Real=0,
    h::Real=0,
)
    return sector_hamiltonian(kinetic, U, Nup, Ndown; h=h)
end

function mul!(
    output::AbstractVector,
    hamiltonian::SectorHamiltonian,
    input::AbstractVector,
)
    length(input) == size(hamiltonian, 2) ||
        throw(DimensionMismatch("input vector has the wrong length"))
    length(output) == size(hamiltonian, 1) ||
        throw(DimensionMismatch("output vector has the wrong length"))
    source = Base.mightalias(output, input) ? copy(input) : input
    fill!(output, zero(eltype(output)))
    @inbounds for column in eachindex(hamiltonian.diagonal)
        coefficient = source[column]
        output[column] += hamiltonian.diagonal[column] * coefficient
        for pointer in hamiltonian.column_pointers[column]:(
            hamiltonian.column_pointers[column + 1] - 1
        )
            output[hamiltonian.row_indices[pointer]] +=
                hamiltonian.hopping_values[pointer] * coefficient
        end
    end
    return output
end

function Base.:*(hamiltonian::SectorHamiltonian, input::AbstractVector)
    T = promote_type(eltype(hamiltonian), eltype(input))
    output = zeros(T, size(hamiltonian, 1))
    return mul!(output, hamiltonian, input)
end

"""Materialize a sector operator as a dense symmetric matrix."""
function dense_sector_hamiltonian(hamiltonian::SectorHamiltonian{T}) where {T}
    dense = zeros(T, size(hamiltonian))
    for column in eachindex(hamiltonian.diagonal)
        dense[column, column] = hamiltonian.diagonal[column]
        for pointer in hamiltonian.column_pointers[column]:(
            hamiltonian.column_pointers[column + 1] - 1
        )
            dense[hamiltonian.row_indices[pointer], column] +=
                hamiltonian.hopping_values[pointer]
        end
    end
    return dense
end

Base.Matrix(hamiltonian::SectorHamiltonian) = dense_sector_hamiltonian(hamiltonian)

function _dense_ground(hamiltonian::SectorHamiltonian)
    decomposition = eigen(Symmetric(dense_sector_hamiltonian(hamiltonian)), 1:1)
    vector = copy(decomposition.vectors[:, 1])
    work = hamiltonian * vector
    energy = dot(vector, work)
    residual = norm(work - energy * vector)
    return (
        energy=energy,
        vector=vector,
        residual=residual,
        iterations=1,
        converged=true,
        method=:dense,
    )
end

function _deterministic_start(T::Type{<:AbstractFloat}, dimension::Int)
    vector = Vector{T}(undef, dimension)
    irrational_scale = sqrt(convert(T, 2))
    for index in eachindex(vector)
        argument = convert(T, index)
        vector[index] = sin(argument) + cos(irrational_scale * argument)
    end
    vector ./= norm(vector)
    return vector
end

function _operator_scale(hamiltonian::SectorHamiltonian{T}) where {T}
    bound = zero(T)
    for column in eachindex(hamiltonian.diagonal)
        column_sum = abs(hamiltonian.diagonal[column])
        for pointer in hamiltonian.column_pointers[column]:(
            hamiltonian.column_pointers[column + 1] - 1
        )
            column_sum += abs(hamiltonian.hopping_values[pointer])
        end
        bound = max(bound, column_sum)
    end
    return max(bound, one(T))
end

# Deterministic Lanczos with two full modified-Gram-Schmidt passes.  Keeping the
# Krylov vectors costs O(dimension*maxiter), rather than O(dimension^2), and is
# important here because the largest 2-by-2 canonical block has dimension 4900.
function _lanczos_ground(
    hamiltonian::SectorHamiltonian{T};
    atol::Real=1.0e-11,
    rtol::Real=1.0e-10,
    maxiter::Integer=min(size(hamiltonian, 1), 300),
    require_convergence::Bool=true,
) where {T}
    dimension = size(hamiltonian, 1)
    iterations_limit = min(_checked_extent(maxiter, "maxiter"), dimension)
    absolute_tolerance = convert(T, _checked_real(atol, "atol"))
    relative_tolerance = convert(T, _checked_real(rtol, "rtol"))
    absolute_tolerance >= zero(T) || throw(ArgumentError("atol must be nonnegative"))
    relative_tolerance >= zero(T) || throw(ArgumentError("rtol must be nonnegative"))

    krylov = Matrix{T}(undef, dimension, iterations_limit)
    diagonal = Vector{T}(undef, iterations_limit)
    off_diagonal = Vector{T}(undef, max(0, iterations_limit - 1))
    previous = zeros(T, dimension)
    current = _deterministic_start(T, dimension)
    work = similar(current)
    scale = _operator_scale(hamiltonian)
    best_value = zero(T)
    best_coefficients = ones(T, 1)
    completed = 0
    converged = false

    for iteration in 1:iterations_limit
        completed = iteration
        view(krylov, :, iteration) .= current
        mul!(work, hamiltonian, current)
        if iteration > 1
            work .-= off_diagonal[iteration - 1] .* previous
        end
        diagonal[iteration] = dot(current, work)
        work .-= diagonal[iteration] .* current

        # Full reorthogonalization, repeated once to control loss of
        # orthogonality near converged or degenerate Ritz values.
        active_basis = view(krylov, :, 1:iteration)
        for _ in 1:2
            coefficients = transpose(active_basis) * work
            work .-= active_basis * coefficients
        end
        next_norm = norm(work)

        tridiagonal = SymTridiagonal(
            diagonal[1:iteration],
            off_diagonal[1:max(0, iteration - 1)],
        )
        decomposition = eigen(tridiagonal, 1:1)
        best_value = decomposition.values[1]
        best_coefficients = copy(decomposition.vectors[:, 1])
        threshold = absolute_tolerance + relative_tolerance * max(abs(best_value), scale)
        ritz_residual = next_norm * abs(best_coefficients[end])
        if ritz_residual <= threshold || next_norm <= eps(T) * scale
            converged = true
            break
        end
        iteration == iterations_limit && break

        off_diagonal[iteration] = next_norm
        previous .= current
        current .= work ./ next_norm
    end

    vector = view(krylov, :, 1:completed) * best_coefficients
    vector ./= norm(vector)
    mul!(work, hamiltonian, vector)
    energy = dot(vector, work)
    residual = norm(work - energy * vector)
    threshold = absolute_tolerance + relative_tolerance * max(abs(energy), scale)
    converged = converged && residual <= 8threshold
    if require_convergence && !converged
        throw(ErrorException(
            "Lanczos ground solve did not converge in $completed iterations " *
            "(residual=$residual, threshold=$threshold)",
        ))
    end
    return (
        energy=energy,
        vector=vector,
        residual=residual,
        iterations=completed,
        converged=converged,
        method=:lanczos,
    )
end

"""
    ground_state_observables(hamiltonian; method=:auto, dense_cutoff=256, ...)
    ground_state_observables(kinetic, U, Nup, Ndown; h=0, ...)

Return the requested sector's ground eigensystem and diagonal observables.
`:auto` uses a dense solve only through `dense_cutoff`; larger blocks use the
deterministic, fully reorthogonalized matrix-free Lanczos implementation.

The observables belong to one normalized ground eigenvector. If the ground
space is degenerate, they can depend on the eigensolver's basis and are not a
ground-subspace trace. Energies remain basis-independent; callers that need a
symmetry-preserving ensemble must resolve or average the full ground space.
"""
function ground_state_observables(
    hamiltonian::SectorHamiltonian;
    method::Symbol=:auto,
    dense_cutoff::Integer=256,
    atol::Real=1.0e-11,
    rtol::Real=1.0e-10,
    maxiter::Integer=min(size(hamiltonian, 1), 300),
    require_convergence::Bool=true,
)
    cutoff = dense_cutoff < 0 ? throw(ArgumentError("dense_cutoff must be nonnegative")) :
             Int(dense_cutoff)
    selected_method = method == :auto ?
        (size(hamiltonian, 1) <= cutoff ? :dense : :lanczos) : method
    solution = if selected_method == :dense
        _dense_ground(hamiltonian)
    elseif selected_method == :lanczos
        _lanczos_ground(
            hamiltonian;
            atol=atol,
            rtol=rtol,
            maxiter=maxiter,
            require_convergence=require_convergence,
        )
    else
        throw(ArgumentError("method must be :auto, :dense, or :lanczos"))
    end

    probabilities = abs2.(solution.vector)
    probabilities ./= sum(probabilities)
    oalm = dot(probabilities, hamiltonian.oalm_values)
    oalm2 = dot(probabilities, abs2.(hamiltonian.oalm_values))
    double_occupancy = dot(probabilities, hamiltonian.double_occupancies)
    centered_interaction = hamiltonian.U * (
        double_occupancy -
        (hamiltonian.basis.Nup + hamiltonian.basis.Ndown) / 2 +
        hamiltonian.basis.nmodes / 4
    )
    field_energy = -hamiltonian.h * oalm
    variance = max(zero(oalm2), oalm2 - abs2(oalm))
    return (
        energy=solution.energy,
        eigenvector=solution.vector,
        Nup=hamiltonian.basis.Nup,
        Ndown=hamiltonian.basis.Ndown,
        N=hamiltonian.basis.Nup + hamiltonian.basis.Ndown,
        Sz=(hamiltonian.basis.Nup - hamiltonian.basis.Ndown) / 2,
        dimension=length(hamiltonian.basis),
        oalm=oalm,
        oalm2=oalm2,
        oalm_variance=variance,
        double_occupancy=double_occupancy,
        kinetic_energy=solution.energy - centered_interaction - field_energy,
        centered_interaction_energy=centered_interaction,
        field_energy=field_energy,
        residual=solution.residual,
        iterations=solution.iterations,
        converged=solution.converged,
        method=solution.method,
        ground_subspace_averaged=false,
    )
end

function ground_state_observables(
    kinetic::AbstractMatrix{<:Real},
    U::Real,
    Nup::Integer,
    Ndown::Integer;
    h::Real=0,
    kwargs...,
)
    hamiltonian = sector_hamiltonian(kinetic, U, Nup, Ndown; h=h)
    return ground_state_observables(hamiltonian; kwargs...)
end

function _allowed_spin_sectors(nmodes::Int, N::Int)
    minimum_up = max(0, N - nmodes)
    maximum_up = min(nmodes, N)
    return [(Nup, N - Nup) for Nup in minimum_up:maximum_up]
end

"""
    canonical_spin_sector_energies(kinetic, U, N; h=0, spin_sectors=nothing, ...)

Compute ground information only for the requested spin sectors at fixed total
particle number `N`.  With `spin_sectors=nothing`, all allowed `(Nup,Ndown)`
partitions of `N` are used.  Eigenvectors are intentionally omitted from the
returned summary rows. The diagonal observables are single-ground-vector
values and can be basis-dependent when a sector ground space is degenerate.
"""
function canonical_spin_sector_energies(
    kinetic::AbstractMatrix{<:Real},
    U::Real,
    N::Integer;
    h::Real=0,
    spin_sectors=nothing,
    kwargs...,
)
    nmodes = size(kinetic, 1)
    size(kinetic, 2) == nmodes || throw(DimensionMismatch("kinetic matrix must be square"))
    (N isa Bool) && throw(ArgumentError("N must be a particle count"))
    0 <= N <= 2nmodes || throw(ArgumentError("N must lie in 0:2nmodes"))
    particle_number = Int(N)
    requested = spin_sectors === nothing ?
        _allowed_spin_sectors(nmodes, particle_number) :
        [(Int(pair[1]), Int(pair[2])) for pair in spin_sectors]
    isempty(requested) && throw(ArgumentError("at least one spin sector is required"))
    length(unique(requested)) == length(requested) ||
        throw(ArgumentError("spin_sectors contains a duplicate"))

    rows = NamedTuple[]
    for (Nup, Ndown) in requested
        Nup + Ndown == particle_number || throw(ArgumentError(
            "each spin sector must sum to N",
        ))
        result = ground_state_observables(
            kinetic,
            U,
            Nup,
            Ndown;
            h=h,
            kwargs...,
        )
        push!(rows, (
            U=U,
            h=h,
            N=result.N,
            Nup=result.Nup,
            Ndown=result.Ndown,
            Sz=result.Sz,
            energy=result.energy,
            dimension=result.dimension,
            oalm=result.oalm,
            oalm2=result.oalm2,
            oalm_variance=result.oalm_variance,
            double_occupancy=result.double_occupancy,
            residual=result.residual,
            iterations=result.iterations,
            converged=result.converged,
            method=result.method,
            ground_subspace_averaged=result.ground_subspace_averaged,
        ))
    end
    return rows
end

"""
    spin_flip_error(rows)

Return `maximum(abs(E[Nup,Ndown]-E[Ndown,Nup]))`.  A complete set of
spin-flipped partners is required.  At `h=0` this is a numerical symmetry
check, not an order parameter.
"""
function spin_flip_error(rows)
    isempty(rows) && throw(ArgumentError("rows must not be empty"))
    energies = Dict((row.Nup, row.Ndown) => row.energy for row in rows)
    errors = Float64[]
    for row in rows
        partner = (row.Ndown, row.Nup)
        haskey(energies, partner) || throw(ArgumentError(
            "missing spin-flipped sector $partner",
        ))
        push!(errors, abs(Float64(row.energy - energies[partner])))
    end
    return maximum(errors)
end

function spin_flip_error(
    kinetic::AbstractMatrix{<:Real},
    U::Real,
    N::Integer;
    kwargs...,
)
    rows = canonical_spin_sector_energies(kinetic, U, N; h=0, kwargs...)
    return spin_flip_error(rows)
end

function _scan_u(
    kinetic::AbstractMatrix{<:Real},
    U_values;
    N::Integer=2,
    h::Real=0,
    spin_sectors=nothing,
    degeneracy_tolerance::Real=1.0e-9,
    kwargs...,
)
    interactions = collect(U_values)
    isempty(interactions) && throw(ArgumentError("U_values must not be empty"))
    all(value -> value isa Real && !(value isa Bool) && isfinite(value) && value >= 0,
        interactions) ||
        throw(ArgumentError(
            "all U_values must be finite, nonnegative real numbers",
        ))
    degeneracy_tolerance >= 0 && isfinite(degeneracy_tolerance) ||
        throw(ArgumentError("degeneracy_tolerance must be finite and nonnegative"))

    scan_rows = NamedTuple[]
    for U in interactions
        sectors = canonical_spin_sector_energies(
            kinetic,
            U,
            N;
            h=h,
            spin_sectors=spin_sectors,
            kwargs...,
        )
        ground_energy = minimum(row.energy for row in sectors)
        ground_sectors = [
            (Nup=row.Nup, Ndown=row.Ndown, Sz=row.Sz, oalm=row.oalm)
            for row in sectors
            if abs(row.energy - ground_energy) <= degeneracy_tolerance
        ]
        complete_pairs = all(
            row -> any(
                partner -> partner.Nup == row.Ndown && partner.Ndown == row.Nup,
                sectors,
            ),
            sectors,
        )
        pairing_error = iszero(h) && complete_pairs ? spin_flip_error(sectors) : missing
        push!(scan_rows, (
            U=U,
            N=Int(N),
            h=h,
            ground_energy=ground_energy,
            ground_sectors=ground_sectors,
            sector_energies=sectors,
            spin_flip_error=pairing_error,
        ))
    end
    return (
        rows=scan_rows,
        N=Int(N),
        h=h,
        kinetic=Matrix(kinetic),
        boundary=:open,
        note="Finite-cluster spin-flip pairing at h=0 is a symmetry check; " *
             "a finite cluster does not spontaneously select a magnetic domain.",
    )
end

"""
    scan_u(kinetic, U_values; N=2, ...)
    scan_u(U_values; Lx=2, Ly=Lx, parameters=DEFAULT_PARAMETERS, N=2, ...)

Scan canonical ground energies over `U`.  The default `N=2` keeps exploratory
scans inexpensive.  At half filling on 2-by-2 (`N=8`), the central dimension
is 4900 and `method=:auto` uses the matrix-free Lanczos solver.
"""
function scan_u(
    kinetic::AbstractMatrix{<:Real},
    U_values;
    kwargs...,
)
    return _scan_u(kinetic, U_values; kwargs...)
end

function scan_u(
    U_values;
    Lx::Integer=2,
    Ly::Integer=Lx,
    parameters::ModelParameters=DEFAULT_PARAMETERS,
    kwargs...,
)
    kinetic = open_cluster_kinetic_matrix(Lx, Ly, parameters)
    return _scan_u(kinetic, U_values; kwargs...)
end

end # module

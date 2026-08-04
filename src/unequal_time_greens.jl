import JDQMCFramework as jdqmcf

raw"""
    UnequalTimeGreensSnapshot{T}

One configuration's sign-unweighted, spin-resolved imaginary-time Green
function.  `real_space` and `momentum_space` both use the axes

```text
(spin, destination_orbital, source_orbital, x, y, tau_index)
```

where `spin = 1, 2` denotes up and down, the orbital axes use
`(d_xz, d_yz)`, and `tau[tau_index] = (tau_index-1)*dtau` spans the closed
grid `0:dtau:beta`.  In `real_space`, `x-1` and `y-1` are periodic cell
displacements.  In `momentum_space`, those axes index `momenta`, the allowed
values `2*pi*n/L` for `n=0,...,L-1`.

The convention has no leading minus sign:

```math
G_{ab,\sigma}(\mathbf r,\tau) = \frac{1}{L^2}\sum_{\mathbf i}
\langle c_{\mathbf i+\mathbf r,a,\sigma}(\tau)
c^\dagger_{\mathbf i,b,\sigma}(0)\rangle_{\mathcal C}.
```

`configuration_sign` is returned separately; neither array has been
multiplied by it.  `maximum_stabilization_error[spin]` reports the largest
correction made by JDQMCFramework while traversing imaginary time.
"""
struct UnequalTimeGreensSnapshot{T<:AbstractFloat}
    tau::Vector{T}
    momenta::Vector{T}
    real_space::Array{T,6}
    momentum_space::Array{Complex{T},6}
    configuration_sign::T
    maximum_stabilization_error::NTuple{2,T}
end

function _validate_unequal_time_measurement_boundary(state::DQMCState)
    up = state.greens_calculator_up
    down = state.greens_calculator_down
    up.Lτ == down.Lτ || throw(ArgumentError(
        "spin-sector Green calculators have different imaginary-time grids",
    ))
    up.l == down.l && up.forward == down.forward || throw(ArgumentError(
        "spin-sector Green calculators are not synchronized",
    ))
    complete_sweep_boundary =
        (up.forward && up.l == 1) || (!up.forward && up.l == up.Lτ)
    complete_sweep_boundary || throw(ArgumentError(
        "unequal-time measurement requires a complete-sweep boundary",
    ))
    return nothing
end

function _translation_average_unequal_time!(
    output::AbstractArray{T,4},
    green::AbstractMatrix{T},
    L::Int,
) where {T<:AbstractFloat}
    size(output) == (2, 2, L, L) || throw(DimensionMismatch(
        "a real-space unequal-time slice must have size (2, 2, L, L)",
    ))
    _, number_of_cells, single_particle_dimension = _lattice_dimensions(L)
    size(green) == (single_particle_dimension, single_particle_dimension) ||
        throw(DimensionMismatch(
            "an unequal-time Green matrix must have size (2L^2, 2L^2)",
        ))

    fill!(output, zero(T))
    inverse_number_of_cells = inv(convert(T, number_of_cells))
    for displacement_y in 0:(L - 1), displacement_x in 0:(L - 1)
        for source_y in 1:L, source_x in 1:L
            destination_x = periodic_coordinate(
                source_x + displacement_x,
                L,
            )
            destination_y = periodic_coordinate(
                source_y + displacement_y,
                L,
            )
            for source_orbital in 1:2, destination_orbital in 1:2
                destination = single_particle_index(
                    destination_x,
                    destination_y,
                    destination_orbital,
                    L,
                )
                source = single_particle_index(
                    source_x,
                    source_y,
                    source_orbital,
                    L,
                )
                output[
                    destination_orbital,
                    source_orbital,
                    displacement_x + 1,
                    displacement_y + 1,
                ] += green[destination, source] * inverse_number_of_cells
            end
        end
    end
    return output
end

function _measure_one_spin_unequal_time!(
    real_space_spin::AbstractArray{T,5},
    propagators,
    parameters::DQMCParameters{T},
) where {T<:AbstractFloat}
    (; Nsp, Ltau) = validate_dqmc_parameters(parameters)
    L = parameters.L
    size(real_space_spin) == (2, 2, L, L, Ltau + 1) ||
        throw(DimensionMismatch(
            "a one-spin result must have size (2, 2, L, L, Ltau+1)",
        ))

    # JDQMCFramework's calculator and its LDR workspace are mutable.  Build a
    # fresh calculator from private propagator copies rather than advancing a
    # calculator or touching a propagator owned by the sampling state.
    private_propagators = deepcopy(propagators)
    calculator = jdqmcf.FermionGreensCalculator(
        private_propagators,
        parameters.beta,
        parameters.dtau,
        parameters.stabilization_period,
    )
    calculator.forward && calculator.l == 1 || error(
        "the pinned JDQMCFramework calculator did not initialize at tau=0",
    )

    green_00 = zeros(T, Nsp, Nsp)
    logdet_green, signdet_green =
        jdqmcf.calculate_equaltime_greens!(green_00, calculator)
    green_tau_0 = similar(green_00)
    green_0_tau = similar(green_00)
    green_tau_tau = similar(green_00)
    jdqmcf.initialize_unequaltime_greens!(
        green_tau_0,
        green_0_tau,
        green_tau_tau,
        green_00,
    )
    _translation_average_unequal_time!(
        view(real_space_spin, :, :, :, :, 1),
        green_tau_0,
        L,
    )

    maximum_error = zero(T)
    visited = falses(Ltau)
    for time_slice in calculator
        jdqmcf.propagate_unequaltime_greens!(
            green_tau_0,
            green_0_tau,
            green_tau_tau,
            calculator,
            private_propagators,
        )
        logdet_green, signdet_green, correction, _ =
            jdqmcf.stabilize_unequaltime_greens!(
                green_tau_0,
                green_0_tau,
                green_tau_tau,
                logdet_green,
                signdet_green,
                calculator,
                private_propagators;
                update_B̄=false,
            )
        maximum_error = max(maximum_error, convert(T, correction))
        _translation_average_unequal_time!(
            view(real_space_spin, :, :, :, :, time_slice + 1),
            green_tau_0,
            L,
        )
        visited[time_slice] = true
    end
    all(visited) || error("the unequal-time traversal missed a time slice")
    all(isfinite, real_space_spin) || throw(OverflowError(
        "unequal-time propagation produced a non-finite Green function",
    ))
    return (
        signdet_green=convert(T, signdet_green),
        maximum_error=maximum_error,
    )
end

raw"""
    fourier_transform_unequal_time_greens(real_space)

Transform a real-space array with axes
`(spin, a, b, dx, dy, tau_index)` according to

```math
G_{ab,\sigma}(\mathbf k,\tau) = \sum_{\mathbf r}
e^{-i\mathbf k\cdot\mathbf r}G_{ab,\sigma}(\mathbf r,\tau).
```

Return `(momenta, values)`, where the fourth and fifth axes of `values` index
`kx` and `ky`.  There is no forward normalization; the inverse transform is
`G(r,tau) = sum_k exp(+i*k*r) G(k,tau) / L^2`.
"""
function fourier_transform_unequal_time_greens(
    real_space::AbstractArray{S,6},
) where {S<:Number}
    size(real_space, 1) == 2 || throw(DimensionMismatch(
        "the spin axis must have length two",
    ))
    size(real_space, 2) == 2 && size(real_space, 3) == 2 ||
        throw(DimensionMismatch("both orbital axes must have length two"))
    L = size(real_space, 4)
    L > 0 && size(real_space, 5) == L || throw(DimensionMismatch(
        "the displacement axes must describe a nonempty square lattice",
    ))
    all(isfinite, real_space) || throw(ArgumentError(
        "real_space must contain only finite entries",
    ))

    T = typeof(float(real(zero(S))))
    momenta = T.(allowed_momenta(L))
    number_of_times = size(real_space, 6)
    momentum_space = zeros(Complex{T}, 2, 2, 2, L, L, number_of_times)
    for (ky_index, ky) in pairs(momenta), (kx_index, kx) in pairs(momenta)
        for displacement_y in 0:(L - 1), displacement_x in 0:(L - 1)
            angle = kx * displacement_x + ky * displacement_y
            phase = Complex{T}(cos(angle), -sin(angle))
            for time_index in 1:number_of_times,
                source_orbital in 1:2,
                destination_orbital in 1:2,
                spin in 1:2
                momentum_space[
                    spin,
                    destination_orbital,
                    source_orbital,
                    kx_index,
                    ky_index,
                    time_index,
                ] += phase * real_space[
                    spin,
                    destination_orbital,
                    source_orbital,
                    displacement_x + 1,
                    displacement_y + 1,
                    time_index,
                ]
            end
        end
    end
    return (momenta=momenta, values=momentum_space)
end

raw"""
    measure_unequal_time_greens(state)

Measure one configuration's stabilized `G_ab_sigma(r,tau)` on the closed
imaginary-time grid `tau = 0:dtau:beta` and its finite-lattice momentum
transform.  The state must be at a synchronized complete-sweep boundary.

The routine is observational: it reconstructs each spin sector using private
propagator copies and a fresh JDQMCFramework calculator/workspace, and does not
advance or otherwise mutate `state`.  Returned Green functions are not sign
weighted; use `snapshot.configuration_sign` exactly once when accumulating a
Monte Carlo estimator.
"""
function measure_unequal_time_greens(
    state::DQMCState{T},
) where {T<:AbstractFloat}
    _validate_unequal_time_measurement_boundary(state)
    (; Ltau) = validate_dqmc_parameters(state.parameters)
    L = state.parameters.L
    real_space = Array{T}(undef, 2, 2, 2, L, L, Ltau + 1)
    up = _measure_one_spin_unequal_time!(
        view(real_space, 1, :, :, :, :, :),
        state.propagators_up,
        state.parameters,
    )
    down = _measure_one_spin_unequal_time!(
        view(real_space, 2, :, :, :, :, :),
        state.propagators_down,
        state.parameters,
    )
    transformed = fourier_transform_unequal_time_greens(real_space)
    tau = [convert(T, time_slice) * state.parameters.dtau
           for time_slice in 0:Ltau]
    configuration_sign = convert(T, sign(inv(up.signdet_green) *
                                           inv(down.signdet_green)))
    return UnequalTimeGreensSnapshot(
        tau,
        transformed.momenta,
        real_space,
        transformed.values,
        configuration_sign,
        (up.maximum_error, down.maximum_error),
    )
end

function _unequal_time_block_ranges(measurement_count::Int, block_size::Int)
    return [
        first_index:min(first_index + block_size - 1, measurement_count)
        for first_index in 1:block_size:measurement_count
    ]
end

function _checked_unequal_time_channel_index(
    value::Integer,
    upper::Int,
    name::AbstractString,
)
    index = try
        Int(value)
    catch error
        error isa InexactError || error isa OverflowError || rethrow()
        throw(ArgumentError("$name must fit in an Int"))
    end
    1 <= index <= upper || throw(ArgumentError(
        "$name must lie in 1:$upper",
    ))
    return index
end

function _validate_unequal_time_snapshot_layout(snapshot, measurement::Int)
    number_of_times = length(snapshot.tau)
    number_of_times > 0 || throw(ArgumentError(
        "snapshots[$measurement].tau must be nonempty",
    ))
    all(isfinite, snapshot.tau) || throw(ArgumentError(
        "snapshots[$measurement].tau must contain only finite entries",
    ))
    L = length(snapshot.momenta)
    L > 0 || throw(ArgumentError(
        "snapshots[$measurement].momenta must be nonempty",
    ))
    all(isfinite, snapshot.momenta) || throw(ArgumentError(
        "snapshots[$measurement].momenta must contain only finite entries",
    ))
    expected_shape = (2, 2, 2, L, L, number_of_times)
    size(snapshot.real_space) == expected_shape || throw(DimensionMismatch(
        "snapshots[$measurement].real_space must have size $expected_shape",
    ))
    size(snapshot.momentum_space) == expected_shape ||
        throw(DimensionMismatch(
            "snapshots[$measurement].momentum_space must have size " *
            "$expected_shape",
        ))
    return (L=L, number_of_times=number_of_times)
end

raw"""
    summarize_sign_reweighted_unequal_time_channel(
        snapshots;
        spin,
        destination_orbital,
        source_orbital,
        kx_index,
        ky_index,
        block_size,
        minimum_absolute_average_sign=0.0,
    )

Summarize one spin/orbital/momentum channel from configuration snapshots.
Each raw channel vector is multiplied by its embedded `configuration_sign`
exactly once.  The central estimate is the vector ratio

```math
G(\tau) = \frac{\sum_m s_m G_m(\tau)}{\sum_m s_m},
```

using one full sign denominator shared by every imaginary-time point.
Contiguous blocks retain a final short block.  Jackknife samples use the
common leave-one-block-out denominator for the entire complex time vector.

The returned `covariance` orders components as
`[real(G(tau)); imag(G(tau))]` and uses
`(B-1)/B * sum_b((G_b-Gbar)*(G_b-Gbar)')`.  Its diagonal yields
`real_standard_errors` followed by `imag_standard_errors`.

A full or leave-one-out denominator is usable only if its finite average sign
has absolute value strictly larger than
`minimum_absolute_average_sign`.  An unusable full denominator returns
`green=missing`; any unusable leave-one-out denominator returns missing
standard errors and covariance.  Individual block ratios are deliberately
not formed or returned, so zero or weak block sign sums are never divided by.
The selected raw Green values are not modified and the snapshots are not
mutated.
"""
function summarize_sign_reweighted_unequal_time_channel(
    snapshots::AbstractVector{<:UnequalTimeGreensSnapshot};
    spin::Integer,
    destination_orbital::Integer,
    source_orbital::Integer,
    kx_index::Integer,
    ky_index::Integer,
    block_size::Integer,
    minimum_absolute_average_sign::Real=0.0,
)
    measurement_count = length(snapshots)
    measurement_count > 0 || throw(ArgumentError(
        "snapshots must be nonempty",
    ))
    requested_block_size = try
        Int(block_size)
    catch error
        error isa InexactError || error isa OverflowError || rethrow()
        throw(ArgumentError("block_size must fit in an Int"))
    end
    requested_block_size > 0 || throw(ArgumentError(
        "block_size must be positive",
    ))
    minimum_sign = float(minimum_absolute_average_sign)
    isfinite(minimum_sign) && 0 <= minimum_sign < 1 || throw(ArgumentError(
        "minimum_absolute_average_sign must be finite and in [0, 1)",
    ))

    reference = first(snapshots)
    layout = _validate_unequal_time_snapshot_layout(reference, 1)
    spin_index = _checked_unequal_time_channel_index(spin, 2, "spin")
    destination_index = _checked_unequal_time_channel_index(
        destination_orbital,
        2,
        "destination_orbital",
    )
    source_index = _checked_unequal_time_channel_index(
        source_orbital,
        2,
        "source_orbital",
    )
    checked_kx_index = _checked_unequal_time_channel_index(
        kx_index,
        layout.L,
        "kx_index",
    )
    checked_ky_index = _checked_unequal_time_channel_index(
        ky_index,
        layout.L,
        "ky_index",
    )

    T = typeof(minimum_sign)
    for (measurement, snapshot) in enumerate(snapshots)
        snapshot_layout = _validate_unequal_time_snapshot_layout(
            snapshot,
            measurement,
        )
        snapshot_layout == layout || throw(DimensionMismatch(
            "all snapshots must have identical Green-function shapes",
        ))
        snapshot.tau == reference.tau || throw(ArgumentError(
            "all snapshots must have identical tau grids",
        ))
        snapshot.momenta == reference.momenta || throw(ArgumentError(
            "all snapshots must have identical momentum grids",
        ))
        T = promote_type(
            T,
            eltype(snapshot.tau),
            eltype(snapshot.momenta),
            typeof(float(snapshot.configuration_sign)),
            typeof(float(real(zero(eltype(snapshot.momentum_space))))),
        )
    end
    threshold = convert(T, minimum_sign)
    signs = Vector{T}(undef, measurement_count)
    raw_channel = Matrix{Complex{T}}(
        undef,
        measurement_count,
        layout.number_of_times,
    )
    for (measurement, snapshot) in enumerate(snapshots)
        sign_value = convert(T, snapshot.configuration_sign)
        isfinite(sign_value) &&
            (sign_value == one(T) || sign_value == -one(T)) ||
            throw(ArgumentError(
                "snapshots[$measurement].configuration_sign must be +1 or -1",
            ))
        signs[measurement] = sign_value
        for time_index in 1:layout.number_of_times
            raw_channel[measurement, time_index] = snapshot.momentum_space[
                spin_index,
                destination_index,
                source_index,
                checked_kx_index,
                checked_ky_index,
                time_index,
            ]
        end
    end

    block_ranges = _unequal_time_block_ranges(
        measurement_count,
        requested_block_size,
    )
    block_lengths = length.(block_ranges)
    block_count = length(block_ranges)
    block_sign_denominators = Vector{T}(undef, block_count)
    block_signed_numerators = zeros(
        Complex{T},
        block_count,
        layout.number_of_times,
    )
    for (block, block_range) in enumerate(block_ranges)
        block_sign_denominators[block] = sum(@view signs[block_range])
        for measurement in block_range, time_index in 1:layout.number_of_times
            block_signed_numerators[block, time_index] +=
                signs[measurement] * raw_channel[measurement, time_index]
        end
    end
    block_average_signs = block_sign_denominators ./ block_lengths
    full_sign_denominator = sum(block_sign_denominators)
    full_signed_numerator = vec(sum(block_signed_numerators; dims=1))
    full_average_sign = full_sign_denominator / measurement_count
    full_denominator_usable = isfinite(full_average_sign) &&
        abs(full_average_sign) > threshold
    channel_values_finite = all(isfinite, raw_channel)
    full_numerator_finite = all(isfinite, full_signed_numerator)

    green = if full_denominator_usable && channel_values_finite &&
               full_numerator_finite
        full_signed_numerator ./ full_sign_denominator
    else
        missing
    end
    central_estimate_valid = !ismissing(green) && all(isfinite, green)

    leave_one_out_sign_denominators = Vector{T}(undef, block_count)
    leave_one_out_average_signs =
        Vector{Union{Missing,T}}(undef, block_count)
    fill!(leave_one_out_average_signs, missing)
    leave_one_out_denominator_usable = falses(block_count)
    for block in 1:block_count
        remaining_count = measurement_count - block_lengths[block]
        remaining_sign_denominator =
            full_sign_denominator - block_sign_denominators[block]
        leave_one_out_sign_denominators[block] = remaining_sign_denominator
        if remaining_count > 0
            remaining_average_sign =
                remaining_sign_denominator / remaining_count
            leave_one_out_average_signs[block] = remaining_average_sign
            leave_one_out_denominator_usable[block] =
                isfinite(remaining_average_sign) &&
                abs(remaining_average_sign) > threshold
        end
    end
    leave_one_out_denominators_usable = block_count >= 2 &&
        all(leave_one_out_denominator_usable)
    leave_one_out_denominator_failure_count =
        count(!, leave_one_out_denominator_usable)
    finite_leave_one_out_average_signs = block_count >= 2 && all(
        value -> !ismissing(value) && isfinite(value),
        leave_one_out_average_signs,
    )
    minimum_absolute_leave_one_out_average_sign =
        finite_leave_one_out_average_signs ? minimum(
            abs(value) for value in leave_one_out_average_signs
        ) : missing

    covariance = missing
    real_standard_errors = missing
    imag_standard_errors = missing
    jackknife_valid = false
    if central_estimate_valid && leave_one_out_denominators_usable
        leave_one_out_vectors = Matrix{T}(
            undef,
            block_count,
            2 * layout.number_of_times,
        )
        for block in 1:block_count
            ratio = (
                full_signed_numerator .-
                view(block_signed_numerators, block, :)
            ) ./ leave_one_out_sign_denominators[block]
            leave_one_out_vectors[
                block,
                1:layout.number_of_times,
            ] .= real.(ratio)
            leave_one_out_vectors[
                block,
                (layout.number_of_times + 1):(2 * layout.number_of_times),
            ] .= imag.(ratio)
        end
        if all(isfinite, leave_one_out_vectors)
            jackknife_center = vec(sum(leave_one_out_vectors; dims=1)) ./
                               block_count
            candidate_covariance = zeros(
                T,
                2 * layout.number_of_times,
                2 * layout.number_of_times,
            )
            for block in 1:block_count
                deviation = view(leave_one_out_vectors, block, :) .-
                            jackknife_center
                candidate_covariance .+= deviation * transpose(deviation)
            end
            candidate_covariance .*= convert(T, block_count - 1) /
                                     convert(T, block_count)
            if all(isfinite, candidate_covariance)
                candidate_standard_errors = sqrt.(max.(
                    zero(T),
                    [candidate_covariance[index, index]
                     for index in axes(candidate_covariance, 1)],
                ))
                covariance = candidate_covariance
                real_standard_errors = candidate_standard_errors[
                    1:layout.number_of_times
                ]
                imag_standard_errors = candidate_standard_errors[
                    (layout.number_of_times + 1):(2 * layout.number_of_times)
                ]
                jackknife_valid = true
            end
        end
    end

    return (
        measurement_count=measurement_count,
        requested_block_size=requested_block_size,
        block_count=block_count,
        block_lengths=block_lengths,
        minimum_block_length=minimum(block_lengths),
        maximum_block_length=maximum(block_lengths),
        spin=spin_index,
        destination_orbital=destination_index,
        source_orbital=source_index,
        kx_index=checked_kx_index,
        ky_index=checked_ky_index,
        kx=convert(T, reference.momenta[checked_kx_index]),
        ky=convert(T, reference.momenta[checked_ky_index]),
        tau=T.(reference.tau),
        configured_minimum_absolute_average_sign=threshold,
        block_sign_denominators=block_sign_denominators,
        block_average_signs=block_average_signs,
        block_signed_numerators=block_signed_numerators,
        full_sign_denominator=full_sign_denominator,
        full_average_sign=full_average_sign,
        full_signed_numerator=full_signed_numerator,
        full_denominator_usable=full_denominator_usable,
        leave_one_out_sign_denominators=leave_one_out_sign_denominators,
        leave_one_out_average_signs=leave_one_out_average_signs,
        leave_one_out_denominator_usable=
            leave_one_out_denominator_usable,
        leave_one_out_denominators_usable=
            leave_one_out_denominators_usable,
        leave_one_out_denominator_failure_count=
            leave_one_out_denominator_failure_count,
        minimum_absolute_leave_one_out_average_sign=
            minimum_absolute_leave_one_out_average_sign,
        denominator_guard_passed=full_denominator_usable &&
            leave_one_out_denominators_usable,
        channel_values_finite=channel_values_finite,
        central_estimate_valid=central_estimate_valid,
        jackknife_valid=jackknife_valid,
        covariance_ready=central_estimate_valid && jackknife_valid,
        green=green,
        real_standard_errors=real_standard_errors,
        imag_standard_errors=imag_standard_errors,
        covariance=covariance,
    )
end

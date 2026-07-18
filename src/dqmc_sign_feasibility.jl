"""
    SignFeasibilityCriteria(; kwargs...)

Independent decision criteria for Milestone 6A sign diagnostics. The average-
sign gate is deliberately only `average_sign - 2*SE > average_sign_threshold`.
Numerical validity, sampling sufficiency, and filling resolution are reported
separately and are combined only in `correlation_measurement_ready`.
"""
struct SignFeasibilityCriteria{T<:AbstractFloat}
    average_sign_threshold::T
    maximum_average_sign_standard_error::T
    minimum_block_count::Int
    maximum_block_drift_zscore::T
    minimum_sign_signal_to_noise::T
    minimum_absolute_average_sign::T
    maximum_greens_recompute_error::T
    u0_filling_tolerance::T
end

function SignFeasibilityCriteria(;
    average_sign_threshold::Real=0.2,
    maximum_average_sign_standard_error::Real=0.02,
    minimum_block_count::Integer=8,
    maximum_block_drift_zscore::Real=2.0,
    minimum_sign_signal_to_noise::Real=5.0,
    minimum_absolute_average_sign::Real=1.0e-8,
    maximum_greens_recompute_error::Real=1.0e-8,
    u0_filling_tolerance::Real=1.0e-10,
)
    values = promote(
        float(average_sign_threshold),
        float(maximum_average_sign_standard_error),
        float(maximum_block_drift_zscore),
        float(minimum_sign_signal_to_noise),
        float(minimum_absolute_average_sign),
        float(maximum_greens_recompute_error),
        float(u0_filling_tolerance),
    )
    T = typeof(first(values))
    criteria = SignFeasibilityCriteria{T}(
        values[1],
        values[2],
        Int(minimum_block_count),
        values[3],
        values[4],
        values[5],
        values[6],
        values[7],
    )
    -one(T) <= criteria.average_sign_threshold < one(T) || throw(
        ArgumentError("average_sign_threshold must be finite and in [-1, 1)"),
    )
    isfinite(criteria.average_sign_threshold) || throw(
        ArgumentError("average_sign_threshold must be finite and in [-1, 1)"),
    )
    isfinite(criteria.maximum_average_sign_standard_error) &&
        criteria.maximum_average_sign_standard_error >= zero(T) || throw(
        ArgumentError(
            "maximum_average_sign_standard_error must be finite and nonnegative",
        ),
    )
    criteria.minimum_block_count > 0 || throw(
        ArgumentError("minimum_block_count must be positive"),
    )
    isfinite(criteria.maximum_block_drift_zscore) &&
        criteria.maximum_block_drift_zscore >= zero(T) || throw(
        ArgumentError("maximum_block_drift_zscore must be finite and nonnegative"),
    )
    isfinite(criteria.minimum_sign_signal_to_noise) &&
        criteria.minimum_sign_signal_to_noise >= zero(T) || throw(
        ArgumentError("minimum_sign_signal_to_noise must be finite and nonnegative"),
    )
    isfinite(criteria.minimum_absolute_average_sign) &&
        zero(T) <= criteria.minimum_absolute_average_sign < one(T) || throw(
        ArgumentError(
            "minimum_absolute_average_sign must be finite and in [0, 1)",
        ),
    )
    isfinite(criteria.maximum_greens_recompute_error) &&
        criteria.maximum_greens_recompute_error >= zero(T) || throw(
        ArgumentError(
            "maximum_greens_recompute_error must be finite and nonnegative",
        ),
    )
    isfinite(criteria.u0_filling_tolerance) &&
        criteria.u0_filling_tolerance >= zero(T) || throw(
        ArgumentError("u0_filling_tolerance must be finite and nonnegative"),
    )
    return criteria
end

const DEFAULT_SIGN_FEASIBILITY_CRITERIA = SignFeasibilityCriteria()

"""
    average_sign_gate_passes(average_sign, standard_error, threshold)

Return the Milestone 6A average-sign gate. For finite inputs and a nonnegative
standard error, the returned Boolean contains no condition other than the
strict inequality `average_sign - 2*standard_error > threshold`.
"""
function average_sign_gate_passes(
    average_sign::Real,
    standard_error::Real,
    threshold::Real,
)
    return average_sign - 2 * standard_error > threshold
end

function _standard_error(values::AbstractVector{T}) where {T<:AbstractFloat}
    number_of_values = length(values)
    number_of_values >= 2 || return missing
    center = sum(values) / number_of_values
    squared_deviations = sum((value - center)^2 for value in values)
    return sqrt(squared_deviations / (number_of_values * (number_of_values - 1)))
end

const _THERMODYNAMIC_REWEIGHTING_FIELDS = (
    :density,
    :doping,
    :kinetic_energy_per_cell,
    :shifted_interaction_energy_per_cell,
    :internal_energy_per_cell,
    :grand_energy_per_cell,
)

function _contiguous_block_ranges(measurement_count::Int, block_size::Int)
    return [
        first_index:min(first_index + block_size - 1, measurement_count)
        for first_index in 1:block_size:measurement_count
    ]
end

function _jackknife_ratio_standard_errors(
    block_sign_sums::AbstractVector{T},
    block_signed_observable_sums::AbstractMatrix{T},
    block_lengths::AbstractVector{Int},
    sign_sum::T,
    signed_observable_sums::AbstractVector{T},
    measurement_count::Int,
    minimum_absolute_average_sign::T,
    full_denominator_usable::Bool,
) where {T<:AbstractFloat}
    number_of_blocks = length(block_sign_sums)
    number_of_observables = length(signed_observable_sums)
    size(block_signed_observable_sums) ==
        (number_of_blocks, number_of_observables) || throw(DimensionMismatch(
            "block signed-observable sums have incompatible dimensions",
        ))
    length(block_lengths) == number_of_blocks || throw(DimensionMismatch(
        "block sign sums and block lengths must have equal length",
    ))

    standard_errors = Vector{Union{Missing,T}}(undef, number_of_observables)
    fill!(standard_errors, missing)
    valid = falses(number_of_observables)
    leave_one_out_sign_sums = Vector{T}(undef, number_of_blocks)
    leave_one_out_average_signs =
        Vector{Union{Missing,T}}(undef, number_of_blocks)
    fill!(leave_one_out_average_signs, missing)
    denominator_usable = falses(number_of_blocks)

    for block in 1:number_of_blocks
        remaining_count = measurement_count - block_lengths[block]
        remaining_sign_sum = sign_sum - block_sign_sums[block]
        leave_one_out_sign_sums[block] = remaining_sign_sum
        if remaining_count > 0
            remaining_average_sign = remaining_sign_sum / remaining_count
            leave_one_out_average_signs[block] = remaining_average_sign
            denominator_usable[block] = isfinite(remaining_average_sign) &&
                abs(remaining_average_sign) > minimum_absolute_average_sign
        end
    end

    all_leave_one_out_denominators_usable = number_of_blocks >= 2 &&
        all(denominator_usable)
    denominator_failure_count = count(!, denominator_usable)
    finite_leave_one_out_average_signs = number_of_blocks >= 2 && all(
        value -> !ismissing(value) && isfinite(value),
        leave_one_out_average_signs,
    )
    minimum_absolute_leave_one_out_average_sign =
        finite_leave_one_out_average_signs ? minimum(
            abs(value) for value in leave_one_out_average_signs
        ) : missing

    if full_denominator_usable && all_leave_one_out_denominators_usable
        estimates = Matrix{T}(undef, number_of_blocks, number_of_observables)
        for block in 1:number_of_blocks, observable in 1:number_of_observables
            estimates[block, observable] = (
                signed_observable_sums[observable] -
                block_signed_observable_sums[block, observable]
            ) / leave_one_out_sign_sums[block]
        end
        for observable in 1:number_of_observables
            observable_estimates = @view estimates[:, observable]
            all(isfinite, observable_estimates) || continue
            center = sum(observable_estimates) / number_of_blocks
            squared_deviations = sum(
                (estimate - center)^2 for estimate in observable_estimates
            )
            standard_error = sqrt(
                convert(T, number_of_blocks - 1) /
                convert(T, number_of_blocks) * squared_deviations,
            )
            isfinite(standard_error) || continue
            standard_errors[observable] = standard_error
            valid[observable] = true
        end
    end

    return (
        standard_errors=standard_errors,
        valid=valid,
        leave_one_out_sign_sums=leave_one_out_sign_sums,
        leave_one_out_average_signs=leave_one_out_average_signs,
        leave_one_out_denominators_usable=
            all_leave_one_out_denominators_usable,
        leave_one_out_denominator_failure_count=denominator_failure_count,
        minimum_absolute_leave_one_out_average_sign=
            minimum_absolute_leave_one_out_average_sign,
    )
end

function _sign_reweighted_ratio_statistics(
    signs::AbstractVector{T},
    raw_observables::AbstractMatrix{T},
    block_ranges::AbstractVector{<:UnitRange{Int}},
    minimum_absolute_average_sign::T,
) where {T<:AbstractFloat}
    measurement_count = length(signs)
    size(raw_observables, 1) == measurement_count || throw(DimensionMismatch(
        "configuration signs and raw observables must have equal length",
    ))
    number_of_observables = size(raw_observables, 2)
    number_of_blocks = length(block_ranges)
    block_lengths = length.(block_ranges)
    block_sign_sums = Vector{T}(undef, number_of_blocks)
    block_signed_observable_sums =
        zeros(T, number_of_blocks, number_of_observables)
    raw_observable_sums = zeros(T, number_of_observables)
    signed_observable_sums = zeros(T, number_of_observables)

    for measurement in 1:measurement_count
        for observable in 1:number_of_observables
            raw_value = raw_observables[measurement, observable]
            raw_observable_sums[observable] += raw_value
        end
    end
    for (block, range) in enumerate(block_ranges)
        block_sign_sums[block] = sum(@view signs[range])
        for measurement in range, observable in 1:number_of_observables
            block_signed_observable_sums[block, observable] +=
                signs[measurement] * raw_observables[measurement, observable]
        end
    end

    sign_sum = sum(block_sign_sums)
    for observable in 1:number_of_observables
        signed_observable_sums[observable] = sum(
            @view block_signed_observable_sums[:, observable]
        )
    end

    average_sign = sign_sum / measurement_count
    sign_denominator_usable = isfinite(average_sign) &&
        abs(average_sign) > minimum_absolute_average_sign
    estimates = Vector{Union{Missing,T}}(undef, number_of_observables)
    fill!(estimates, missing)
    if sign_denominator_usable
        for observable in 1:number_of_observables
            if isfinite(signed_observable_sums[observable])
                estimates[observable] =
                    signed_observable_sums[observable] / sign_sum
            end
        end
    end
    jackknife = _jackknife_ratio_standard_errors(
        block_sign_sums,
        block_signed_observable_sums,
        block_lengths,
        sign_sum,
        signed_observable_sums,
        measurement_count,
        minimum_absolute_average_sign,
        sign_denominator_usable,
    )
    return (
        block_lengths=block_lengths,
        block_sign_sums=block_sign_sums,
        block_signed_observable_sums=block_signed_observable_sums,
        raw_observable_sums=raw_observable_sums,
        signed_observable_sums=signed_observable_sums,
        sign_sum=sign_sum,
        average_sign=average_sign,
        sign_denominator_usable=sign_denominator_usable,
        estimates=estimates,
        jackknife=jackknife,
    )
end

function _block_drift_diagnostic(
    block_means::AbstractVector{T},
    maximum_zscore::T,
) where {T<:AbstractFloat}
    number_of_blocks = length(block_means)
    number_of_blocks >= 4 || return (
        evaluable=false,
        zscore=missing,
        stable=false,
    )
    split = fld(number_of_blocks, 2)
    early = block_means[1:split]
    late = block_means[(split + 1):end]
    early_mean = sum(early) / length(early)
    late_mean = sum(late) / length(late)
    early_se = _standard_error(early)
    late_se = _standard_error(late)
    combined_se = hypot(early_se, late_se)
    difference = abs(early_mean - late_mean)
    zscore = if iszero(combined_se)
        iszero(difference) ? zero(T) : T(Inf)
    else
        difference / combined_se
    end
    stable = !isnan(zscore) && zscore <= maximum_zscore
    return (evaluable=true, zscore=zscore, stable=stable)
end

function _jackknife_ratio_standard_error(
    block_sign_sums::AbstractVector{T},
    block_signed_filling_sums::AbstractVector{T},
    block_lengths::AbstractVector{Int},
    sign_sum::T,
    signed_filling_sum::T,
    measurement_count::Int,
    minimum_absolute_average_sign::T,
) where {T<:AbstractFloat}
    full_average_sign = sign_sum / measurement_count
    full_denominator_usable = isfinite(full_average_sign) &&
        abs(full_average_sign) > minimum_absolute_average_sign
    jackknife = _jackknife_ratio_standard_errors(
        block_sign_sums,
        reshape(block_signed_filling_sums, :, 1),
        block_lengths,
        sign_sum,
        T[signed_filling_sum],
        measurement_count,
        minimum_absolute_average_sign,
        full_denominator_usable,
    )
    return (
        valid=jackknife.valid[1],
        standard_error=jackknife.standard_errors[1],
    )
end

"""
    summarize_sign_reweighted_thermodynamics(
        configuration_signs, raw_observables;
        block_size, minimum_absolute_average_sign=0.0,
    )

Sign reweight the six fields returned by `measure_thermodynamic_observables`.
Every raw configuration estimator is multiplied by its configuration sign
exactly once. The central estimate for observable `a` is the ratio of total
sums `sum(s .* O_a) / sum(s)`, never an average of divided block ratios.

All observables use the same contiguous blocks. Their standard errors use the
leave-one-block-out ratios
`(X_a - X_b,a) / (S - S_b)`. The full denominator and every leave-one-out
denominator must have a finite average sign whose absolute value is strictly
larger than `minimum_absolute_average_sign`; otherwise the affected estimates
or errors are returned as `missing` and the corresponding validity flags are
false. A final short block is retained.
"""
function summarize_sign_reweighted_thermodynamics(
    configuration_signs::AbstractVector{<:Real},
    raw_observables::AbstractVector;
    block_size::Integer,
    minimum_absolute_average_sign::Real=0.0,
)
    measurement_count = length(configuration_signs)
    measurement_count > 0 || throw(ArgumentError(
        "configuration_signs and raw_observables must be nonempty",
    ))
    length(raw_observables) == measurement_count || throw(DimensionMismatch(
        "configuration_signs and raw_observables must have equal length",
    ))
    block_size > 0 || throw(ArgumentError("block_size must be positive"))
    minimum_sign = float(minimum_absolute_average_sign)
    isfinite(minimum_sign) && 0 <= minimum_sign < 1 || throw(ArgumentError(
        "minimum_absolute_average_sign must be finite and in [0, 1)",
    ))

    T = typeof(minimum_sign)
    for configuration_sign in configuration_signs
        T = promote_type(T, typeof(float(configuration_sign)))
    end
    for (measurement, observation) in enumerate(raw_observables)
        for name in _THERMODYNAMIC_REWEIGHTING_FIELDS
            hasproperty(observation, name) || throw(ArgumentError(
                "raw_observables[$measurement] has no field $name",
            ))
            value = getproperty(observation, name)
            value isa Real || throw(ArgumentError(
                "raw_observables[$measurement].$name must be real",
            ))
            T = promote_type(T, typeof(float(value)))
        end
    end

    signs = T.(configuration_signs)
    number_of_observables = length(_THERMODYNAMIC_REWEIGHTING_FIELDS)
    raw_values = Matrix{T}(undef, measurement_count, number_of_observables)
    for measurement in 1:measurement_count,
        (observable, name) in enumerate(_THERMODYNAMIC_REWEIGHTING_FIELDS)
        raw_values[measurement, observable] =
            getproperty(raw_observables[measurement], name)
    end
    block_ranges = _contiguous_block_ranges(measurement_count, Int(block_size))
    statistics = _sign_reweighted_ratio_statistics(
        signs,
        raw_values,
        block_ranges,
        convert(T, minimum_sign),
    )
    block_sign_means = statistics.block_sign_sums ./ statistics.block_lengths
    average_sign_standard_error = _standard_error(block_sign_means)
    jackknife = statistics.jackknife

    configuration_signs_finite = all(isfinite, signs)
    configuration_signs_valid = configuration_signs_finite && all(
        sign_value -> sign_value == one(T) || sign_value == -one(T),
        signs,
    )
    raw_observables_finite = all(isfinite, raw_values)
    finite_results = all(isfinite, (
        statistics.sign_sum,
        statistics.average_sign,
        statistics.block_sign_sums...,
        statistics.raw_observable_sums...,
        statistics.signed_observable_sums...,
    )) &&
        (ismissing(average_sign_standard_error) ||
         isfinite(average_sign_standard_error)) &&
        (ismissing(jackknife.minimum_absolute_leave_one_out_average_sign) ||
         isfinite(jackknife.minimum_absolute_leave_one_out_average_sign)) &&
        all(
            value -> ismissing(value) || isfinite(value),
            statistics.estimates,
        ) && all(
            value -> ismissing(value) || isfinite(value),
            jackknife.standard_errors,
        )

    return (
        measurement_count=measurement_count,
        requested_block_size=Int(block_size),
        block_count=length(block_ranges),
        minimum_block_length=minimum(statistics.block_lengths),
        maximum_block_length=maximum(statistics.block_lengths),
        configuration_sign_sum=statistics.sign_sum,
        average_sign=statistics.average_sign,
        average_sign_standard_error=average_sign_standard_error,
        configured_minimum_absolute_average_sign=convert(T, minimum_sign),
        sign_denominator_usable=statistics.sign_denominator_usable,
        leave_one_out_denominators_usable=
            jackknife.leave_one_out_denominators_usable,
        leave_one_out_denominator_failure_count=
            jackknife.leave_one_out_denominator_failure_count,
        minimum_absolute_leave_one_out_average_sign=
            jackknife.minimum_absolute_leave_one_out_average_sign,
        denominator_guard_passed=statistics.sign_denominator_usable &&
            jackknife.leave_one_out_denominators_usable,
        density=statistics.estimates[1],
        density_standard_error=jackknife.standard_errors[1],
        density_jackknife_valid=jackknife.valid[1],
        doping=statistics.estimates[2],
        doping_standard_error=jackknife.standard_errors[2],
        doping_jackknife_valid=jackknife.valid[2],
        kinetic_energy_per_cell=statistics.estimates[3],
        kinetic_energy_per_cell_standard_error=jackknife.standard_errors[3],
        kinetic_energy_per_cell_jackknife_valid=jackknife.valid[3],
        shifted_interaction_energy_per_cell=statistics.estimates[4],
        shifted_interaction_energy_per_cell_standard_error=
            jackknife.standard_errors[4],
        shifted_interaction_energy_per_cell_jackknife_valid=
            jackknife.valid[4],
        internal_energy_per_cell=statistics.estimates[5],
        internal_energy_per_cell_standard_error=jackknife.standard_errors[5],
        internal_energy_per_cell_jackknife_valid=jackknife.valid[5],
        grand_energy_per_cell=statistics.estimates[6],
        grand_energy_per_cell_standard_error=jackknife.standard_errors[6],
        grand_energy_per_cell_jackknife_valid=jackknife.valid[6],
        configuration_signs_finite=configuration_signs_finite,
        configuration_signs_valid=configuration_signs_valid,
        raw_observables_finite=raw_observables_finite,
        finite_results=finite_results,
    )
end

"""
    summarize_sign_feasibility(configuration_signs, raw_total_densities;
                               block_size, criteria, ...)

Analyze one chain without constructing or changing a DQMC state.
`raw_total_densities[i]` is the unsigned configuration estimator returned by
`measure_basic_observables`; consequently `signed_filling_sum` is formed
exactly once as `sum(configuration_signs .* raw_total_densities)`.

Contiguous block means estimate the average-sign standard error. Block
stability compares the early and late halves of those block means using the
configured drift z-score. A final short block is retained and its length is
reported. The filling is omitted when the configured sign-denominator guard
fails rather than returning a misleading ratio.
"""
function summarize_sign_feasibility(
    configuration_signs::AbstractVector{<:Real},
    raw_total_densities::AbstractVector{<:Real};
    block_size::Integer,
    criteria::SignFeasibilityCriteria=DEFAULT_SIGN_FEASIBILITY_CRITERIA,
    maximum_greens_recompute_error::Real=0.0,
    acceptance_rate::Union{Missing,Real}=missing,
)
    measurement_count = length(configuration_signs)
    measurement_count > 0 || throw(ArgumentError(
        "configuration_signs and raw_total_densities must be nonempty",
    ))
    length(raw_total_densities) == measurement_count || throw(DimensionMismatch(
        "configuration_signs and raw_total_densities must have equal length",
    ))
    block_size > 0 || throw(ArgumentError("block_size must be positive"))

    T = promote_type(
        typeof(float(first(configuration_signs))),
        typeof(float(first(raw_total_densities))),
        typeof(float(maximum_greens_recompute_error)),
        typeof(criteria.average_sign_threshold),
    )
    signs = T.(configuration_signs)
    raw_fillings = T.(raw_total_densities)
    correction = convert(T, maximum_greens_recompute_error)
    acceptance = ismissing(acceptance_rate) ? missing : convert(T, acceptance_rate)

    block_ranges = _contiguous_block_ranges(measurement_count, Int(block_size))
    ratio_statistics = _sign_reweighted_ratio_statistics(
        signs,
        reshape(raw_fillings, :, 1),
        block_ranges,
        convert(T, criteria.minimum_absolute_average_sign),
    )
    block_lengths = ratio_statistics.block_lengths
    block_sign_sums = ratio_statistics.block_sign_sums
    block_signed_filling_sums = vec(
        ratio_statistics.block_signed_observable_sums,
    )
    block_sign_means = block_sign_sums ./ block_lengths

    sign_sum = ratio_statistics.sign_sum
    signed_filling_sum = only(ratio_statistics.signed_observable_sums)
    raw_filling_sum = only(ratio_statistics.raw_observable_sums)
    average_sign = ratio_statistics.average_sign
    average_sign_standard_error = _standard_error(block_sign_means)
    average_sign_lower_2se = ismissing(average_sign_standard_error) ? missing :
        average_sign - 2 * average_sign_standard_error
    average_sign_gate_passed = !ismissing(average_sign_standard_error) &&
        average_sign_gate_passes(
            average_sign,
            average_sign_standard_error,
            criteria.average_sign_threshold,
        )

    block_drift = _block_drift_diagnostic(
        block_sign_means,
        convert(T, criteria.maximum_block_drift_zscore),
    )
    block_count = length(block_ranges)
    block_count_sufficient = block_count >= criteria.minimum_block_count
    average_sign_precision_sufficient =
        !ismissing(average_sign_standard_error) &&
        isfinite(average_sign_standard_error) &&
        average_sign_standard_error <=
            criteria.maximum_average_sign_standard_error
    sampling_sufficient = block_count_sufficient && block_drift.stable &&
        average_sign_precision_sufficient

    sign_denominator_usable = ratio_statistics.sign_denominator_usable
    total_density = only(ratio_statistics.estimates)
    sign_signal_to_noise = if ismissing(average_sign_standard_error) ||
                              !isfinite(average_sign_standard_error) ||
                              !isfinite(average_sign)
        missing
    elseif iszero(average_sign_standard_error)
        iszero(average_sign) ? zero(T) : T(Inf)
    else
        abs(average_sign) / average_sign_standard_error
    end
    sign_signal_to_noise_sufficient =
        !ismissing(sign_signal_to_noise) &&
        sign_signal_to_noise >= criteria.minimum_sign_signal_to_noise
    filling_resolution_sufficient = sign_denominator_usable &&
        sign_signal_to_noise_sufficient
    filling_jackknife = (
        valid=ratio_statistics.jackknife.valid[1],
        standard_error=ratio_statistics.jackknife.standard_errors[1],
    )

    configuration_signs_finite = all(isfinite, signs)
    configuration_signs_valid = configuration_signs_finite && all(
        sign_value -> sign_value == one(T) || sign_value == -one(T),
        signs,
    )
    raw_fillings_finite = all(isfinite, raw_fillings)
    acceptance_valid = ismissing(acceptance) ||
        (isfinite(acceptance) && zero(T) <= acceptance <= one(T))
    correction_finite_nonnegative = isfinite(correction) && correction >= zero(T)
    correction_guard_passed = correction_finite_nonnegative &&
        correction <= criteria.maximum_greens_recompute_error
    finite_results = all(isfinite, (
        sign_sum,
        signed_filling_sum,
        raw_filling_sum,
        average_sign,
    )) &&
        (ismissing(average_sign_standard_error) ||
         isfinite(average_sign_standard_error)) &&
        (ismissing(average_sign_lower_2se) || isfinite(average_sign_lower_2se)) &&
        (ismissing(total_density) || isfinite(total_density)) &&
        (ismissing(filling_jackknife.standard_error) ||
         isfinite(filling_jackknife.standard_error)) &&
        (ismissing(block_drift.zscore) ||
         isfinite(block_drift.zscore) || block_drift.zscore == T(Inf))
    numerical_valid = configuration_signs_valid && raw_fillings_finite &&
        acceptance_valid && correction_guard_passed && finite_results
    correlation_measurement_ready = average_sign_gate_passed &&
        numerical_valid && sampling_sufficient && filling_resolution_sufficient

    return (
        measurement_count=measurement_count,
        requested_block_size=Int(block_size),
        block_count=block_count,
        minimum_block_length=minimum(block_lengths),
        maximum_block_length=maximum(block_lengths),
        configuration_sign_sum=sign_sum,
        raw_filling_sum=raw_filling_sum,
        signed_filling_sum=signed_filling_sum,
        average_sign=average_sign,
        average_sign_standard_error=average_sign_standard_error,
        average_sign_lower_2se=average_sign_lower_2se,
        configured_average_sign_threshold=criteria.average_sign_threshold,
        average_sign_gate_passed=average_sign_gate_passed,
        block_count_sufficient=block_count_sufficient,
        block_stability_evaluable=block_drift.evaluable,
        block_drift_zscore=block_drift.zscore,
        block_stable=block_drift.stable,
        average_sign_precision_sufficient=average_sign_precision_sufficient,
        sampling_sufficient=sampling_sufficient,
        sign_denominator_usable=sign_denominator_usable,
        sign_signal_to_noise=sign_signal_to_noise,
        sign_signal_to_noise_sufficient=sign_signal_to_noise_sufficient,
        filling_resolution_sufficient=filling_resolution_sufficient,
        total_density=total_density,
        total_density_standard_error=filling_jackknife.standard_error,
        filling_jackknife_valid=filling_jackknife.valid,
        configuration_signs_finite=configuration_signs_finite,
        configuration_signs_valid=configuration_signs_valid,
        raw_fillings_finite=raw_fillings_finite,
        acceptance_rate=acceptance,
        acceptance_valid=acceptance_valid,
        maximum_greens_recompute_error=correction,
        correction_guard_passed=correction_guard_passed,
        finite_results=finite_results,
        numerical_valid=numerical_valid,
        correlation_measurement_ready=correlation_measurement_ready,
    )
end

"""
    run_dqmc_sign_feasibility(parameters; case_name, phase, block_size, ...)

Run one explicitly requested Milestone 6A chain using the unchanged Milestone
5 update kernel. Per-configuration total density is raw when measured and is
multiplied by the configuration sign exactly once by
`summarize_sign_feasibility`.
"""
function run_dqmc_sign_feasibility(
    parameters::DQMCParameters{T};
    case_name::AbstractString="sign_feasibility",
    phase::Union{Symbol,AbstractString}=:fixed_mu_diagnostics,
    block_size::Integer,
    criteria::SignFeasibilityCriteria=DEFAULT_SIGN_FEASIBILITY_CRITERIA,
    model_parameters::ModelParameters=DEFAULT_PARAMETERS,
) where {T}
    block_size > 0 || throw(ArgumentError("block_size must be positive"))
    dimensions = validate_dqmc_parameters(parameters)
    state = initialize_dqmc_state(
        parameters;
        model_parameters=model_parameters,
    )
    start_time = time_ns()
    for _ in 1:parameters.burnin_sweeps
        sweep_all_time_slices!(state)
    end

    configuration_signs = Vector{T}(undef, parameters.measurement_sweeps)
    raw_total_densities = Vector{T}(undef, parameters.measurement_sweeps)
    for measurement in 1:parameters.measurement_sweeps
        for _ in 1:parameters.sweeps_between_measurements
            sweep_all_time_slices!(state)
        end
        observables = measure_basic_observables(state)
        configuration_signs[measurement] = observables.configuration_sign
        raw_total_densities[measurement] = observables.total_density
    end
    recompute_equal_time_greens!(state)
    runtime_seconds = (time_ns() - start_time) / 1.0e9
    acceptance_rate = state.accepted_updates / state.attempted_updates
    diagnostics = summarize_sign_feasibility(
        configuration_signs,
        raw_total_densities;
        block_size=block_size,
        criteria=criteria,
        maximum_greens_recompute_error=state.maximum_greens_recompute_error,
        acceptance_rate=acceptance_rate,
    )

    u0_control_applicable = iszero(parameters.U)
    u0_expected_total_density = if u0_control_applicable
        free_density = free_density_matrix(state.kinetic_mu, parameters.beta)
        total_filling_per_physical_site(free_density, parameters.L)
    else
        missing
    end
    u0_filling_error = if u0_control_applicable &&
                          !ismissing(diagnostics.total_density)
        abs(diagnostics.total_density - u0_expected_total_density)
    else
        missing
    end
    u0_sign_control_passed = u0_control_applicable ?
        all(==(one(T)), configuration_signs) : missing
    u0_filling_control_passed = u0_control_applicable ?
        (!ismissing(u0_filling_error) &&
         u0_filling_error <= criteria.u0_filling_tolerance) : missing
    u0_control_passed = u0_control_applicable ?
        (u0_sign_control_passed && u0_filling_control_passed) : missing

    return merge(
        (
            phase=String(phase),
            case=String(case_name),
            L=parameters.L,
            beta=parameters.beta,
            dtau=parameters.dtau,
            Ltau=dimensions.Ltau,
            U=parameters.U,
            mu=parameters.mu,
            matrix_dimension=dimensions.Nsp,
            burnin_sweeps=parameters.burnin_sweeps,
            measurement_sweeps=parameters.measurement_sweeps,
            sweeps_between_measurements=parameters.sweeps_between_measurements,
            seed=parameters.seed,
            runtime_seconds=runtime_seconds,
        ),
        diagnostics,
        (
            u0_control_applicable=u0_control_applicable,
            u0_expected_total_density=u0_expected_total_density,
            u0_filling_error=u0_filling_error,
            u0_sign_control_passed=u0_sign_control_passed,
            u0_filling_control_passed=u0_filling_control_passed,
            u0_control_passed=u0_control_passed,
        ),
    )
end

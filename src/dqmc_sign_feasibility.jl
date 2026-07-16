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
    average_sign_threshold::Real=0.1,
    maximum_average_sign_standard_error::Real=0.05,
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
    number_of_blocks = length(block_sign_sums)
    number_of_blocks >= 2 || return (valid=false, standard_error=missing)
    estimates = Vector{T}(undef, number_of_blocks)
    for block in 1:number_of_blocks
        remaining_count = measurement_count - block_lengths[block]
        remaining_count > 0 || return (valid=false, standard_error=missing)
        remaining_sign_sum = sign_sum - block_sign_sums[block]
        remaining_average_sign = remaining_sign_sum / remaining_count
        isfinite(remaining_average_sign) &&
            abs(remaining_average_sign) > minimum_absolute_average_sign ||
            return (valid=false, standard_error=missing)
        estimates[block] = (
            signed_filling_sum - block_signed_filling_sums[block]
        ) / remaining_sign_sum
    end
    all(isfinite, estimates) || return (valid=false, standard_error=missing)
    center = sum(estimates) / number_of_blocks
    squared_deviations = sum((estimate - center)^2 for estimate in estimates)
    standard_error = sqrt(
        (number_of_blocks - 1) / number_of_blocks * squared_deviations,
    )
    return (valid=isfinite(standard_error), standard_error=standard_error)
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

    block_ranges = collect(1:Int(block_size):measurement_count)
    block_ranges = [
        first_index:min(first_index + Int(block_size) - 1, measurement_count)
        for first_index in block_ranges
    ]
    block_lengths = length.(block_ranges)
    block_sign_sums = T[sum(@view signs[range]) for range in block_ranges]
    block_signed_filling_sums = T[
        sum(signs[index] * raw_fillings[index] for index in range)
        for range in block_ranges
    ]
    block_sign_means = block_sign_sums ./ block_lengths

    sign_sum = sum(signs)
    signed_filling_sum = sum(
        signs[index] * raw_fillings[index] for index in eachindex(signs)
    )
    raw_filling_sum = sum(raw_fillings)
    average_sign = sign_sum / measurement_count
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

    sign_denominator_usable = isfinite(average_sign) &&
        abs(average_sign) > criteria.minimum_absolute_average_sign
    total_density = sign_denominator_usable && isfinite(signed_filling_sum) ?
        signed_filling_sum / sign_sum : missing
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
    filling_jackknife = sign_denominator_usable ?
        _jackknife_ratio_standard_error(
            block_sign_sums,
            block_signed_filling_sums,
            block_lengths,
            sign_sum,
            signed_filling_sum,
            measurement_count,
            convert(T, criteria.minimum_absolute_average_sign),
        ) : (valid=false, standard_error=missing)

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

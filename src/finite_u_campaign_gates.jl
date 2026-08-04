const FIXED_DENSITY_TARGET_FILLING = 1.46021875

const FIXED_DENSITY_CALIBRATION_WORKLOAD = (
    burnin_sweeps=200,
    measurement_sweeps=800,
    sweeps_between_measurements=2,
    block_size=100,
    stabilization_period=5,
)

const FIXED_DENSITY_PRODUCTION_WORKLOAD = (
    burnin_sweeps=400,
    measurement_sweeps=1600,
    sweeps_between_measurements=2,
    block_size=100,
    stabilization_period=5,
)

const RESOURCE_PROJECTED_TWO_CHAIN_SECONDS_LIMIT = 86400.0
const RESOURCE_PEAK_RSS_BYTES_LIMIT = 5_905_580_032
const RESOURCE_GREENS_RECOMPUTE_ERROR_LIMIT = 1.0e-8

_finite_real(value) =
    value isa Real && !(value isa Bool) && isfinite(value)
_positive_finite_real(value) = _finite_real(value) && value > 0
_nonnegative_finite_real(value) = _finite_real(value) && value >= 0

_indexed_reason(noun, index, reason) = Symbol(
    string(noun, "_", index, "_", reason),
)

function _ordered_unique_reasons(reasons)
    return Tuple(unique(reasons))
end

function _calibration_decision(
    status,
    advance,
    selected_mu,
    unrounded_mu,
    compatible_point_mu,
    bracket_mus,
    sorted_mus,
    reasons,
)
    return (
        status=status,
        advance=advance,
        selected_mu=selected_mu,
        interpolated_mu_unrounded=unrounded_mu,
        compatible_point_mu=compatible_point_mu,
        bracket_mus=bracket_mus,
        sorted_mus=sorted_mus,
        reasons=_ordered_unique_reasons(reasons),
        automatic_extension_authorized=false,
    )
end

function _validate_expected_calibration_workload(expected_workload)
    for field in propertynames(FIXED_DENSITY_CALIBRATION_WORKLOAD)
        hasproperty(expected_workload, field) || throw(
            ArgumentError(
                string("expected workload is missing field ", field),
            ),
        )
        value = getproperty(expected_workload, field)
        value isa Integer && !(value isa Bool) && value >= 0 || throw(
            ArgumentError(
                string("expected workload field ", field, " is invalid"),
            ),
        )
    end
    expected_workload.measurement_sweeps > 0 ||
        throw(ArgumentError("measurement sweeps must be positive"))
    expected_workload.sweeps_between_measurements > 0 ||
        throw(ArgumentError("sweeps between measurements must be positive"))
    expected_workload.block_size > 0 ||
        throw(ArgumentError("block size must be positive"))
    expected_workload.stabilization_period > 0 ||
        throw(ArgumentError("stabilization period must be positive"))
    return nothing
end

"""
    assess_three_point_fixed_density_calibration(points;
        target_filling=FIXED_DENSITY_TARGET_FILLING,
        expected_workload=FIXED_DENSITY_CALIBRATION_WORKLOAD)

Make a deterministic fixed-density calibration decision from exactly three
completed result records. Each record must provide mu, total_density,
filling_jackknife_standard_error, filling_jackknife_valid, completed,
the five workload fields, measurement_count,
full_sign_denominator_usable, and
leave_one_out_denominators_usable.

Valid records are sorted by increasing chemical potential. Densities must be
monotone nondecreasing. A point whose density differs from the target by at
most two standard errors is accepted directly. Otherwise the target must lie
between adjacent unrounded density means; their linear interpolation is
rounded exactly once to four decimal places. Invalid or unbracketed inputs
stop and never authorize automatic bracket extension.
"""
function assess_three_point_fixed_density_calibration(
    points;
    target_filling=FIXED_DENSITY_TARGET_FILLING,
    expected_workload=FIXED_DENSITY_CALIBRATION_WORKLOAD,
)
    _validate_expected_calibration_workload(expected_workload)
    records = collect(points)
    reasons = Symbol[]
    if length(records) != 3
        push!(reasons, :requires_exactly_three_points)
        return _calibration_decision(
            :STOP_INVALID,
            false,
            nothing,
            nothing,
            nothing,
            nothing,
            (),
            reasons,
        )
    end
    if !_finite_real(target_filling)
        push!(reasons, :target_filling_not_finite)
    end

    required_fields = (
        :mu,
        :total_density,
        :filling_jackknife_standard_error,
        :filling_jackknife_valid,
        :completed,
        :burnin_sweeps,
        :measurement_sweeps,
        :sweeps_between_measurements,
        :block_size,
        :stabilization_period,
        :measurement_count,
        :full_sign_denominator_usable,
        :leave_one_out_denominators_usable,
    )
    workload_fields = propertynames(FIXED_DENSITY_CALIBRATION_WORKLOAD)

    for (index, record) in pairs(records)
        missing_fields = Symbol[]
        for field in required_fields
            if !hasproperty(record, field)
                push!(missing_fields, field)
                push!(
                    reasons,
                    _indexed_reason("point", index, string("missing_", field)),
                )
            end
        end
        isempty(missing_fields) || continue

        _finite_real(record.mu) || push!(
            reasons,
            _indexed_reason("point", index, "mu_not_finite"),
        )
        _finite_real(record.total_density) || push!(
            reasons,
            _indexed_reason("point", index, "density_not_finite"),
        )
        _positive_finite_real(
            record.filling_jackknife_standard_error,
        ) || push!(
            reasons,
            _indexed_reason(
                "point",
                index,
                "density_standard_error_not_positive_finite",
            ),
        )
        record.filling_jackknife_valid === true || push!(
            reasons,
            _indexed_reason("point", index, "density_jackknife_not_valid"),
        )
        record.completed === true || push!(
            reasons,
            _indexed_reason("point", index, "not_completed"),
        )
        record.full_sign_denominator_usable === true || push!(
            reasons,
            _indexed_reason("point", index, "full_denominator_not_usable"),
        )
        record.leave_one_out_denominators_usable === true || push!(
            reasons,
            _indexed_reason(
                "point",
                index,
                "leave_one_out_denominators_not_usable",
            ),
        )
        for field in workload_fields
            getproperty(record, field) ==
                getproperty(expected_workload, field) || push!(
                    reasons,
                    _indexed_reason(
                        "point",
                        index,
                        string(field, "_mismatch"),
                    ),
                )
        end
        record.measurement_count ==
            expected_workload.measurement_sweeps || push!(
                reasons,
                _indexed_reason(
                    "point",
                    index,
                    "measurement_count_mismatch",
                ),
            )
    end

    if !isempty(reasons)
        return _calibration_decision(
            :STOP_INVALID,
            false,
            nothing,
            nothing,
            nothing,
            nothing,
            (),
            reasons,
        )
    end

    sorted_records = sort(records; by=record -> record.mu)
    sorted_mus = Tuple(record.mu for record in sorted_records)
    if any(
        sorted_records[index].mu >= sorted_records[index + 1].mu
        for index in 1:2
    )
        push!(reasons, :mu_values_not_strictly_increasing)
    end
    if any(
        sorted_records[index].total_density >
            sorted_records[index + 1].total_density
        for index in 1:2
    )
        push!(reasons, :density_not_monotone_nondecreasing)
    end
    if !isempty(reasons)
        return _calibration_decision(
            :STOP_INVALID,
            false,
            nothing,
            nothing,
            nothing,
            nothing,
            sorted_mus,
            reasons,
        )
    end

    compatible_indices = [
        index for index in eachindex(sorted_records)
        if abs(sorted_records[index].total_density - target_filling) <=
           2 * sorted_records[index].filling_jackknife_standard_error
    ]
    if !isempty(compatible_indices)
        best_index = first(sort(
            compatible_indices;
            by=index -> (
                abs(
                    sorted_records[index].total_density -
                    target_filling,
                ),
                sorted_records[index].mu,
            ),
        ))
        selected_mu = sorted_records[best_index].mu
        return _calibration_decision(
            :TARGET_COMPATIBLE,
            true,
            selected_mu,
            nothing,
            selected_mu,
            nothing,
            sorted_mus,
            reasons,
        )
    end

    bracket_index = findfirst(1:2) do index
        lower_density = sorted_records[index].total_density
        upper_density = sorted_records[index + 1].total_density
        lower_density < target_filling < upper_density
    end
    if isnothing(bracket_index)
        push!(reasons, :target_not_compatible_or_adjacent_mean_bracketed)
        return _calibration_decision(
            :STOP_UNBRACKETED,
            false,
            nothing,
            nothing,
            nothing,
            nothing,
            sorted_mus,
            reasons,
        )
    end

    lower = sorted_records[bracket_index]
    upper = sorted_records[bracket_index + 1]
    unrounded_mu =
        lower.mu +
        (target_filling - lower.total_density) *
        (upper.mu - lower.mu) /
        (upper.total_density - lower.total_density)
    selected_mu = round(unrounded_mu; digits=4)
    return _calibration_decision(
        :INTERPOLATED,
        true,
        selected_mu,
        unrounded_mu,
        nothing,
        (lower.mu, upper.mu),
        sorted_mus,
        reasons,
    )
end

function _two_chain_decision(
    status,
    advance,
    reasons,
    expected_seeds,
    observed_seeds,
    density_compatible_seeds,
    conservative_sign_lower_2se,
)
    return (
        status=status,
        advance=advance,
        reasons=_ordered_unique_reasons(reasons),
        expected_seeds=Tuple(expected_seeds),
        observed_seeds=Tuple(observed_seeds),
        density_compatible_seeds=Tuple(density_compatible_seeds),
        conservative_average_sign_lower_2se=conservative_sign_lower_2se,
    )
end

function _normalized_strict_status(value)
    if value isa Symbol || value isa AbstractString
        return Symbol(uppercase(String(value)))
    end
    return nothing
end

"""
    assess_two_chain_gate(results; expected_seeds,
        target_filling=FIXED_DENSITY_TARGET_FILLING,
        expected_workload=FIXED_DENSITY_PRODUCTION_WORKLOAD)

Adjudicate one registered two-chain campaign point. The expected seed set must
contain exactly two distinct integer seeds, and results must contain exactly
those two records. Each record must be complete, have usable sign
denominators and a valid positive-error density jackknife, and be compatible
with the target at two standard errors. Only two PASS classifications advance.
A valid pair containing any DEFINITIVE_SIGN_FAIL is reported separately;
missing, malformed, density-incompatible, or inconclusive data are not
adjudicable.
"""
function assess_two_chain_gate(
    results;
    expected_seeds,
    target_filling=FIXED_DENSITY_TARGET_FILLING,
    expected_workload=FIXED_DENSITY_PRODUCTION_WORKLOAD,
)
    _validate_expected_calibration_workload(expected_workload)
    records = collect(results)
    expected = collect(expected_seeds)
    reasons = Symbol[]

    if length(expected) != 2
        push!(reasons, :requires_exactly_two_expected_seeds)
    elseif !all(seed -> seed isa Integer && !(seed isa Bool), expected)
        push!(reasons, :expected_seeds_must_be_integers)
    elseif length(unique(expected)) != 2
        push!(reasons, :expected_seeds_must_be_distinct)
    end
    _finite_real(target_filling) ||
        push!(reasons, :target_filling_not_finite)
    if !isempty(reasons)
        return _two_chain_decision(
            :NOT_ADJUDICABLE,
            false,
            reasons,
            expected,
            (),
            (),
            nothing,
        )
    end

    if length(records) != 2
        push!(reasons, :requires_exactly_two_chain_results)
        return _two_chain_decision(
            :NOT_ADJUDICABLE,
            false,
            reasons,
            expected,
            (),
            (),
            nothing,
        )
    end

    required_fields = (
        :seed,
        :completed,
        :total_density,
        :filling_jackknife_standard_error,
        :filling_jackknife_valid,
        :full_sign_denominator_usable,
        :leave_one_out_denominators_usable,
        :average_sign,
        :average_sign_standard_error,
        :strict_pilot_status,
        :burnin_sweeps,
        :measurement_sweeps,
        :sweeps_between_measurements,
        :block_size,
        :stabilization_period,
        :measurement_count,
    )
    for (index, record) in pairs(records)
        for field in required_fields
            hasproperty(record, field) || push!(
                reasons,
                _indexed_reason("chain", index, string("missing_", field)),
            )
        end
    end
    if !isempty(reasons)
        return _two_chain_decision(
            :NOT_ADJUDICABLE,
            false,
            reasons,
            expected,
            (),
            (),
            nothing,
        )
    end

    observed = [record.seed for record in records]
    if !all(seed -> seed isa Integer && !(seed isa Bool), observed)
        push!(reasons, :result_seeds_must_be_integers)
    else
        length(unique(observed)) == 2 ||
            push!(reasons, :result_seeds_must_be_distinct)
        any(seed -> !(seed in expected), observed) &&
            push!(reasons, :unexpected_result_seed)
        any(seed -> !(seed in observed), expected) &&
            push!(reasons, :missing_expected_seed)
    end

    conservative_sign_lower_2se = if all(
        record ->
            _finite_real(record.average_sign) &&
            _nonnegative_finite_real(record.average_sign_standard_error),
        records,
    )
        minimum(
            record.average_sign -
            2 * record.average_sign_standard_error
            for record in records
        )
    else
        nothing
    end

    density_compatible_seeds = Int[]
    normalized_statuses = Union{Nothing,Symbol}[]
    for (index, record) in pairs(records)
        record.completed === true || push!(
            reasons,
            _indexed_reason("chain", index, "not_completed"),
        )
        for field in propertynames(FIXED_DENSITY_PRODUCTION_WORKLOAD)
            getproperty(record, field) ==
                getproperty(expected_workload, field) || push!(
                    reasons,
                    _indexed_reason(
                        "chain",
                        index,
                        string(field, "_mismatch"),
                    ),
                )
        end
        record.measurement_count ==
            expected_workload.measurement_sweeps || push!(
                reasons,
                _indexed_reason(
                    "chain",
                    index,
                    "measurement_count_mismatch",
                ),
            )
        density_valid =
            _finite_real(record.total_density) &&
            _positive_finite_real(
                record.filling_jackknife_standard_error,
            )
        _finite_real(record.total_density) || push!(
            reasons,
            _indexed_reason("chain", index, "density_not_finite"),
        )
        _positive_finite_real(
            record.filling_jackknife_standard_error,
        ) || push!(
            reasons,
            _indexed_reason(
                "chain",
                index,
                "density_standard_error_not_positive_finite",
            ),
        )
        record.filling_jackknife_valid === true || push!(
            reasons,
            _indexed_reason("chain", index, "density_jackknife_not_valid"),
        )
        record.full_sign_denominator_usable === true || push!(
            reasons,
            _indexed_reason("chain", index, "full_denominator_not_usable"),
        )
        record.leave_one_out_denominators_usable === true || push!(
            reasons,
            _indexed_reason(
                "chain",
                index,
                "leave_one_out_denominators_not_usable",
            ),
        )
        _finite_real(record.average_sign) || push!(
            reasons,
            _indexed_reason("chain", index, "average_sign_not_finite"),
        )
        _nonnegative_finite_real(
            record.average_sign_standard_error,
        ) || push!(
            reasons,
            _indexed_reason(
                "chain",
                index,
                "average_sign_standard_error_not_nonnegative_finite",
            ),
        )

        normalized_status = _normalized_strict_status(
            record.strict_pilot_status,
        )
        push!(normalized_statuses, normalized_status)
        normalized_status in (
            :PASS,
            :DEFINITIVE_SIGN_FAIL,
            :INCONCLUSIVE,
        ) || push!(
            reasons,
            _indexed_reason(
                "chain",
                index,
                "strict_pilot_status_unrecognized",
            ),
        )

        if density_valid && record.filling_jackknife_valid === true
            if abs(record.total_density - target_filling) <=
               2 * record.filling_jackknife_standard_error
                record.seed isa Integer &&
                    push!(density_compatible_seeds, record.seed)
            else
                push!(
                    reasons,
                    _indexed_reason(
                        "chain",
                        index,
                        "density_not_target_compatible",
                    ),
                )
            end
        end
    end

    if !isempty(reasons)
        return _two_chain_decision(
            :NOT_ADJUDICABLE,
            false,
            reasons,
            expected,
            observed,
            density_compatible_seeds,
            conservative_sign_lower_2se,
        )
    end

    if all(
        status -> status in (:PASS, :DEFINITIVE_SIGN_FAIL),
        normalized_statuses,
    ) && any(==(:DEFINITIVE_SIGN_FAIL), normalized_statuses)
        return _two_chain_decision(
            :DEFINITIVE_SIGN_FAIL,
            false,
            [:strict_pilot_definitive_sign_fail],
            expected,
            observed,
            density_compatible_seeds,
            conservative_sign_lower_2se,
        )
    end
    if all(==(:PASS), normalized_statuses)
        return _two_chain_decision(
            :ADVANCE,
            true,
            reasons,
            expected,
            observed,
            density_compatible_seeds,
            conservative_sign_lower_2se,
        )
    end
    return _two_chain_decision(
        :NOT_ADJUDICABLE,
        false,
        [:strict_pilot_inconclusive],
        expected,
        observed,
        density_compatible_seeds,
        conservative_sign_lower_2se,
    )
end

function _resource_decision(status, feasible, reasons, observed)
    return (
        status=status,
        feasible=feasible,
        reasons=_ordered_unique_reasons(reasons),
        observed=observed,
        limits=(
            projected_two_chain_seconds=
                RESOURCE_PROJECTED_TWO_CHAIN_SECONDS_LIMIT,
            peak_rss_bytes=RESOURCE_PEAK_RSS_BYTES_LIMIT,
            maximum_greens_recompute_error=
                RESOURCE_GREENS_RECOMPUTE_ERROR_LIMIT,
        ),
    )
end

"""
    assess_resource_gate(record)

Classify a resource benchmark record. Feasibility requires projected time for
two 3600-sweep chains at or below 86400 seconds, peak RSS at or below 5.5 GiB,
zero swap, no OOM or numerical failure, and maximum Green-function correction
at or below 1e-8. Missing or malformed metrics are not adjudicable;
well-formed threshold failures are resource-infeasible.
"""
function assess_resource_gate(record)
    required_fields = (
        :projected_two_chain_seconds,
        :peak_rss_bytes,
        :swap_bytes,
        :oom_occurred,
        :numerical_failure,
        :maximum_greens_recompute_error,
    )
    reasons = Symbol[]
    for field in required_fields
        hasproperty(record, field) ||
            push!(reasons, Symbol(string("missing_", field)))
    end
    if !isempty(reasons)
        return _resource_decision(
            :NOT_ADJUDICABLE,
            false,
            reasons,
            nothing,
        )
    end

    _nonnegative_finite_real(record.projected_two_chain_seconds) || push!(
        reasons,
        :projected_two_chain_seconds_not_nonnegative_finite,
    )
    _nonnegative_finite_real(record.peak_rss_bytes) ||
        push!(reasons, :peak_rss_bytes_not_nonnegative_finite)
    _nonnegative_finite_real(record.swap_bytes) ||
        push!(reasons, :swap_bytes_not_nonnegative_finite)
    record.oom_occurred isa Bool ||
        push!(reasons, :oom_occurred_not_boolean)
    record.numerical_failure isa Bool ||
        push!(reasons, :numerical_failure_not_boolean)
    _nonnegative_finite_real(
        record.maximum_greens_recompute_error,
    ) || push!(
        reasons,
        :maximum_greens_recompute_error_not_nonnegative_finite,
    )

    observed = (
        projected_two_chain_seconds=record.projected_two_chain_seconds,
        peak_rss_bytes=record.peak_rss_bytes,
        swap_bytes=record.swap_bytes,
        oom_occurred=record.oom_occurred,
        numerical_failure=record.numerical_failure,
        maximum_greens_recompute_error=
            record.maximum_greens_recompute_error,
    )
    if !isempty(reasons)
        return _resource_decision(
            :NOT_ADJUDICABLE,
            false,
            reasons,
            observed,
        )
    end

    record.projected_two_chain_seconds <=
        RESOURCE_PROJECTED_TWO_CHAIN_SECONDS_LIMIT || push!(
            reasons,
            :projected_two_chain_time_exceeds_limit,
        )
    record.peak_rss_bytes <= RESOURCE_PEAK_RSS_BYTES_LIMIT ||
        push!(reasons, :peak_rss_exceeds_limit)
    iszero(record.swap_bytes) || push!(reasons, :swap_used)
    record.oom_occurred && push!(reasons, :oom_occurred)
    record.numerical_failure && push!(reasons, :numerical_failure)
    record.maximum_greens_recompute_error <=
        RESOURCE_GREENS_RECOMPUTE_ERROR_LIMIT || push!(
            reasons,
            :greens_recompute_error_exceeds_limit,
        )

    if isempty(reasons)
        return _resource_decision(
            :RESOURCE_FEASIBLE,
            true,
            reasons,
            observed,
        )
    end
    return _resource_decision(
        :RESOURCE_INFEASIBLE,
        false,
        reasons,
        observed,
    )
end

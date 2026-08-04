"""
    _altermagnetic_observables_from_greens(
        greens_up, greens_down, L, configuration_sign,
    )

Return raw altermagnetic observables for one auxiliary-field configuration.
The equal-time convention is `G_sigma[a,b] = <c_sigma,a c_sigma,b^dagger>`,
so `C_sigma = I - G_sigma` and

```text
<n_sigma,a n_sigma,b> =
    C_sigma[a,a] * C_sigma[b,b]
    + delta(a,b) * C_sigma[a,a]
    - C_sigma[a,b] * C_sigma[b,a].
```

The Kronecker-delta term is the same-mode contact contribution.  Opposite
spins factorize for a fixed auxiliary field.  Consequently the returned
`S_alm_q0` is the complete conditional estimator of
`(1/N) * sum(i,j) <O_i O_j>`, not the square of the conditional one-point
estimator.

The values are deliberately not multiplied by `configuration_sign`; a sampler
must sign-weight each raw value exactly once when accumulating measurements.
Orbital-fast indices use odd entries for `d_xz` and even entries for `d_yz`.
"""
function _altermagnetic_observables_from_greens(
    greens_up::AbstractMatrix{<:Real},
    greens_down::AbstractMatrix{<:Real},
    L::Integer,
    configuration_sign::Real,
)
    _, number_of_sites, dimension = _lattice_dimensions(L)
    size(greens_up) == (dimension, dimension) || throw(DimensionMismatch(
        "greens_up must have size (2L^2, 2L^2) = ($dimension, $dimension)",
    ))
    size(greens_down) == (dimension, dimension) || throw(DimensionMismatch(
        "greens_down must have size (2L^2, 2L^2) = ($dimension, $dimension)",
    ))

    up_zero, down_zero, sign_value = promote(
        float(zero(eltype(greens_up))),
        float(zero(eltype(greens_down))),
        float(configuration_sign),
    )
    T = typeof(up_zero)
    unit = one(T)
    half = unit / 2
    occupations_up = Vector{T}(undef, dimension)
    occupations_down = Vector{T}(undef, dimension)
    for index in 1:dimension
        occupations_up[index] = unit - convert(T, greens_up[index, index])
        occupations_down[index] = unit - convert(T, greens_down[index, index])
    end

    magnetization_x_sum = zero(T)
    magnetization_y_sum = zero(T)
    for site in 0:(number_of_sites - 1)
        index_x = 2 * site + 1
        index_y = index_x + 1
        magnetization_x_sum +=
            occupations_up[index_x] - occupations_down[index_x]
        magnetization_y_sum +=
            occupations_up[index_y] - occupations_down[index_y]
    end

    number_of_sites_typed = convert(T, number_of_sites)
    mx = magnetization_x_sum / number_of_sites_typed
    my = magnetization_y_sum / number_of_sites_typed
    m_fm = half * (mx + my)
    delta_alm = half * (mx - my)
    total_alm_mean = number_of_sites_typed * delta_alm

    total_alm_covariance = zero(T)
    for row in 1:dimension, column in 1:dimension
        orbital_weight_row = isodd(row) ? half : -half
        orbital_weight_column = isodd(column) ? half : -half

        density_up_row_column =
            (row == column ? unit : zero(T)) -
            convert(T, greens_up[row, column])
        density_up_column_row =
            (row == column ? unit : zero(T)) -
            convert(T, greens_up[column, row])
        density_down_row_column =
            (row == column ? unit : zero(T)) -
            convert(T, greens_down[row, column])
        density_down_column_row =
            (row == column ? unit : zero(T)) -
            convert(T, greens_down[column, row])

        contact_up = row == column ? occupations_up[row] : zero(T)
        contact_down = row == column ? occupations_down[row] : zero(T)
        total_alm_covariance += orbital_weight_row * orbital_weight_column * (
            contact_up - density_up_row_column * density_up_column_row +
            contact_down - density_down_row_column * density_down_column_row
        )
    end

    total_alm_second_moment = total_alm_mean^2 + total_alm_covariance
    S_alm_q0 = total_alm_second_moment / number_of_sites_typed
    m_alm_squared = S_alm_q0 / number_of_sites_typed
    return (
        mx=mx,
        my=my,
        m_fm=m_fm,
        delta_alm=delta_alm,
        S_alm_q0=S_alm_q0,
        m_alm_squared=m_alm_squared,
        configuration_sign=convert(T, sign_value),
    )
end

"""
    measure_altermagnetic_observables(state)

Measure the raw Eq. 4 altermagnetic one- and two-point observables for the
current DQMC auxiliary-field configuration.  With

`O_i = ((n_ix_up - n_ix_down) - (n_iy_up - n_iy_down)) / 2`,

the result contains the per-cell orbital magnetizations `mx` and `my`,
`m_fm = (mx + my)/2`, `delta_alm = (mx - my)/2`,
`S_alm_q0 = sum(i,j)<O_i O_j>/N`, `m_alm_squared = S_alm_q0/N`, and the
configuration sign.  Observable values are raw conditional estimators and
must be sign-weighted exactly once by the caller.
"""
function measure_altermagnetic_observables(state::DQMCState)
    return _altermagnetic_observables_from_greens(
        state.greens_up,
        state.greens_down,
        state.parameters.L,
        _configuration_sign(state),
    )
end

const _ALMAGNETIC_REWEIGHTING_FIELDS = (
    :mx,
    :my,
    :m_fm,
    :delta_alm,
    :S_alm_q0,
    :m_alm_squared,
)

"""
    summarize_sign_reweighted_altermagnetism(
        configuration_signs, raw_observables;
        block_size, minimum_absolute_average_sign=0.0,
    )

Sign reweight the six physical fields returned by
`measure_altermagnetic_observables`. Every raw configuration estimator is
multiplied by its separately supplied configuration sign exactly once. The
central estimate for observable `a` is `sum(s .* O_a) / sum(s)`, never an
average of block ratios. The `configuration_sign` field in a raw observation
is not used as an additional weight.

All observables use common contiguous blocks and the same full and
leave-one-block-out sign denominators. Standard errors use the jackknife
ratios `(X_a - X_b,a) / (S - S_b)`. A final short block is retained. A
denominator is usable only when its average sign is finite and has absolute
value strictly larger than `minimum_absolute_average_sign`; unusable central
or jackknife denominators produce `missing` estimates or errors and false
validity flags.
"""
function summarize_sign_reweighted_altermagnetism(
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
        for name in _ALMAGNETIC_REWEIGHTING_FIELDS
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
    number_of_observables = length(_ALMAGNETIC_REWEIGHTING_FIELDS)
    raw_values = Matrix{T}(undef, measurement_count, number_of_observables)
    for measurement in 1:measurement_count,
        (observable, name) in enumerate(_ALMAGNETIC_REWEIGHTING_FIELDS)
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
        mx=statistics.estimates[1],
        mx_standard_error=jackknife.standard_errors[1],
        mx_jackknife_valid=jackknife.valid[1],
        my=statistics.estimates[2],
        my_standard_error=jackknife.standard_errors[2],
        my_jackknife_valid=jackknife.valid[2],
        m_fm=statistics.estimates[3],
        m_fm_standard_error=jackknife.standard_errors[3],
        m_fm_jackknife_valid=jackknife.valid[3],
        delta_alm=statistics.estimates[4],
        delta_alm_standard_error=jackknife.standard_errors[4],
        delta_alm_jackknife_valid=jackknife.valid[4],
        S_alm_q0=statistics.estimates[5],
        S_alm_q0_standard_error=jackknife.standard_errors[5],
        S_alm_q0_jackknife_valid=jackknife.valid[5],
        m_alm_squared=statistics.estimates[6],
        m_alm_squared_standard_error=jackknife.standard_errors[6],
        m_alm_squared_jackknife_valid=jackknife.valid[6],
        configuration_signs_finite=configuration_signs_finite,
        configuration_signs_valid=configuration_signs_valid,
        raw_observables_finite=raw_observables_finite,
        finite_results=finite_results,
    )
end

using AltermagnetDQMCTrial

const THERMODYNAMIC_NAMES = (
    :density,
    :doping,
    :kinetic_energy_per_cell,
    :shifted_interaction_energy_per_cell,
    :internal_energy_per_cell,
    :grand_energy_per_cell,
)

const L = 2
const BETA = 0.2
const DTAU = 0.1
const INTERACTION = 2.0
const CHEMICAL_POTENTIAL = -1.897
const STABILIZATION_PERIOD = 2
const WARMUP_SWEEPS = 200
const BIN_COUNT = 32
const MEASUREMENTS_PER_BIN = 100
const MEASUREMENT_SWEEPS = BIN_COUNT * MEASUREMENTS_PER_BIN
const SWEEPS_BETWEEN_MEASUREMENTS = 1
const SEED = Int(0x5A17C0DE)

const EXPECTED_SINGLE_PARTICLE_DIMENSION = 8
const EXPECTED_TIME_SLICES = 2
const EXPECTED_FIELD_COUNT = 16
const EXPECTED_CONFIGURATION_COUNT = 65_536
const EXPECTED_COMPLETE_SWEEPS = WARMUP_SWEEPS + MEASUREMENT_SWEEPS
const EXPECTED_ATTEMPTED_UPDATES =
    EXPECTED_COMPLETE_SWEEPS * EXPECTED_FIELD_COUNT
const EXPECTED_DIRECTION_SWEEPS = div(EXPECTED_COMPLETE_SWEEPS, 2)
const EXPECTED_DIRECTION_ATTEMPTS =
    EXPECTED_DIRECTION_SWEEPS * EXPECTED_FIELD_COUNT

const GREEN_CORRECTION_TOLERANCE = 1.0e-9
const THERMODYNAMIC_IDENTITY_TOLERANCE = 1.0e-10
const STATISTICAL_FLOOR = 1.0e-10
const ALL_POSITIVE_SIGN_TOLERANCE = 1.0e-12

function require_condition(condition::Bool, message::AbstractString)
    condition || error(message)
    return nothing
end

function fields_from_mask(
    mask::Integer,
    single_particle_dimension::Int,
    time_slice_count::Int,
)
    field_count = Base.checked_mul(
        single_particle_dimension,
        time_slice_count,
    )
    configuration_count = 1 << field_count
    0 <= mask < configuration_count || throw(ArgumentError(
        "HS mask must lie in 0:$(configuration_count - 1)",
    ))

    fields = Matrix{Int8}(
        undef,
        single_particle_dimension,
        time_slice_count,
    )
    # Julia linear indexing is column-major: bits 0:7 fill slice 1 and bits
    # 8:15 fill slice 2. Every bit index in 0:(field_count-1) occurs once.
    for bit_index in 0:(field_count - 1)
        bit = (mask >> bit_index) & 1
        fields[bit_index + 1] = iszero(bit) ? Int8(-1) : Int8(1)
    end
    return fields
end

function verify_field_decoder(dimensions)
    field_count = dimensions.Nsp * dimensions.Ltau
    require_condition(
        field_count == EXPECTED_FIELD_COUNT,
        "the decoder audit expected exactly 16 HS variables",
    )

    zero_mask_fields = fields_from_mask(0, dimensions.Nsp, dimensions.Ltau)
    full_mask_fields = fields_from_mask(
        (1 << field_count) - 1,
        dimensions.Nsp,
        dimensions.Ltau,
    )
    require_condition(
        all(==(-Int8(1)), zero_mask_fields),
        "the zero HS mask did not decode entirely to -1",
    )
    require_condition(
        all(==(Int8(1)), full_mask_fields),
        "the full HS mask did not decode entirely to +1",
    )

    visited_positions = falses(field_count)
    for bit_index in 0:(field_count - 1)
        one_bit_fields = fields_from_mask(
            1 << bit_index,
            dimensions.Nsp,
            dimensions.Ltau,
        )
        positive_positions = findall(==(Int8(1)), one_bit_fields)
        require_condition(
            length(positive_positions) == 1,
            "HS bit $bit_index did not select exactly one field",
        )
        linear_index = LinearIndices(one_bit_fields)[only(positive_positions)]
        require_condition(
            linear_index == bit_index + 1,
            "HS bit $bit_index decoded to the wrong field",
        )
        require_condition(
            !visited_positions[linear_index],
            "HS field position $linear_index was decoded more than once",
        )
        visited_positions[linear_index] = true
    end
    require_condition(
        all(visited_positions),
        "the HS decoder did not visit all 16 field positions",
    )

    println("field_decoder_verified_bits=", count(visited_positions))
    println("field_decoder_layout=bits_0_to_7_slice_1_bits_8_to_15_slice_2")
    return nothing
end

function thermodynamic_named_tuple(values)
    return NamedTuple{THERMODYNAMIC_NAMES}(Tuple(values))
end

function exact_auxiliary_field_average(parameters, dimensions)
    field_count = dimensions.Nsp * dimensions.Ltau
    configuration_count = 1 << field_count
    log_abs_weights = Vector{Float64}(undef, configuration_count)
    configuration_signs = Vector{Float64}(undef, configuration_count)
    raw_observations = Vector{NamedTuple}(undef, configuration_count)
    finite_configuration_count = 0
    valid_sign_count = 0
    matching_sign_count = 0
    positive_weight_count = 0

    # This positive, field-independent discrete-HS prefactor cancels from
    # normalized ratios. Including it records the complete finite-dtau
    # auxiliary-field weight, not a Hamiltonian exact-diagonalization weight.
    log_hs_prefactor = field_count * (
        -log(2.0) - parameters.dtau * parameters.U / 4
    )

    for mask in 0:(configuration_count - 1)
        fields = fields_from_mask(mask, dimensions.Nsp, dimensions.Ltau)
        state = initialize_dqmc_state(parameters; hs_fields=fields)

        determinant_product_sign = Float64(sign(
            inv(state.signdet_green_up) * inv(state.signdet_green_down),
        ))
        recorded_sign = measure_basic_observables(state).configuration_sign
        log_abs_fermion_weight =
            -state.logdet_green_up - state.logdet_green_down
        log_abs_weight = log_abs_fermion_weight + log_hs_prefactor
        absolute_weight = exp(log_abs_weight)
        observation = measure_thermodynamic_observables(state)

        sign_valid = determinant_product_sign in (-1.0, 1.0)
        sign_matches = recorded_sign == determinant_product_sign
        positive_weight = isfinite(absolute_weight) && absolute_weight > 0.0
        finite_configuration =
            isfinite(log_abs_weight) &&
            positive_weight &&
            all(isfinite, state.greens_up) &&
            all(isfinite, state.greens_down) &&
            all(
                name -> isfinite(getproperty(observation, name)),
                THERMODYNAMIC_NAMES,
            )

        require_condition(
            sign_valid,
            "invalid determinant-product sign for HS mask $mask",
        )
        require_condition(
            sign_matches,
            "recorded sign disagrees with determinant-product sign for mask $mask",
        )
        require_condition(
            positive_weight,
            "absolute sampling weight is not finite and positive for mask $mask",
        )
        require_condition(
            finite_configuration,
            "non-finite exact Green matrix, weight, or observable for mask $mask",
        )

        index = mask + 1
        log_abs_weights[index] = log_abs_weight
        configuration_signs[index] = determinant_product_sign
        raw_observations[index] = observation
        valid_sign_count += 1
        matching_sign_count += 1
        positive_weight_count += 1
        finite_configuration_count += 1
    end

    log_weight_shift = maximum(log_abs_weights)
    scaled_abs_weights = exp.(log_abs_weights .- log_weight_shift)
    require_condition(
        all(weight -> isfinite(weight) && weight > 0.0, scaled_abs_weights),
        "a scaled exact absolute weight is not finite and positive",
    )
    absolute_denominator = sum(scaled_abs_weights)
    signed_denominator = sum(configuration_signs .* scaled_abs_weights)
    positive_absolute_weight = mapreduce(
        identity,
        +,
        (
            scaled_abs_weights[index]
            for index in eachindex(configuration_signs)
            if configuration_signs[index] == 1.0
        );
        init=zero(eltype(scaled_abs_weights)),
    )
    negative_absolute_weight = mapreduce(
        identity,
        +,
        (
            scaled_abs_weights[index]
            for index in eachindex(configuration_signs)
            if configuration_signs[index] == -1.0
        );
        init=zero(eltype(scaled_abs_weights)),
    )
    require_condition(
        isfinite(absolute_denominator) && absolute_denominator > 0.0,
        "the exact absolute-weight denominator is unusable",
    )
    require_condition(
        isfinite(signed_denominator) && !iszero(signed_denominator),
        "the exact signed-weight denominator is unusable",
    )
    require_condition(
        isfinite(positive_absolute_weight) &&
        isfinite(negative_absolute_weight) &&
        positive_absolute_weight >= 0.0 &&
        negative_absolute_weight >= 0.0,
        "an exact sign-sector absolute weight is invalid",
    )

    exact_average_sign = signed_denominator / absolute_denominator
    negative_absolute_weight_probability =
        negative_absolute_weight / absolute_denominator
    require_condition(
        isfinite(exact_average_sign) && -1.0 <= exact_average_sign <= 1.0,
        "the exact average sign is invalid",
    )
    require_condition(
        isfinite(negative_absolute_weight_probability) &&
        0.0 <= negative_absolute_weight_probability <= 1.0,
        "the exact negative absolute-weight probability is invalid",
    )
    require_condition(
        abs(
            negative_absolute_weight_probability -
            (1.0 - exact_average_sign) / 2,
        ) <= 1.0e-12,
        "the two exact negative-weight probability formulas disagree",
    )

    exact_values = thermodynamic_named_tuple(
        sum(
            configuration_signs[index] * scaled_abs_weights[index] *
            getproperty(raw_observations[index], name)
            for index in eachindex(configuration_signs)
        ) / signed_denominator
        for name in THERMODYNAMIC_NAMES
    )
    require_condition(
        all(
            name -> isfinite(getproperty(exact_values, name)),
            THERMODYNAMIC_NAMES,
        ),
        "an exact sign-reweighted observable is non-finite",
    )

    return (
        values=exact_values,
        average_sign=exact_average_sign,
        negative_absolute_weight_probability=
            negative_absolute_weight_probability,
        configuration_count=configuration_count,
        finite_configuration_count=finite_configuration_count,
        valid_sign_count=valid_sign_count,
        matching_sign_count=matching_sign_count,
        positive_weight_count=positive_weight_count,
        positive_sign_count=count(==(1.0), configuration_signs),
        negative_sign_count=count(==(-1.0), configuration_signs),
        minimum_log_abs_weight=minimum(log_abs_weights),
        maximum_log_abs_weight=maximum(log_abs_weights),
        log_weight_shift=log_weight_shift,
        scaled_positive_absolute_weight=positive_absolute_weight,
        scaled_negative_absolute_weight=negative_absolute_weight,
        scaled_absolute_denominator=absolute_denominator,
        scaled_signed_denominator=signed_denominator,
    )
end

function audited_complete_sweep!(state, dimensions)
    up = state.greens_calculator_up
    down = state.greens_calculator_down
    require_condition(
        up.forward == down.forward && up.l == down.l,
        "spin calculators are not synchronized before a complete sweep",
    )
    require_condition(
        up.Lτ == dimensions.Ltau && down.Lτ == dimensions.Ltau,
        "a Green calculator has the wrong imaginary-time length",
    )

    direction_before = up.forward
    expected_start = direction_before ? 1 : dimensions.Ltau
    require_condition(
        up.l == expected_start,
        "a complete sweep did not start at its expected boundary",
    )

    attempted_before = state.attempted_updates
    accepted_before = state.accepted_updates
    sweep = sweep_all_time_slices!(state)
    expected_attempts = dimensions.Nsp * dimensions.Ltau

    require_condition(
        up.forward == down.forward && up.l == down.l,
        "spin calculators are not synchronized after a complete sweep",
    )
    require_condition(
        up.forward == !direction_before,
        "a complete sweep did not toggle traversal direction",
    )
    expected_end = direction_before ? dimensions.Ltau : 1
    require_condition(
        up.l == expected_end,
        "a complete sweep did not reach its expected end boundary",
    )
    require_condition(
        sweep.attempted == expected_attempts &&
        state.attempted_updates - attempted_before == expected_attempts,
        "a complete sweep attempted the wrong number of local updates",
    )
    require_condition(
        state.accepted_updates - accepted_before == sweep.accepted,
        "a complete sweep reported the wrong accepted-update count",
    )

    return (result=sweep, forward=direction_before)
end

function run_production_chain(parameters, dimensions)
    # No fixed field is supplied: initialization, shuffles, and acceptance
    # draws consume one seeded production RNG trajectory.
    state = initialize_dqmc_state(parameters)
    forward_sweep_count = 0
    reverse_sweep_count = 0
    forward_attempted_updates = 0
    reverse_attempted_updates = 0

    for _ in 1:WARMUP_SWEEPS
        audited = audited_complete_sweep!(state, dimensions)
        if audited.forward
            forward_sweep_count += 1
            forward_attempted_updates += audited.result.attempted
        else
            reverse_sweep_count += 1
            reverse_attempted_updates += audited.result.attempted
        end
    end
    warmup_attempted_updates = state.attempted_updates
    warmup_accepted_updates = state.accepted_updates

    configuration_signs = Vector{Float64}(undef, MEASUREMENT_SWEEPS)
    raw_observations = Vector{NamedTuple}(undef, MEASUREMENT_SWEEPS)
    measurements_accumulated = 0
    for measurement in 1:MEASUREMENT_SWEEPS
        for _ in 1:SWEEPS_BETWEEN_MEASUREMENTS
            audited = audited_complete_sweep!(state, dimensions)
            if audited.forward
                forward_sweep_count += 1
                forward_attempted_updates += audited.result.attempted
            else
                reverse_sweep_count += 1
                reverse_attempted_updates += audited.result.attempted
            end
        end

        recorded_sign = measure_basic_observables(state).configuration_sign
        observation = measure_thermodynamic_observables(state)
        require_condition(
            isfinite(recorded_sign) && recorded_sign in (-1.0, 1.0),
            "a sampled configuration sign is invalid",
        )
        require_condition(
            all(
                name -> isfinite(getproperty(observation, name)),
                THERMODYNAMIC_NAMES,
            ),
            "a sampled raw thermodynamic observable is non-finite",
        )
        configuration_signs[measurement] = recorded_sign
        raw_observations[measurement] = observation
        measurements_accumulated += 1
    end

    # Direction counts are complete here. The explicit recomputation below
    # iterates the calculators internally and is intentionally not counted as
    # a Metropolis sweep.
    final_recomputation = recompute_equal_time_greens!(state)
    acceptance_rate = state.accepted_updates / state.attempted_updates
    summary = summarize_sign_reweighted_thermodynamics(
        configuration_signs,
        raw_observations;
        block_size=MEASUREMENTS_PER_BIN,
        minimum_absolute_average_sign=0.0,
    )

    return (
        state=state,
        summary=summary,
        configuration_signs=configuration_signs,
        raw_observations=raw_observations,
        measurements_accumulated=measurements_accumulated,
        warmup_attempted_updates=warmup_attempted_updates,
        warmup_accepted_updates=warmup_accepted_updates,
        measurement_attempted_updates=
            state.attempted_updates - warmup_attempted_updates,
        measurement_accepted_updates=
            state.accepted_updates - warmup_accepted_updates,
        acceptance_rate=acceptance_rate,
        final_recomputation=final_recomputation,
        forward_sweep_count=forward_sweep_count,
        reverse_sweep_count=reverse_sweep_count,
        forward_attempted_updates=forward_attempted_updates,
        reverse_attempted_updates=reverse_attempted_updates,
    )
end

function thermodynamic_identity_errors(values, mu)
    return (
        doping=abs(values.doping - (values.density - 2) / 4),
        internal=abs(
            values.internal_energy_per_cell -
            values.kinetic_energy_per_cell -
            values.shifted_interaction_energy_per_cell,
        ),
        grand=abs(
            values.grand_energy_per_cell -
            values.internal_energy_per_cell + mu * values.density,
        ),
    )
end

function comparison_rows(exact_values, dqmc_summary)
    rows = NamedTuple[]
    for name in THERMODYNAMIC_NAMES
        exact_value = getproperty(exact_values, name)
        estimate = getproperty(dqmc_summary, name)
        standard_error = getproperty(
            dqmc_summary,
            Symbol(name, "_standard_error"),
        )
        jackknife_valid = getproperty(
            dqmc_summary,
            Symbol(name, "_jackknife_valid"),
        )
        require_condition(
            jackknife_valid &&
            !ismissing(estimate) &&
            !ismissing(standard_error),
            "missing or invalid DQMC ratio-jackknife result for $name",
        )
        require_condition(
            isfinite(estimate) &&
            isfinite(standard_error) &&
            standard_error > 0.0,
            "DQMC ratio-jackknife SE for $name is not finite and positive",
        )

        absolute_deviation = abs(estimate - exact_value)
        se_normalized_deviation = absolute_deviation / standard_error
        allowed_deviation = max(5 * standard_error, STATISTICAL_FLOOR)
        push!(rows, (
            observable=name,
            exact=exact_value,
            dqmc=estimate,
            standard_error=standard_error,
            absolute_deviation=absolute_deviation,
            se_normalized_deviation=se_normalized_deviation,
            allowed_deviation=allowed_deviation,
            passed=absolute_deviation <= allowed_deviation,
        ))
    end
    return rows
end

function validate_average_sign(exact, run)
    summary = run.summary
    sign_standard_error = summary.average_sign_standard_error
    require_condition(
        !ismissing(sign_standard_error) &&
        isfinite(summary.average_sign) &&
        isfinite(sign_standard_error) &&
        sign_standard_error >= 0.0,
        "the sampled average sign or blocked sign SE is invalid",
    )

    sampled_positive_sign_count = count(==(1.0), run.configuration_signs)
    sampled_negative_sign_count = count(==(-1.0), run.configuration_signs)
    println("sampled_positive_sign_count=", sampled_positive_sign_count)
    println("sampled_negative_sign_count=", sampled_negative_sign_count)

    if exact.negative_sign_count == 0
        require_condition(
            exact.positive_sign_count == EXPECTED_CONFIGURATION_COUNT,
            "the exact all-positive branch has the wrong positive-sign count",
        )
        require_condition(
            abs(exact.average_sign - 1.0) <= ALL_POSITIVE_SIGN_TOLERANCE,
            "the all-positive exact average sign is not one",
        )
        require_condition(
            sampled_positive_sign_count == MEASUREMENT_SWEEPS &&
            sampled_negative_sign_count == 0,
            "an all-positive exact ensemble produced a negative sampled sign",
        )
        require_condition(
            abs(summary.average_sign - 1.0) <=
                ALL_POSITIVE_SIGN_TOLERANCE,
            "the all-positive sampled average sign is not one",
        )
        require_condition(
            iszero(sign_standard_error),
            "the all-positive sampled sign SE is not zero",
        )
        println("average_sign_gate_branch=all_exact_configurations_positive")
        println("average_sign_allowed_deviation=", ALL_POSITIVE_SIGN_TOLERANCE)
        println("average_sign_gate_passed=true")
        return (
            branch=:all_positive,
            sampled_positive_sign_count=sampled_positive_sign_count,
            sampled_negative_sign_count=sampled_negative_sign_count,
            absolute_deviation=abs(summary.average_sign - exact.average_sign),
            allowed_deviation=ALL_POSITIVE_SIGN_TOLERANCE,
            passed=true,
        )
    end

    require_condition(
        exact.positive_sign_count > 0 && exact.negative_sign_count > 0,
        "the exact sign ensemble is neither all-positive nor mixed-sign",
    )
    println("average_sign_gate_branch=mixed_exact_signs")
    println(
        "exact_negative_absolute_weight_probability=",
        exact.negative_absolute_weight_probability,
    )
    if sampled_positive_sign_count == 0 || sampled_negative_sign_count == 0
        println("average_sign_gate_status=STATISTICALLY_INCONCLUSIVE")
        error(
            "the exact ensemble has both signs but the sampled chain has only one",
        )
    end
    require_condition(
        sign_standard_error > 0.0,
        "the mixed-sign blocked sign SE is not strictly positive",
    )
    absolute_deviation = abs(summary.average_sign - exact.average_sign)
    allowed_deviation = max(5 * sign_standard_error, STATISTICAL_FLOOR)
    require_condition(
        absolute_deviation <= allowed_deviation,
        "the sampled average sign disagrees with exact enumeration",
    )
    println("average_sign_allowed_deviation=", allowed_deviation)
    println("average_sign_gate_passed=true")
    return (
        branch=:mixed,
        sampled_positive_sign_count=sampled_positive_sign_count,
        sampled_negative_sign_count=sampled_negative_sign_count,
        absolute_deviation=absolute_deviation,
        allowed_deviation=allowed_deviation,
        passed=true,
    )
end

function print_parameters(parameters, dimensions, model_parameters)
    field_count = dimensions.Nsp * dimensions.Ltau
    println("L=", parameters.L)
    println("boundary_conditions=PBC")
    println("orbital_basis=(d_xz,d_yz)")
    println("t1_over_t=", model_parameters.t1)
    println("t2_over_t=", model_parameters.t2)
    println("t3_over_t=", model_parameters.t3)
    println("t4_over_t=", model_parameters.t4)
    println("beta_t=", parameters.beta)
    println("delta_tau_t=", parameters.dtau)
    println("L_tau=", dimensions.Ltau)
    println("U_over_t=", parameters.U)
    println("mu_over_t=", parameters.mu)
    println("matrix_dimension=", dimensions.Nsp)
    println("hs_variable_count=", field_count)
    println("exact_configuration_count=", 1 << field_count)
    println("seed_hex=0x5A17C0DE")
    println("seed_decimal=", parameters.seed)
    println("stabilization_period_slices=", parameters.stabilization_period)
    println("warmup_complete_sweeps=", WARMUP_SWEEPS)
    println("sweeps_between_measurements=", SWEEPS_BETWEEN_MEASUREMENTS)
    println("contiguous_bin_count=", BIN_COUNT)
    println("measurements_per_bin=", MEASUREMENTS_PER_BIN)
    println("total_measurements=", MEASUREMENT_SWEEPS)
    println("expected_complete_sweeps=", EXPECTED_COMPLETE_SWEEPS)
    println("expected_attempted_updates=", EXPECTED_ATTEMPTED_UPDATES)
    println("production_average_sign_threshold_applied=false")
end

function print_exact_diagnostics(exact)
    println("exact_finite_configuration_count=", exact.finite_configuration_count)
    println("exact_valid_sign_count=", exact.valid_sign_count)
    println("exact_matching_sign_count=", exact.matching_sign_count)
    println("exact_positive_weight_count=", exact.positive_weight_count)
    println("exact_positive_sign_count=", exact.positive_sign_count)
    println("exact_negative_sign_count=", exact.negative_sign_count)
    println("exact_minimum_log_abs_weight=", exact.minimum_log_abs_weight)
    println("exact_maximum_log_abs_weight=", exact.maximum_log_abs_weight)
    println("exact_log_weight_shift=", exact.log_weight_shift)
    println(
        "exact_scaled_positive_absolute_weight=",
        exact.scaled_positive_absolute_weight,
    )
    println(
        "exact_scaled_negative_absolute_weight=",
        exact.scaled_negative_absolute_weight,
    )
    println(
        "exact_scaled_absolute_denominator=",
        exact.scaled_absolute_denominator,
    )
    println(
        "exact_scaled_signed_denominator=",
        exact.scaled_signed_denominator,
    )
    println("exact_average_sign=", exact.average_sign)
    println(
        "exact_negative_absolute_weight_probability=",
        exact.negative_absolute_weight_probability,
    )
end

function print_sampling_diagnostics(run)
    state = run.state
    summary = run.summary
    println("measurements_accumulated=", run.measurements_accumulated)
    println("warmup_attempted_updates=", run.warmup_attempted_updates)
    println("warmup_accepted_updates=", run.warmup_accepted_updates)
    println("measurement_attempted_updates=", run.measurement_attempted_updates)
    println("measurement_accepted_updates=", run.measurement_accepted_updates)
    println("attempted_updates=", state.attempted_updates)
    println("accepted_updates=", state.accepted_updates)
    println("rejected_updates=", state.attempted_updates - state.accepted_updates)
    println("acceptance_rate=", run.acceptance_rate)
    println("average_sign=", summary.average_sign)
    println("average_sign_blocked_standard_error=", summary.average_sign_standard_error)
    println("block_count=", summary.block_count)
    println("minimum_block_length=", summary.minimum_block_length)
    println("maximum_block_length=", summary.maximum_block_length)
    println("configuration_sign_sum=", summary.configuration_sign_sum)
    println("sign_denominator_usable=", summary.sign_denominator_usable)
    println(
        "leave_one_out_denominators_usable=",
        summary.leave_one_out_denominators_usable,
    )
    println(
        "leave_one_out_denominator_failure_count=",
        summary.leave_one_out_denominator_failure_count,
    )
    println(
        "minimum_absolute_leave_one_out_average_sign=",
        summary.minimum_absolute_leave_one_out_average_sign,
    )
    println("denominator_guard_passed=", summary.denominator_guard_passed)
    println("forward_complete_sweeps=", run.forward_sweep_count)
    println("reverse_complete_sweeps=", run.reverse_sweep_count)
    println("forward_attempted_updates=", run.forward_attempted_updates)
    println("reverse_attempted_updates=", run.reverse_attempted_updates)
    println("both_sweep_directions_verified=true")
    println(
        "final_green_recomputation_error=",
        run.final_recomputation.maximum,
    )
    println(
        "maximum_greens_recompute_error=",
        state.maximum_greens_recompute_error,
    )
end

function print_identity_diagnostics(exact_errors, sampled_errors)
    for name in propertynames(exact_errors)
        println(
            "exact_thermodynamic_identity_error_", name, "=",
            getproperty(exact_errors, name),
        )
        println(
            "sampled_thermodynamic_identity_error_", name, "=",
            getproperty(sampled_errors, name),
        )
    end
    println("maximum_exact_thermodynamic_identity_error=", maximum(exact_errors))
    println(
        "maximum_sampled_thermodynamic_identity_error=",
        maximum(sampled_errors),
    )
end

function print_comparison(rows)
    println(
        "observable\texact\tdqmc\tratio_jackknife_se\tabs_deviation\t" *
        "se_normalized_deviation\tallowed_max_5se_or_1e-10\tpassed",
    )
    for row in rows
        println(
            row.observable, '\t',
            row.exact, '\t',
            row.dqmc, '\t',
            row.standard_error, '\t',
            row.absolute_deviation, '\t',
            row.se_normalized_deviation, '\t',
            row.allowed_deviation, '\t',
            row.passed,
        )
    end
end

function run_validation()
    parameters = DQMCParameters(
        L=L,
        beta=BETA,
        dtau=DTAU,
        U=INTERACTION,
        mu=CHEMICAL_POTENTIAL,
        stabilization_period=STABILIZATION_PERIOD,
        burnin_sweeps=WARMUP_SWEEPS,
        measurement_sweeps=MEASUREMENT_SWEEPS,
        sweeps_between_measurements=SWEEPS_BETWEEN_MEASUREMENTS,
        seed=SEED,
    )
    dimensions = validate_dqmc_parameters(parameters)
    model_parameters = DEFAULT_PARAMETERS

    require_condition(parameters.L == L, "this validation is restricted to L=2")
    require_condition(
        dimensions.Nsp == EXPECTED_SINGLE_PARTICLE_DIMENSION,
        "this validation requires an 8-dimensional single-spin matrix",
    )
    require_condition(
        dimensions.Ltau == EXPECTED_TIME_SLICES,
        "this validation is restricted to L_tau=2",
    )
    require_condition(
        parameters.U == INTERACTION,
        "this validation is restricted to U/t=2",
    )
    require_condition(
        parameters.mu == CHEMICAL_POTENTIAL,
        "this validation is restricted to mu/t=-1.897",
    )
    require_condition(
        parameters.dtau == DTAU && parameters.beta == BETA,
        "this validation has the wrong imaginary-time grid",
    )
    require_condition(
        parameters.stabilization_period == STABILIZATION_PERIOD,
        "this validation requires stabilization cadence 2",
    )
    require_condition(
        parameters.burnin_sweeps == WARMUP_SWEEPS &&
        parameters.measurement_sweeps == MEASUREMENT_SWEEPS &&
        parameters.sweeps_between_measurements ==
            SWEEPS_BETWEEN_MEASUREMENTS,
        "this validation has the wrong fixed sampling workload",
    )
    require_condition(
        parameters.seed == SEED,
        "this validation has the wrong fixed seed",
    )
    require_condition(
        (model_parameters.t1, model_parameters.t2,
         model_parameters.t3, model_parameters.t4) ==
        (-1.0, -1.75, -0.85, -0.65),
        "the hopping parameters differ from the authorized model",
    )

    println("validation_name=L2_Ltau2_U2_exact_HS_vs_Metropolis")
    println(
        "oracle_scope=exact enumeration of the finite-delta-tau " *
        "auxiliary-field ensemble",
    )
    println("oracle_is_hamiltonian_exact_diagonalization=false")
    println("output_mode=print_only")
    print_parameters(parameters, dimensions, model_parameters)
    verify_field_decoder(dimensions)

    exact_started = time_ns()
    exact = exact_auxiliary_field_average(parameters, dimensions)
    exact_runtime_seconds = (time_ns() - exact_started) / 1.0e9
    require_condition(
        exact.configuration_count == EXPECTED_CONFIGURATION_COUNT,
        "exact enumeration did not contain 65,536 configurations",
    )
    require_condition(
        exact.finite_configuration_count == EXPECTED_CONFIGURATION_COUNT,
        "not all 65,536 exact configurations were finite",
    )
    require_condition(
        exact.valid_sign_count == EXPECTED_CONFIGURATION_COUNT,
        "not all exact determinant-product signs were valid",
    )
    require_condition(
        exact.matching_sign_count == EXPECTED_CONFIGURATION_COUNT,
        "not all exact recorded signs matched determinant-product signs",
    )
    require_condition(
        exact.positive_weight_count == EXPECTED_CONFIGURATION_COUNT,
        "not all exact absolute sampling weights were positive",
    )
    require_condition(
        exact.positive_sign_count + exact.negative_sign_count ==
            EXPECTED_CONFIGURATION_COUNT,
        "exact positive and negative sign counts do not sum to 65,536",
    )
    print_exact_diagnostics(exact)
    println("exact_enumeration_runtime_seconds=", exact_runtime_seconds)

    chain_started = time_ns()
    run = run_production_chain(parameters, dimensions)
    metropolis_runtime_seconds = (time_ns() - chain_started) / 1.0e9
    summary = run.summary
    expected_warmup_attempts = WARMUP_SWEEPS * EXPECTED_FIELD_COUNT
    expected_measurement_attempts =
        MEASUREMENT_SWEEPS * SWEEPS_BETWEEN_MEASUREMENTS *
        EXPECTED_FIELD_COUNT

    require_condition(
        run.measurements_accumulated == MEASUREMENT_SWEEPS &&
        summary.measurement_count == MEASUREMENT_SWEEPS,
        "the production chain did not accumulate exactly 3,200 measurements",
    )
    require_condition(
        summary.block_count == BIN_COUNT &&
        summary.minimum_block_length == MEASUREMENTS_PER_BIN &&
        summary.maximum_block_length == MEASUREMENTS_PER_BIN,
        "the production summary does not contain 32 bins of 100",
    )
    require_condition(
        run.warmup_attempted_updates == expected_warmup_attempts,
        "the warm-up did not attempt exactly 3,200 local updates",
    )
    require_condition(
        run.measurement_attempted_updates == expected_measurement_attempts,
        "measurement sweeps did not attempt exactly 51,200 local updates",
    )
    require_condition(
        run.state.attempted_updates == EXPECTED_ATTEMPTED_UPDATES,
        "the production chain did not attempt exactly 54,400 local updates",
    )
    require_condition(
        0 < run.state.accepted_updates < run.state.attempted_updates,
        "the production chain did not contain both accepted and rejected updates",
    )
    require_condition(
        isfinite(run.acceptance_rate) && 0.0 < run.acceptance_rate < 1.0,
        "the acceptance rate is not finite and strictly between zero and one",
    )
    require_condition(
        summary.configuration_signs_valid,
        "the sampled configuration signs are invalid",
    )
    require_condition(
        summary.raw_observables_finite,
        "a sampled raw observable is non-finite",
    )
    require_condition(
        summary.finite_results,
        "a sampled sign-reweighting result is non-finite",
    )
    require_condition(
        summary.sign_denominator_usable &&
        summary.leave_one_out_denominators_usable &&
        summary.denominator_guard_passed &&
        summary.leave_one_out_denominator_failure_count == 0,
        "a full or leave-one-bin-out denominator failed the existing guard",
    )
    require_condition(
        run.forward_sweep_count == EXPECTED_DIRECTION_SWEEPS &&
        run.reverse_sweep_count == EXPECTED_DIRECTION_SWEEPS &&
        run.forward_sweep_count + run.reverse_sweep_count ==
            EXPECTED_COMPLETE_SWEEPS,
        "both complete-sweep directions were not exercised as expected",
    )
    require_condition(
        run.forward_attempted_updates == EXPECTED_DIRECTION_ATTEMPTS &&
        run.reverse_attempted_updates == EXPECTED_DIRECTION_ATTEMPTS,
        "the two sweep directions did not each attempt 27,200 updates",
    )
    require_condition(
        isfinite(run.final_recomputation.maximum) &&
        run.final_recomputation.maximum <= GREEN_CORRECTION_TOLERANCE,
        "the final Green correction exceeded 1e-9",
    )
    require_condition(
        isfinite(run.state.maximum_greens_recompute_error) &&
        run.state.maximum_greens_recompute_error <=
            GREEN_CORRECTION_TOLERANCE,
        "the cumulative maximum Green correction exceeded 1e-9",
    )
    print_sampling_diagnostics(run)
    println("metropolis_runtime_seconds=", metropolis_runtime_seconds)

    exact_identity_errors = thermodynamic_identity_errors(
        exact.values,
        parameters.mu,
    )
    sampled_identity_errors = thermodynamic_identity_errors(
        summary,
        parameters.mu,
    )
    print_identity_diagnostics(exact_identity_errors, sampled_identity_errors)
    require_condition(
        maximum(exact_identity_errors) <= THERMODYNAMIC_IDENTITY_TOLERANCE,
        "an exact thermodynamic identity failed",
    )
    require_condition(
        maximum(sampled_identity_errors) <=
            THERMODYNAMIC_IDENTITY_TOLERANCE,
        "a sampled thermodynamic identity failed",
    )

    println(
        "average_sign_comparison exact=", exact.average_sign,
        " dqmc=", summary.average_sign,
        " blocked_sign_se=", summary.average_sign_standard_error,
        " absolute_deviation=", abs(summary.average_sign - exact.average_sign),
    )
    validate_average_sign(exact, run)

    rows = comparison_rows(exact.values, summary)
    print_comparison(rows)
    failed_observables = [row.observable for row in rows if !row.passed]
    require_condition(
        isempty(failed_observables),
        "DQMC-exact comparison exceeded max(5*SE, 1e-10) for " *
        join(string.(failed_observables), ", "),
    )

    println("thermodynamic_statistical_comparison_passed=true")
    println("sampled_scope=L2_Ltau2_U2_only")
    println("files_written=false")
    return nothing
end

function main()
    validation_started = time_ns()
    try
        run_validation()
    catch
        println("overall_validation_status=FAIL")
        println(
            "validation_runtime_seconds=",
            (time_ns() - validation_started) / 1.0e9,
        )
        rethrow()
    end
    println("overall_validation_status=PASS")
    println(
        "validation_runtime_seconds=",
        (time_ns() - validation_started) / 1.0e9,
    )
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

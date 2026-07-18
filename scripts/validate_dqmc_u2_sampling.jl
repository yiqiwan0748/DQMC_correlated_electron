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
const BETA = 0.1
const DTAU = 0.1
const INTERACTION = 2.0
const CHEMICAL_POTENTIAL = -1.897
const STABILIZATION_PERIOD = 1
const WARMUP_SWEEPS = 200
const BIN_COUNT = 32
const MEASUREMENTS_PER_BIN = 100
const MEASUREMENT_SWEEPS = BIN_COUNT * MEASUREMENTS_PER_BIN
const SEED = Int(0x5A17C0DE)

function require_condition(condition::Bool, message::AbstractString)
    condition || error(message)
    return nothing
end

function fields_from_mask(mask::Integer, single_particle_dimension::Int)
    fields = Matrix{Int8}(undef, single_particle_dimension, 1)
    for index in 1:single_particle_dimension
        bit = (mask >> (index - 1)) & 1
        fields[index, 1] = iszero(bit) ? Int8(-1) : Int8(1)
    end
    return fields
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

    # The discrete-HS prefactor is positive and identical for all fields. It
    # cancels from normalized ratios, but including it makes log_abs_weights
    # the log of the complete Trotterized auxiliary-field weight.
    log_hs_prefactor = field_count * (
        -log(2.0) - parameters.dtau * parameters.U / 4
    )

    for mask in 0:(configuration_count - 1)
        fields = fields_from_mask(mask, dimensions.Nsp)
        state = initialize_dqmc_state(parameters; hs_fields=fields)

        determinant_sign = Float64(sign(
            inv(state.signdet_green_up) * inv(state.signdet_green_down),
        ))
        measured_sign = measure_basic_observables(state).configuration_sign
        log_abs_fermion_weight =
            -state.logdet_green_up - state.logdet_green_down
        log_abs_weight = log_abs_fermion_weight + log_hs_prefactor
        absolute_weight = exp(log_abs_weight)
        observation = measure_thermodynamic_observables(state)

        finite = determinant_sign in (-1.0, 1.0) &&
            measured_sign == determinant_sign &&
            isfinite(log_abs_weight) &&
            isfinite(absolute_weight) && absolute_weight > 0.0 &&
            all(isfinite, state.greens_up) &&
            all(isfinite, state.greens_down) &&
            all(name -> isfinite(getproperty(observation, name)),
                THERMODYNAMIC_NAMES)
        finite || error(
            "non-finite exact weight, Green matrix, sign, or observable " *
            "for auxiliary-field mask $mask",
        )

        index = mask + 1
        log_abs_weights[index] = log_abs_weight
        configuration_signs[index] = determinant_sign
        raw_observations[index] = observation
        finite_configuration_count += 1
    end

    log_weight_shift = maximum(log_abs_weights)
    scaled_abs_weights = exp.(log_abs_weights .- log_weight_shift)
    absolute_denominator = sum(scaled_abs_weights)
    signed_denominator = sum(configuration_signs .* scaled_abs_weights)
    require_condition(
        isfinite(absolute_denominator) && absolute_denominator > 0.0,
        "the exact absolute-weight denominator is unusable",
    )
    require_condition(
        isfinite(signed_denominator) && !iszero(signed_denominator),
        "the exact signed-weight denominator is unusable",
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
        all(name -> isfinite(getproperty(exact_values, name)),
            THERMODYNAMIC_NAMES),
        "an exact sign-reweighted observable is non-finite",
    )

    return (
        values=exact_values,
        average_sign=signed_denominator / absolute_denominator,
        configuration_count=configuration_count,
        finite_configuration_count=finite_configuration_count,
        positive_sign_count=count(==(1.0), configuration_signs),
        negative_sign_count=count(==(-1.0), configuration_signs),
        minimum_log_abs_weight=minimum(log_abs_weights),
        maximum_log_abs_weight=maximum(log_abs_weights),
        log_weight_shift=log_weight_shift,
        scaled_absolute_denominator=absolute_denominator,
        scaled_signed_denominator=signed_denominator,
    )
end

function run_production_chain(parameters, dimensions)
    # Supplying no fixed field is intentional: initialization and all later
    # shuffles/acceptance tests consume one seeded production RNG trajectory.
    state = initialize_dqmc_state(parameters)
    for _ in 1:WARMUP_SWEEPS
        sweep = sweep_all_time_slices!(state)
        require_condition(
            sweep.attempted == dimensions.Nsp * dimensions.Ltau,
            "a warm-up sweep attempted the wrong number of updates",
        )
    end
    warmup_attempted_updates = state.attempted_updates
    warmup_accepted_updates = state.accepted_updates

    configuration_signs = Vector{Float64}(undef, MEASUREMENT_SWEEPS)
    raw_observations = Vector{NamedTuple}(undef, MEASUREMENT_SWEEPS)
    measurements_accumulated = 0
    for measurement in 1:MEASUREMENT_SWEEPS
        sweep = sweep_all_time_slices!(state)
        require_condition(
            sweep.attempted == dimensions.Nsp * dimensions.Ltau,
            "a measurement sweep attempted the wrong number of updates",
        )
        configuration_signs[measurement] =
            measure_basic_observables(state).configuration_sign
        raw_observations[measurement] =
            measure_thermodynamic_observables(state)
        measurements_accumulated += 1
    end

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
        measurements_accumulated=measurements_accumulated,
        warmup_attempted_updates=warmup_attempted_updates,
        warmup_accepted_updates=warmup_accepted_updates,
        measurement_attempted_updates=
            state.attempted_updates - warmup_attempted_updates,
        measurement_accepted_updates=
            state.accepted_updates - warmup_accepted_updates,
        acceptance_rate=acceptance_rate,
        final_recomputation=final_recomputation,
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
            jackknife_valid && !ismissing(estimate) && !ismissing(standard_error),
            "missing or invalid DQMC estimate/error for $name",
        )
        require_condition(
            isfinite(estimate) && isfinite(standard_error),
            "non-finite DQMC estimate/error for $name",
        )

        absolute_deviation = abs(estimate - exact_value)
        se_normalized_deviation = if iszero(standard_error)
            iszero(absolute_deviation) ? 0.0 : Inf
        else
            absolute_deviation / standard_error
        end
        allowed_deviation = max(5 * standard_error, 1.0e-10)
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

function print_parameters(parameters, dimensions, model_parameters)
    field_count = dimensions.Nsp * dimensions.Ltau
    println("Milestone 6B Phase 1 exact-enumeration sampling validation")
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
    println("warmup_sweeps=", WARMUP_SWEEPS)
    println("contiguous_bin_count=", BIN_COUNT)
    println("measurement_sweeps_per_bin=", MEASUREMENTS_PER_BIN)
    println("measurements_per_sweep=1")
    println("total_measurements=", MEASUREMENT_SWEEPS)
    println(
        "expected_attempted_updates=",
        (WARMUP_SWEEPS + MEASUREMENT_SWEEPS) * field_count,
    )
end

function print_exact_diagnostics(exact)
    println("exact_finite_configuration_count=", exact.finite_configuration_count)
    println("exact_positive_sign_count=", exact.positive_sign_count)
    println("exact_negative_sign_count=", exact.negative_sign_count)
    println("exact_minimum_log_abs_weight=", exact.minimum_log_abs_weight)
    println("exact_maximum_log_abs_weight=", exact.maximum_log_abs_weight)
    println("exact_log_weight_shift=", exact.log_weight_shift)
    println(
        "exact_scaled_absolute_denominator=",
        exact.scaled_absolute_denominator,
    )
    println(
        "exact_scaled_signed_denominator=",
        exact.scaled_signed_denominator,
    )
    println("exact_average_sign=", exact.average_sign)
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
    println("average_sign_standard_error=", summary.average_sign_standard_error)
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
    println(
        "final_green_recomputation_error=",
        run.final_recomputation.maximum,
    )
    println(
        "maximum_greens_recompute_error=",
        state.maximum_greens_recompute_error,
    )
end

function print_comparison(rows)
    println(
        "observable\texact\tdqmc\tjackknife_se\tabs_deviation\t" *
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

function main()
    parameters = DQMCParameters(
        L=L,
        beta=BETA,
        dtau=DTAU,
        U=INTERACTION,
        mu=CHEMICAL_POTENTIAL,
        stabilization_period=STABILIZATION_PERIOD,
        burnin_sweeps=WARMUP_SWEEPS,
        measurement_sweeps=MEASUREMENT_SWEEPS,
        sweeps_between_measurements=1,
        seed=SEED,
    )
    dimensions = validate_dqmc_parameters(parameters)
    model_parameters = DEFAULT_PARAMETERS
    require_condition(parameters.L == 2, "this validation is restricted to L=2")
    require_condition(dimensions.Ltau == 1,
                      "this validation is restricted to L_tau=1")
    require_condition(parameters.U == 2.0,
                      "this validation is restricted to U/t=2")
    require_condition(parameters.stabilization_period == 1,
                      "L_tau=1 requires stabilization_period=1")
    require_condition(
        (model_parameters.t1, model_parameters.t2,
         model_parameters.t3, model_parameters.t4) ==
        (-1.0, -1.75, -0.85, -0.65),
        "the hopping parameters differ from the authorized model",
    )

    print_parameters(parameters, dimensions, model_parameters)

    exact = exact_auxiliary_field_average(parameters, dimensions)
    require_condition(exact.configuration_count == 256,
                      "exact enumeration did not contain 256 configurations")
    require_condition(exact.finite_configuration_count == 256,
                      "not all 256 exact configurations were finite")
    print_exact_diagnostics(exact)

    run = run_production_chain(parameters, dimensions)
    expected_attempts = (WARMUP_SWEEPS + MEASUREMENT_SWEEPS) *
                        dimensions.Nsp * dimensions.Ltau
    require_condition(run.measurements_accumulated == MEASUREMENT_SWEEPS,
                      "not all 3,200 measurements were accumulated")
    require_condition(run.state.attempted_updates == expected_attempts,
                      "the production chain did not attempt 27,200 updates")
    require_condition(run.state.accepted_updates > 0,
                      "the production chain accepted no updates")
    require_condition(
        run.state.accepted_updates < run.state.attempted_updates,
        "the production chain rejected no updates",
    )
    print_sampling_diagnostics(run)

    summary = run.summary
    require_condition(summary.measurement_count == MEASUREMENT_SWEEPS,
                      "the reweighting summary has the wrong measurement count")
    require_condition(summary.block_count == BIN_COUNT,
                      "the reweighting summary has the wrong block count")
    require_condition(summary.configuration_signs_valid,
                      "the production configuration signs are invalid")
    require_condition(summary.raw_observables_finite,
                      "a production raw observable is non-finite")
    require_condition(summary.finite_results,
                      "a production reweighting result is non-finite")
    require_condition(summary.denominator_guard_passed,
                      "a full or leave-one-out denominator guard failed")
    require_condition(
        !ismissing(summary.average_sign_standard_error) &&
        isfinite(summary.average_sign) &&
        isfinite(summary.average_sign_standard_error),
        "the average sign or its uncertainty is non-finite",
    )

    exact_identity_errors = thermodynamic_identity_errors(
        exact.values,
        parameters.mu,
    )
    dqmc_identity_errors = thermodynamic_identity_errors(summary, parameters.mu)
    maximum_exact_identity_error = maximum(exact_identity_errors)
    maximum_dqmc_identity_error = maximum(dqmc_identity_errors)
    println("maximum_exact_thermodynamic_identity_error=",
            maximum_exact_identity_error)
    println("maximum_dqmc_thermodynamic_identity_error=",
            maximum_dqmc_identity_error)
    require_condition(maximum_exact_identity_error <= 1.0e-10,
                      "an exact thermodynamic identity failed")
    require_condition(maximum_dqmc_identity_error <= 1.0e-10,
                      "a DQMC thermodynamic identity failed")

    println(
        "average_sign_comparison exact=", exact.average_sign,
        " dqmc=", summary.average_sign,
        " jackknife_se=", summary.average_sign_standard_error,
        " absolute_deviation=", abs(summary.average_sign - exact.average_sign),
    )
    rows = comparison_rows(exact.values, summary)
    print_comparison(rows)
    failed_observables = [row.observable for row in rows if !row.passed]
    if !isempty(failed_observables)
        println("statistical_comparison_failed=", join(failed_observables, ','))
        error(
            "DQMC-exact comparison exceeded max(5*SE, 1e-10) for " *
            join(string.(failed_observables), ", "),
        )
    end

    println("statistical_comparison_passed=true")
    println("sampled_scope=L2_Ltau1_U2_only")
    println("files_written=false")
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

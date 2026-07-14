using LinearAlgebra: I
using Random: AbstractRNG, Xoshiro, rand, shuffle!

import JDQMCFramework as jdqmcf

"""Parameters for the correctness-first two-orbital DQMC kernel."""
struct DQMCParameters{T<:AbstractFloat}
    L::Int
    beta::T
    dtau::T
    U::T
    mu::T
    stabilization_period::Int
    burnin_sweeps::Int
    measurement_sweeps::Int
    sweeps_between_measurements::Int
    seed::Int
end

"""
    DQMCParameters(; L=2, beta=1.0, dtau=0.1, U=2.0, mu=0.0,
                     stabilization_period=5, burnin_sweeps=20,
                     measurement_sweeps=50,
                     sweeps_between_measurements=1, seed=1234)

Construct and validate simulation parameters. The four real inputs are
promoted to one floating-point type; ordinary inputs use `Float64`.
"""
function DQMCParameters(;
    L::Integer=2,
    beta::Real=1.0,
    dtau::Real=0.1,
    U::Real=2.0,
    mu::Real=0.0,
    stabilization_period::Integer=5,
    burnin_sweeps::Integer=20,
    measurement_sweeps::Integer=50,
    sweeps_between_measurements::Integer=1,
    seed::Integer=1234,
)
    beta_value, dtau_value, interaction, chemical_potential =
        promote(float(beta), float(dtau), float(U), float(mu))
    T = typeof(beta_value)
    parameters = DQMCParameters{T}(
        Int(L), beta_value, dtau_value, interaction, chemical_potential,
        Int(stabilization_period), Int(burnin_sweeps),
        Int(measurement_sweeps), Int(sweeps_between_measurements), Int(seed),
    )
    validate_dqmc_parameters(parameters)
    return parameters
end

"""
    validate_dqmc_parameters(parameters)

Validate a parameter set and return `(Nsites, Nsp, Ltau)`. Integrality of
`beta/dtau` uses zero relative tolerance and absolute tolerance
`64*eps(T)*max(1, abs(beta/dtau))`. This covers representation roundoff while
rejecting a physically different imaginary-time grid.
"""
function validate_dqmc_parameters(parameters::DQMCParameters{T}) where {T}
    (; L, beta, dtau, U, mu, stabilization_period, burnin_sweeps,
       measurement_sweeps, sweeps_between_measurements) = parameters
    _, number_of_sites, single_particle_dimension = _lattice_dimensions(L)
    isfinite(beta) && beta > zero(T) ||
        throw(ArgumentError("beta must be finite and positive"))
    isfinite(dtau) && dtau > zero(T) ||
        throw(ArgumentError("dtau must be finite and positive"))
    isfinite(U) && U >= zero(T) ||
        throw(ArgumentError("U must be finite and nonnegative"))
    isfinite(mu) || throw(ArgumentError("mu must be finite"))
    stabilization_period > 0 ||
        throw(ArgumentError("stabilization_period must be positive"))
    burnin_sweeps >= 0 ||
        throw(ArgumentError("burnin_sweeps must be nonnegative"))
    measurement_sweeps > 0 ||
        throw(ArgumentError("measurement_sweeps must be positive"))
    sweeps_between_measurements > 0 || throw(ArgumentError(
        "sweeps_between_measurements must be positive",
    ))

    ratio = beta / dtau
    isfinite(ratio) || throw(ArgumentError("beta/dtau must be finite"))
    Ltau = try
        round(Int, ratio)
    catch error
        error isa InexactError || error isa OverflowError || rethrow()
        throw(ArgumentError("beta/dtau must fit in an Int"))
    end
    Ltau > 0 || throw(ArgumentError("beta/dtau must be positive"))
    tolerance = 64 * eps(T) * max(one(T), abs(ratio))
    abs(ratio - convert(T, Ltau)) <= tolerance || throw(ArgumentError(
        "beta/dtau must be integral within tolerance $tolerance",
    ))
    return (
        Nsites=number_of_sites,
        Nsp=single_particle_dimension,
        Ltau=Ltau,
    )
end

"""Return `alpha = acosh(exp(dtau*U/2))` for the real repulsive HS field."""
function hubbard_stratonovich_alpha(U::Real, dtau::Real)
    interaction, time_step = promote(float(U), float(dtau))
    T = typeof(interaction)
    isfinite(interaction) && interaction >= zero(T) ||
        throw(ArgumentError("U must be finite and nonnegative"))
    isfinite(time_step) && time_step > zero(T) ||
        throw(ArgumentError("dtau must be finite and positive"))
    alpha = acosh(exp(time_step * interaction / 2))
    isfinite(alpha) || throw(OverflowError("the HS coupling is not finite"))
    return alpha
end

hubbard_stratonovich_alpha(parameters::DQMCParameters) =
    hubbard_stratonovich_alpha(parameters.U, parameters.dtau)

"""
    initialize_hs_fields(parameters[, rng])

Return an `Nsp x Ltau` `Int8` matrix of independent `+1` and `-1` fields.
The keyword form initializes `Xoshiro(parameters.seed)`.
"""
function initialize_hs_fields(parameters::DQMCParameters, rng::AbstractRNG)
    (; Nsp, Ltau) = validate_dqmc_parameters(parameters)
    fields = Matrix{Int8}(undef, Nsp, Ltau)
    for index in eachindex(fields)
        fields[index] = rand(rng, Bool) ? Int8(1) : Int8(-1)
    end
    return fields
end

initialize_hs_fields(
    parameters::DQMCParameters;
    rng::AbstractRNG=Xoshiro(parameters.seed),
) = initialize_hs_fields(parameters, rng)

function _checked_hs_fields(fields::AbstractMatrix, Nsp::Int, Ltau::Int)
    size(fields) == (Nsp, Ltau) || throw(DimensionMismatch(
        "HS fields must have size (Nsp, Ltau) = ($Nsp, $Ltau)",
    ))
    all(field -> field == -1 || field == 1, fields) ||
        throw(ArgumentError("every HS field must be +1 or -1"))
    return Matrix{Int8}(fields)
end

"""
    build_slice_propagators(parameters, K_mu, hs_fields)

Build dense asymmetric propagators
`B_up[l] = diag(exp(+alpha*s[:,l]))*exp(-dtau*K_mu)` and
`B_down[l] = diag(exp(-alpha*s[:,l]))*exp(-dtau*K_mu)`.
"""
function build_slice_propagators(
    parameters::DQMCParameters{T},
    kinetic_mu::AbstractMatrix{<:Real},
    hs_fields::AbstractMatrix,
) where {T}
    (; Nsp, Ltau) = validate_dqmc_parameters(parameters)
    size(kinetic_mu) == (Nsp, Nsp) || throw(DimensionMismatch(
        "K_mu must have size (Nsp, Nsp) = ($Nsp, $Nsp)",
    ))
    kinetic_mu_typed = Matrix{T}(kinetic_mu)
    _checked_real_symmetric_matrix(kinetic_mu_typed, "K_mu")
    fields = _checked_hs_fields(hs_fields, Nsp, Ltau)
    expm_dtau_kinetic, expp_dtau_kinetic =
        dense_slice_propagators(kinetic_mu_typed, parameters.dtau)
    alpha = convert(T, hubbard_stratonovich_alpha(parameters))

    propagators_up = jdqmcf.AsymExactPropagator{T,T}[]
    propagators_down = jdqmcf.AsymExactPropagator{T,T}[]
    sizehint!(propagators_up, Ltau)
    sizehint!(propagators_down, Ltau)
    for time_slice in 1:Ltau
        potential_up = Vector{T}(undef, Nsp)
        potential_down = Vector{T}(undef, Nsp)
        for index in 1:Nsp
            exponent = alpha * fields[index, time_slice]
            potential_up[index] = exp(exponent)
            potential_down[index] = exp(-exponent)
        end
        push!(propagators_up, jdqmcf.AsymExactPropagator(
            potential_up, expm_dtau_kinetic, expp_dtau_kinetic,
        ))
        push!(propagators_down, jdqmcf.AsymExactPropagator(
            potential_down, expm_dtau_kinetic, expp_dtau_kinetic,
        ))
    end
    return (
        up=propagators_up,
        down=propagators_down,
        expm_dtau_kinetic=expm_dtau_kinetic,
        expp_dtau_kinetic=expp_dtau_kinetic,
        alpha=alpha,
    )
end

"""
Mutable state for two spin determinants. Stored log determinants and signs
describe `det(G_sigma)`; the sampled weight contains `det(G_sigma^-1)`.
"""
mutable struct DQMCState{
    T<:AbstractFloat,
    R<:AbstractRNG,
    P<:jdqmcf.AbstractPropagator{T,T},
    C<:jdqmcf.FermionGreensCalculator{T,T},
}
    parameters::DQMCParameters{T}
    kinetic::Matrix{T}
    kinetic_mu::Matrix{T}
    alpha::T
    hs_fields::Matrix{Int8}
    propagators_up::Vector{P}
    propagators_down::Vector{P}
    greens_calculator_up::C
    greens_calculator_down::C
    greens_up::Matrix{T}
    greens_down::Matrix{T}
    logdet_green_up::T
    signdet_green_up::T
    logdet_green_down::T
    signdet_green_down::T
    rng::R
    update_u::Vector{T}
    update_v::Vector{T}
    wrap_workspace::Matrix{T}
    permutation::Vector{Int}
    attempted_updates::Int
    accepted_updates::Int
    maximum_greens_recompute_error::T
end

"""
    initialize_dqmc_state(parameters; model_parameters=DEFAULT_PARAMETERS,
                          hs_fields=nothing, rng=nothing)

Construct `K`, `K_mu`, the HS fields and propagators, then calculate both
equal-time Green matrices with JDQMCFramework's stabilized routine.
"""
function initialize_dqmc_state(
    parameters::DQMCParameters{T};
    model_parameters::ModelParameters=DEFAULT_PARAMETERS,
    hs_fields=nothing,
    rng::Union{Nothing,AbstractRNG}=nothing,
) where {T}
    (; Nsp, Ltau) = validate_dqmc_parameters(parameters)
    simulation_rng = isnothing(rng) ? Xoshiro(parameters.seed) : rng
    kinetic = Matrix{T}(realspace_kinetic_matrix(parameters.L, model_parameters))
    kinetic_mu = Matrix{T}(shifted_kinetic_matrix(kinetic, parameters.mu))
    fields = isnothing(hs_fields) ?
             initialize_hs_fields(parameters, simulation_rng) :
             _checked_hs_fields(hs_fields, Nsp, Ltau)
    slices = build_slice_propagators(parameters, kinetic_mu, fields)
    calculator_up = jdqmcf.FermionGreensCalculator(
        slices.up, parameters.beta, parameters.dtau,
        parameters.stabilization_period,
    )
    calculator_down = jdqmcf.FermionGreensCalculator(
        slices.down, parameters.beta, parameters.dtau,
        parameters.stabilization_period,
    )
    greens_up = zeros(T, Nsp, Nsp)
    greens_down = zeros(T, Nsp, Nsp)
    logdet_green_up, signdet_green_up =
        jdqmcf.calculate_equaltime_greens!(greens_up, calculator_up)
    logdet_green_down, signdet_green_down =
        jdqmcf.calculate_equaltime_greens!(greens_down, calculator_down)
    return DQMCState(
        parameters, kinetic, kinetic_mu, slices.alpha, fields,
        slices.up, slices.down, calculator_up, calculator_down,
        greens_up, greens_down, logdet_green_up, signdet_green_up,
        logdet_green_down, signdet_green_down, simulation_rng,
        zeros(T, Nsp), zeros(T, Nsp), zeros(T, Nsp, Nsp),
        collect(1:Nsp), 0, 0, zero(T),
    )
end

function _checked_update_location(state::DQMCState, index::Integer, time_slice::Integer)
    (; Nsp, Ltau) = validate_dqmc_parameters(state.parameters)
    1 <= index <= Nsp || throw(BoundsError(state.hs_fields, (index, time_slice)))
    1 <= time_slice <= Ltau ||
        throw(BoundsError(state.hs_fields, (index, time_slice)))
    return Int(index), Int(time_slice)
end

"""
    local_hs_flip_ratio(state, index, time_slice)

Return spin-resolved and total determinant ratios. The Green matrices must
correspond to the supplied cyclic slice, as inside `sweep_time_slice!`.
"""
function local_hs_flip_ratio(
    state::DQMCState{T}, index::Integer, time_slice::Integer,
) where {T}
    i, l = _checked_update_location(state, index, time_slice)
    field = state.hs_fields[i, l]
    delta_v_up = 2 * state.alpha * field / state.parameters.dtau
    delta_v_down = -delta_v_up
    ratio_up_raw, delta_up = jdqmcf.local_update_det_ratio(
        state.greens_up, delta_v_up, i, state.parameters.dtau,
    )
    ratio_down_raw, delta_down = jdqmcf.local_update_det_ratio(
        state.greens_down, delta_v_down, i, state.parameters.dtau,
    )
    ratio_up = convert(T, ratio_up_raw)
    ratio_down = convert(T, ratio_down_raw)
    return (
        up=ratio_up,
        down=ratio_down,
        total=ratio_up * ratio_down,
        delta_up=convert(T, delta_up),
        delta_down=convert(T, delta_down),
    )
end

"""
    attempt_local_hs_flip!(state, index, time_slice;
                           acceptance_uniform=nothing)

Attempt a Metropolis update with probability `min(1, abs(R_up*R_down))`.
`acceptance_uniform` in `[0,1)` is a deterministic validation hook.
"""
function attempt_local_hs_flip!(
    state::DQMCState{T}, index::Integer, time_slice::Integer;
    acceptance_uniform::Union{Nothing,Real}=nothing,
) where {T}
    i, l = _checked_update_location(state, index, time_slice)
    ratios = local_hs_flip_ratio(state, i, l)
    probability = min(one(T), abs(ratios.total))
    uniform = if isnothing(acceptance_uniform)
        rand(state.rng)
    else
        value = float(acceptance_uniform)
        isfinite(value) && 0 <= value < 1 || throw(ArgumentError(
            "acceptance_uniform must be finite and lie in [0, 1)",
        ))
        value
    end
    state.attempted_updates += 1
    accepted = uniform < probability
    if accepted
        state.accepted_updates += 1
        state.hs_fields[i, l] = -state.hs_fields[i, l]
        state.logdet_green_up, state.signdet_green_up =
            jdqmcf.local_update_greens!(
                state.greens_up, state.logdet_green_up,
                state.signdet_green_up, state.propagators_up[l], ratios.up,
                ratios.delta_up, i, state.update_u, state.update_v,
            )
        state.logdet_green_down, state.signdet_green_down =
            jdqmcf.local_update_greens!(
                state.greens_down, state.logdet_green_down,
                state.signdet_green_down, state.propagators_down[l],
                ratios.down, ratios.delta_down, i,
                state.update_u, state.update_v,
            )
    end
    return (
        accepted=accepted,
        probability=probability,
        ratio_up=ratios.up,
        ratio_down=ratios.down,
        total_ratio=ratios.total,
    )
end

function _check_calculator_synchronization(state::DQMCState, time_slice::Int)
    up = state.greens_calculator_up
    down = state.greens_calculator_down
    up.l == time_slice && down.l == time_slice || throw(ArgumentError(
        "Green calculators are not positioned at time slice $time_slice",
    ))
    up.forward == down.forward || throw(ArgumentError(
        "spin-sector Green calculators are not synchronized",
    ))
    return nothing
end

"""
    sweep_time_slice!(state, time_slice)

Propagate to a slice yielded by the up-spin calculator, visit every
site-orbital in seeded random order, and invoke framework stabilization.
"""
function sweep_time_slice!(state::DQMCState{T}, time_slice::Integer) where {T}
    l = Int(time_slice)
    _check_calculator_synchronization(state, l)
    attempts_before = state.attempted_updates
    accepted_before = state.accepted_updates
    jdqmcf.propagate_equaltime_greens!(
        state.greens_up, state.greens_calculator_up, state.propagators_up,
    )
    jdqmcf.propagate_equaltime_greens!(
        state.greens_down, state.greens_calculator_down,
        state.propagators_down,
    )
    jdqmcf.partially_wrap_greens_reverse!(
        state.greens_up, state.propagators_up[l], state.wrap_workspace,
    )
    jdqmcf.partially_wrap_greens_reverse!(
        state.greens_down, state.propagators_down[l], state.wrap_workspace,
    )
    shuffle!(state.rng, state.permutation)
    for index in state.permutation
        attempt_local_hs_flip!(state, index, l)
    end
    jdqmcf.partially_wrap_greens_forward!(
        state.greens_up, state.propagators_up[l], state.wrap_workspace,
    )
    jdqmcf.partially_wrap_greens_forward!(
        state.greens_down, state.propagators_down[l], state.wrap_workspace,
    )
    state.logdet_green_up, state.signdet_green_up, error_up, _ =
        jdqmcf.stabilize_equaltime_greens!(
            state.greens_up, state.logdet_green_up,
            state.signdet_green_up, state.greens_calculator_up,
            state.propagators_up; update_B̄=true,
        )
    state.logdet_green_down, state.signdet_green_down, error_down, _ =
        jdqmcf.stabilize_equaltime_greens!(
            state.greens_down, state.logdet_green_down,
            state.signdet_green_down, state.greens_calculator_down,
            state.propagators_down; update_B̄=true,
        )
    corrected_error = max(error_up, error_down)
    state.maximum_greens_recompute_error = max(
        state.maximum_greens_recompute_error, corrected_error,
    )
    # The up calculator advances through its surrounding `for` loop. Match
    # that pending advance for the down calculator, as in the frozen reference.
    iterate(state.greens_calculator_down, state.greens_calculator_up)
    attempted = state.attempted_updates - attempts_before
    accepted = state.accepted_updates - accepted_before
    return (
        attempted=attempted,
        accepted=accepted,
        acceptance_rate=accepted / attempted,
        recompute_error=corrected_error,
    )
end

"""Perform one complete alternating-direction sweep over all time slices."""
function sweep_all_time_slices!(state::DQMCState)
    attempts_before = state.attempted_updates
    accepted_before = state.accepted_updates
    largest_error = zero(eltype(state.greens_up))
    for time_slice in state.greens_calculator_up
        slice_result = sweep_time_slice!(state, time_slice)
        largest_error = max(largest_error, slice_result.recompute_error)
    end
    attempted = state.attempted_updates - attempts_before
    accepted = state.accepted_updates - accepted_before
    return (
        attempted=attempted,
        accepted=accepted,
        acceptance_rate=accepted / attempted,
        recompute_error=largest_error,
    )
end

"""
    recompute_equal_time_greens!(state)

Recalculate both `G(0,0)` matrices and all framework factorizations from the
slice propagators at a complete-sweep boundary. Return elementwise errors.
"""
function recompute_equal_time_greens!(state::DQMCState)
    up = state.greens_calculator_up
    down = state.greens_calculator_down
    boundary_up = up.l == 1 || up.l == up.Lτ
    boundary_down = down.l == 1 || down.l == down.Lτ
    boundary_up && boundary_down || throw(ArgumentError(
        "full Green recomputation requires a complete-sweep boundary",
    ))
    up.l == down.l && up.forward == down.forward || throw(ArgumentError(
        "spin-sector Green calculators are not synchronized",
    ))
    old_up = copy(state.greens_up)
    old_down = copy(state.greens_down)
    state.logdet_green_up, state.signdet_green_up =
        jdqmcf.calculate_equaltime_greens!(
            state.greens_up, up, state.propagators_up,
        )
    state.logdet_green_down, state.signdet_green_down =
        jdqmcf.calculate_equaltime_greens!(
            state.greens_down, down, state.propagators_down,
        )
    error_up = maximum(abs.(old_up .- state.greens_up))
    error_down = maximum(abs.(old_down .- state.greens_down))
    corrected_error = max(error_up, error_down)
    state.maximum_greens_recompute_error = max(
        state.maximum_greens_recompute_error, corrected_error,
    )
    return (up=error_up, down=error_down, maximum=corrected_error)
end

function _configuration_sign(state::DQMCState{T}) where {T}
    return convert(T, sign(
        inv(state.signdet_green_up) * inv(state.signdet_green_down),
    ))
end

"""
    measure_basic_observables(state)

Measure the current sign, total/orbital fillings per physical site, and
orbital double occupancies using `n_sigma,i = 1-G_sigma[i,i]`.
"""
function measure_basic_observables(state::DQMCState{T}) where {T}
    number_of_sites = state.parameters.L^2
    density_x = zero(T)
    density_y = zero(T)
    double_x = zero(T)
    double_y = zero(T)
    for site in 0:(number_of_sites - 1)
        index_x = 2 * site + 1
        index_y = index_x + 1
        occupation_up_x = one(T) - state.greens_up[index_x, index_x]
        occupation_down_x = one(T) - state.greens_down[index_x, index_x]
        occupation_up_y = one(T) - state.greens_up[index_y, index_y]
        occupation_down_y = one(T) - state.greens_down[index_y, index_y]
        density_x += occupation_up_x + occupation_down_x
        density_y += occupation_up_y + occupation_down_y
        double_x += occupation_up_x * occupation_down_x
        double_y += occupation_up_y * occupation_down_y
    end
    density_x /= number_of_sites
    density_y /= number_of_sites
    double_x /= number_of_sites
    double_y /= number_of_sites
    return (
        configuration_sign=_configuration_sign(state),
        total_density=density_x + density_y,
        orbital_x_density=density_x,
        orbital_y_density=density_y,
        double_occupancy_x=double_x,
        double_occupancy_y=double_y,
    )
end

"""
    run_dqmc_smoke(parameters; case_name="smoke")

Run the small burn-in and measurement workload. Interacting observables are
sign reweighted; acceptance includes burn-in and measurement sweeps.
"""
function run_dqmc_smoke(
    parameters::DQMCParameters{T};
    case_name::AbstractString="smoke",
    model_parameters::ModelParameters=DEFAULT_PARAMETERS,
) where {T}
    dimensions = validate_dqmc_parameters(parameters)
    state = initialize_dqmc_state(
        parameters; model_parameters=model_parameters,
    )
    start_time = time_ns()
    for _ in 1:parameters.burnin_sweeps
        sweep_all_time_slices!(state)
    end

    sign_sum = zero(T)
    total_density_sum = zero(T)
    density_x_sum = zero(T)
    density_y_sum = zero(T)
    double_x_sum = zero(T)
    double_y_sum = zero(T)
    for _ in 1:parameters.measurement_sweeps
        for _ in 1:parameters.sweeps_between_measurements
            sweep_all_time_slices!(state)
        end
        observables = measure_basic_observables(state)
        configuration_sign = observables.configuration_sign
        sign_sum += configuration_sign
        total_density_sum += configuration_sign * observables.total_density
        density_x_sum += configuration_sign * observables.orbital_x_density
        density_y_sum += configuration_sign * observables.orbital_y_density
        double_x_sum += configuration_sign * observables.double_occupancy_x
        double_y_sum += configuration_sign * observables.double_occupancy_y
    end
    iszero(sign_sum) && throw(ErrorException(
        "the measured sign sum is zero; sign-reweighted observables are undefined",
    ))
    recompute_equal_time_greens!(state)
    runtime_seconds = (time_ns() - start_time) / 1.0e9
    acceptance_rate = state.accepted_updates / state.attempted_updates
    return (
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
        acceptance_rate=acceptance_rate,
        average_sign=sign_sum / parameters.measurement_sweeps,
        total_density=total_density_sum / sign_sum,
        orbital_x_density=density_x_sum / sign_sum,
        orbital_y_density=density_y_sum / sign_sum,
        double_occupancy_x=double_x_sum / sign_sum,
        double_occupancy_y=double_y_sum / sign_sum,
        maximum_greens_recompute_error=
            state.maximum_greens_recompute_error,
        runtime_seconds=runtime_seconds,
    )
end

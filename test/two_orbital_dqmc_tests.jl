import JDQMCFramework as jdqmcf

const DQMC_U0_GREEN_MAX = Ref(0.0)
const DQMC_DIRECT_GREEN_MAX = Ref(0.0)
const DQMC_RATIO_MAX = Ref(0.0)
const DQMC_ACCEPTED_UPDATE_MAX = Ref(0.0)
const DQMC_SWEEP_RECOMPUTE_MAX = Ref(0.0)

# This tolerance covers stabilized dense propagation and one-body trace
# accumulation while remaining negligible on the order-one energy scale.
const DQMC_THERMODYNAMIC_ATOL = 1.0e-10

function _dense_dqmc_propagator(propagator)
    return Diagonal(propagator.expmΔτV) * propagator.expmΔτK
end

function _direct_cyclic_product(propagators, time_slice::Integer)
    Ltau = length(propagators)
    1 <= time_slice <= Ltau || throw(BoundsError(propagators, time_slice))
    dimension = size(first(propagators), 1)
    product = Matrix{eltype(first(propagators).expmΔτK)}(
        I,
        dimension,
        dimension,
    )
    for slice in time_slice:-1:1
        product = product * _dense_dqmc_propagator(propagators[slice])
    end
    for slice in Ltau:-1:(time_slice + 1)
        product = product * _dense_dqmc_propagator(propagators[slice])
    end
    return product
end

function _direct_cyclic_green(propagators, time_slice::Integer)
    product = _direct_cyclic_product(propagators, time_slice)
    dimension = size(product, 1)
    return inv(Matrix{eltype(product)}(I, dimension, dimension) + product)
end

function _direct_fermion_determinant(propagators, time_slice::Integer)
    product = _direct_cyclic_product(propagators, time_slice)
    dimension = size(product, 1)
    return det(Matrix{eltype(product)}(I, dimension, dimension) + product)
end

function _set_state_greens_to_direct_slice!(state, time_slice)
    direct_up = _direct_cyclic_green(state.propagators_up, time_slice)
    direct_down = _direct_cyclic_green(state.propagators_down, time_slice)
    copyto!(state.greens_up, direct_up)
    copyto!(state.greens_down, direct_down)
    state.logdet_green_up, state.signdet_green_up = logabsdet(direct_up)
    state.logdet_green_down, state.signdet_green_down = logabsdet(direct_down)
    return nothing
end

function _without_runtime(summary)
    names = filter(!=(:runtime_seconds), propertynames(summary))
    return NamedTuple{Tuple(names)}(Tuple(getproperty(summary, name) for name in names))
end

function _test_smoke_invariants(summary)
    @test 0 <= summary.acceptance_rate <= 1
    @test -1 <= summary.average_sign <= 1
    @test 0 <= summary.total_density <= 4
    @test 0 <= summary.orbital_x_density <= 2
    @test 0 <= summary.orbital_y_density <= 2
    @test isapprox(
        summary.total_density,
        summary.orbital_x_density + summary.orbital_y_density;
        atol=2.0e-12,
        rtol=0.0,
    )
    @test 0 <= summary.double_occupancy_x <= 1
    @test 0 <= summary.double_occupancy_y <= 1
    @test summary.maximum_greens_recompute_error >= 0
    @test summary.runtime_seconds >= 0
    @test all(
        isfinite,
        (
            summary.acceptance_rate,
            summary.average_sign,
            summary.total_density,
            summary.orbital_x_density,
            summary.orbital_y_density,
            summary.double_occupancy_x,
            summary.double_occupancy_y,
            summary.maximum_greens_recompute_error,
            summary.runtime_seconds,
        ),
    )
    return nothing
end

@testset "Milestone 5 two-orbital DQMC" begin
    @testset "Parameters, dimensions, and fields" begin
        for L in (2, 3, 4)
            parameters = DQMCParameters(
                L=L,
                beta=0.2,
                dtau=0.1,
                U=1.0,
                stabilization_period=1,
                burnin_sweeps=0,
                measurement_sweeps=1,
                seed=5000 + L,
            )
            dimensions = validate_dqmc_parameters(parameters)
            fields = initialize_hs_fields(parameters)
            state = initialize_dqmc_state(parameters; hs_fields=fields)
            @test dimensions.Nsites == L^2
            @test dimensions.Nsp == 2 * L^2
            @test dimensions.Ltau == 2
            @test size(fields) == (2 * L^2, 2)
            @test all(field -> field == -1 || field == 1, fields)
            @test size(state.greens_up) == (2 * L^2, 2 * L^2)
            @test size(state.greens_down) == (2 * L^2, 2 * L^2)
            @test length(state.propagators_up) == 2
            @test length(state.propagators_down) == 2
        end

        first_fields = initialize_hs_fields(DQMCParameters(seed=9041))
        second_fields = initialize_hs_fields(DQMCParameters(seed=9041))
        @test first_fields == second_fields
        @test_throws ArgumentError DQMCParameters(beta=1.0, dtau=0.3)
        @test_throws ArgumentError DQMCParameters(U=-1.0)
        @test_throws ArgumentError DQMCParameters(stabilization_period=0)
    end

    @testset "Scalar discrete HS identity" begin
        U = 2.3
        dtau = 0.1
        alpha = hubbard_stratonovich_alpha(U, dtau)
        @test alpha == acosh(exp(dtau * U / 2))
        for occupation_up in (0, 1), occupation_down in (0, 1)
            lhs = exp(
                -dtau * U *
                (occupation_up - 0.5) * (occupation_down - 0.5),
            )
            rhs = 0.5 * exp(-dtau * U / 4) * sum(
                exp(
                    alpha * field *
                    (occupation_up - occupation_down),
                ) for field in (-1, 1)
            )
            @test isapprox(lhs, rhs; atol=8eps(Float64), rtol=8eps(Float64))
        end
        @test hubbard_stratonovich_alpha(0.0, dtau) == 0.0
    end

    @testset "U=0 integration" begin
        cases = (
            (L=2, beta=1.0, dtau=0.1, mu=0.0),
            (L=3, beta=0.6, dtau=0.1, mu=-1.897),
            (L=4, beta=0.4, dtau=0.1, mu=0.35),
        )
        for (case_index, case) in enumerate(cases)
            parameters = DQMCParameters(
                L=case.L,
                beta=case.beta,
                dtau=case.dtau,
                U=0.0,
                mu=case.mu,
                stabilization_period=2,
                burnin_sweeps=0,
                measurement_sweeps=1,
                seed=6100 + case_index,
            )
            state = initialize_dqmc_state(parameters)
            expected_green = free_equal_time_greens(
                state.kinetic_mu,
                parameters.beta,
            )
            mismatch = max(
                maximum(abs.(state.greens_up .- expected_green)),
                maximum(abs.(state.greens_down .- expected_green)),
            )
            DQMC_U0_GREEN_MAX[] = max(DQMC_U0_GREEN_MAX[], mismatch)
            @test mismatch < 1.0e-11

            expected_density = Matrix{Float64}(I, size(expected_green)...) .-
                               expected_green
            expected_orbitals = orbital_resolved_density(
                expected_density,
                case.L,
            )
            expected_total = total_filling_per_physical_site(
                expected_density,
                case.L,
            )
            measured = measure_basic_observables(state)
            @test measured.configuration_sign == 1.0
            @test isapprox(
                measured.total_density,
                expected_total;
                atol=1.0e-11,
                rtol=0.0,
            )
            @test isapprox(
                measured.orbital_x_density,
                expected_orbitals.n_x;
                atol=1.0e-11,
                rtol=0.0,
            )
            @test isapprox(
                measured.orbital_y_density,
                expected_orbitals.n_y;
                atol=1.0e-11,
                rtol=0.0,
            )
        end
    end

    @testset "Deterministic U=0 thermodynamic benchmark" begin
        L = 4
        beta = 4.0
        mu = -1.897
        dtau = 0.1
        parameters = DQMCParameters(
            L=L,
            beta=beta,
            dtau=dtau,
            U=0.0,
            mu=mu,
            stabilization_period=5,
            burnin_sweeps=0,
            measurement_sweeps=1,
            seed=6197,
        )
        dimensions = validate_dqmc_parameters(parameters)
        fixed_fields = [
            isodd(index + time_slice) ? Int8(1) : Int8(-1)
            for index in 1:dimensions.Nsp, time_slice in 1:dimensions.Ltau
        ]
        state = initialize_dqmc_state(parameters; hs_fields=fixed_fields)

        exact = finite_lattice_free_thermodynamics(L, beta, mu)
        measured = measure_thermodynamic_observables(state)
        @test measure_basic_observables(state).configuration_sign == 1.0
        @test propertynames(measured) == propertynames(exact)
        for name in propertynames(exact)
            @test isapprox(
                getproperty(measured, name),
                getproperty(exact, name);
                atol=DQMC_THERMODYNAMIC_ATOL,
                rtol=0.0,
            )
        end

        @test measured.shifted_interaction_energy_per_cell == 0.0
        @test measured.internal_energy_per_cell ==
              measured.kinetic_energy_per_cell
        @test measured.grand_energy_per_cell ==
              measured.internal_energy_per_cell - mu * measured.density

        dimension = dimensions.Nsp
        identity_matrix = Matrix{Float64}(I, dimension, dimension)
        density_up = identity_matrix - state.greens_up
        density_down = identity_matrix - state.greens_down
        bare_kinetic_energy = (
            tr(state.kinetic * density_up) +
            tr(state.kinetic * density_down)
        ) / L^2
        shifted_kinetic_energy = (
            tr(state.kinetic_mu * density_up) +
            tr(state.kinetic_mu * density_down)
        ) / L^2
        @test isapprox(
            measured.internal_energy_per_cell,
            bare_kinetic_energy;
            atol=DQMC_THERMODYNAMIC_ATOL,
            rtol=0.0,
        )
        @test isapprox(
            measured.grand_energy_per_cell,
            shifted_kinetic_energy;
            atol=DQMC_THERMODYNAMIC_ATOL,
            rtol=0.0,
        )
        @test !isapprox(
            measured.internal_energy_per_cell,
            shifted_kinetic_energy;
            atol=DQMC_THERMODYNAMIC_ATOL,
            rtol=0.0,
        )
    end

    @testset "Direct cyclic products" begin
        parameters = DQMCParameters(
            L=2,
            beta=0.5,
            dtau=0.1,
            U=1.5,
            mu=-0.3,
            stabilization_period=1,
            burnin_sweeps=0,
            measurement_sweeps=1,
            seed=7201,
        )
        state = initialize_dqmc_state(parameters)
        initial_direct = _direct_cyclic_green(
            state.propagators_up,
            length(state.propagators_up),
        )
        initial_error = maximum(abs.(state.greens_up .- initial_direct))
        DQMC_DIRECT_GREEN_MAX[] = max(DQMC_DIRECT_GREEN_MAX[], initial_error)
        @test initial_error < 2.0e-11

        compared_slices = Int[]
        for time_slice in state.greens_calculator_up
            jdqmcf.propagate_equaltime_greens!(
                state.greens_up,
                state.greens_calculator_up,
                state.propagators_up,
            )
            state.logdet_green_up, state.signdet_green_up,
                stabilization_error, _ = jdqmcf.stabilize_equaltime_greens!(
                    state.greens_up,
                    state.logdet_green_up,
                    state.signdet_green_up,
                    state.greens_calculator_up,
                    state.propagators_up;
                    update_B̄=false,
                )
            direct = _direct_cyclic_green(state.propagators_up, time_slice)
            mismatch = maximum(abs.(state.greens_up .- direct))
            DQMC_DIRECT_GREEN_MAX[] = max(
                DQMC_DIRECT_GREEN_MAX[],
                mismatch,
                stabilization_error,
            )
            @test mismatch < 2.0e-11
            push!(compared_slices, time_slice)
        end
        @test length(compared_slices) == 5
        @test length(unique(compared_slices)) > 1
    end

    @testset "Full determinant ratios" begin
        parameters = DQMCParameters(
            L=2,
            beta=0.5,
            dtau=0.1,
            U=1.2,
            mu=-0.2,
            stabilization_period=2,
            burnin_sweeps=0,
            measurement_sweeps=1,
            seed=7301,
        )
        state = initialize_dqmc_state(parameters)
        proposals = ((1, 1), (2, 4), (4, 7), (5, 2))
        for (time_slice, index) in proposals
            _set_state_greens_to_direct_slice!(state, time_slice)
            ratios = local_hs_flip_ratio(state, index, time_slice)
            determinant_up_before = _direct_fermion_determinant(
                state.propagators_up,
                time_slice,
            )
            determinant_down_before = _direct_fermion_determinant(
                state.propagators_down,
                time_slice,
            )
            proposed_up = deepcopy(state.propagators_up)
            proposed_down = deepcopy(state.propagators_down)
            field = state.hs_fields[index, time_slice]
            proposed_up[time_slice].expmΔτV[index] =
                exp(-state.alpha * field)
            proposed_down[time_slice].expmΔτV[index] =
                exp(+state.alpha * field)
            direct_up = _direct_fermion_determinant(
                proposed_up,
                time_slice,
            ) / determinant_up_before
            direct_down = _direct_fermion_determinant(
                proposed_down,
                time_slice,
            ) / determinant_down_before
            mismatches = (
                abs(ratios.up - direct_up),
                abs(ratios.down - direct_down),
                abs(ratios.total - direct_up * direct_down),
            )
            DQMC_RATIO_MAX[] = max(DQMC_RATIO_MAX[], mismatches...)
            @test isapprox(ratios.up, direct_up; atol=5.0e-11, rtol=5.0e-11)
            @test isapprox(
                ratios.down,
                direct_down;
                atol=5.0e-11,
                rtol=5.0e-11,
            )
            @test isapprox(
                ratios.total,
                direct_up * direct_down;
                atol=1.0e-10,
                rtol=1.0e-10,
            )
        end
    end

    @testset "Accepted rank-one updates" begin
        parameters = DQMCParameters(
            L=2,
            beta=0.5,
            dtau=0.1,
            U=1.4,
            mu=-0.25,
            stabilization_period=2,
            burnin_sweeps=0,
            measurement_sweeps=1,
            seed=7401,
        )
        state = initialize_dqmc_state(parameters)
        for (time_slice, index) in ((2, 1), (4, 6), (1, 8))
            _set_state_greens_to_direct_slice!(state, time_slice)
            old_field = state.hs_fields[index, time_slice]
            update = attempt_local_hs_flip!(
                state,
                index,
                time_slice;
                acceptance_uniform=0.0,
            )
            @test update.accepted
            @test state.hs_fields[index, time_slice] == -old_field
            direct_up = _direct_cyclic_green(
                state.propagators_up,
                time_slice,
            )
            direct_down = _direct_cyclic_green(
                state.propagators_down,
                time_slice,
            )
            mismatch = max(
                maximum(abs.(state.greens_up .- direct_up)),
                maximum(abs.(state.greens_down .- direct_down)),
            )
            DQMC_ACCEPTED_UPDATE_MAX[] = max(
                DQMC_ACCEPTED_UPDATE_MAX[],
                mismatch,
            )
            @test mismatch < 2.0e-10
            direct_log_up, direct_sign_up = logabsdet(direct_up)
            direct_log_down, direct_sign_down = logabsdet(direct_down)
            @test state.signdet_green_up == direct_sign_up
            @test state.signdet_green_down == direct_sign_down
            @test isapprox(
                state.logdet_green_up,
                direct_log_up;
                atol=2.0e-10,
                rtol=2.0e-10,
            )
            @test isapprox(
                state.logdet_green_down,
                direct_log_down;
                atol=2.0e-10,
                rtol=2.0e-10,
            )
            direct_weight_sign = sign(
                _direct_fermion_determinant(
                    state.propagators_up,
                    time_slice,
                ) *
                _direct_fermion_determinant(
                    state.propagators_down,
                    time_slice,
                ),
            )
            @test measure_basic_observables(state).configuration_sign ==
                  direct_weight_sign
        end
    end

    @testset "Sweeps and stabilization" begin
        parameters = DQMCParameters(
            L=2,
            beta=0.5,
            dtau=0.05,
            U=2.0,
            mu=-0.4,
            stabilization_period=2,
            burnin_sweeps=0,
            measurement_sweeps=1,
            seed=7501,
        )
        state = initialize_dqmc_state(parameters)
        for _ in 1:6
            sweep = sweep_all_time_slices!(state)
            recomputation = recompute_equal_time_greens!(state)
            DQMC_SWEEP_RECOMPUTE_MAX[] = max(
                DQMC_SWEEP_RECOMPUTE_MAX[],
                sweep.recompute_error,
                recomputation.maximum,
            )
            @test sweep.attempted == 2 * parameters.L^2 * 10
            @test 0 <= sweep.acceptance_rate <= 1
            @test isfinite(sweep.recompute_error)
            @test recomputation.maximum < 1.0e-9
            @test all(isfinite, state.greens_up)
            @test all(isfinite, state.greens_down)
        end
        @test state.maximum_greens_recompute_error < 1.0e-9
    end

    @testset "Smoke-run invariants and determinism" begin
        smoke_cases = (
            (
                name="u0_integration",
                U=0.0,
                mu=0.0,
                seed=8101,
            ),
            (
                name="interacting_mu0",
                U=2.0,
                mu=0.0,
                seed=8102,
            ),
            (
                name="interacting_first_vhs",
                U=2.0,
                mu=-1.897,
                seed=8103,
            ),
        )
        summaries = NamedTuple[]
        for case in smoke_cases
            parameters = DQMCParameters(
                L=2,
                beta=1.0,
                dtau=0.1,
                U=case.U,
                mu=case.mu,
                stabilization_period=5,
                burnin_sweeps=2,
                measurement_sweeps=5,
                sweeps_between_measurements=1,
                seed=case.seed,
            )
            summary = run_dqmc_smoke(parameters; case_name=case.name)
            _test_smoke_invariants(summary)
            push!(summaries, summary)
        end
        @test summaries[1].acceptance_rate == 1.0
        @test summaries[1].average_sign == 1.0

        deterministic_parameters = DQMCParameters(
            L=2,
            beta=1.0,
            dtau=0.1,
            U=2.0,
            mu=-1.897,
            stabilization_period=5,
            burnin_sweeps=2,
            measurement_sweeps=5,
            sweeps_between_measurements=1,
            seed=8201,
        )
        first_summary = run_dqmc_smoke(
            deterministic_parameters;
            case_name="deterministic",
        )
        second_summary = run_dqmc_smoke(
            deterministic_parameters;
            case_name="deterministic",
        )
        _test_smoke_invariants(first_summary)
        _test_smoke_invariants(second_summary)
        @test _without_runtime(first_summary) == _without_runtime(second_summary)
    end
end

println("DQMC_METRIC maximum_u0_green_mismatch=", DQMC_U0_GREEN_MAX[])
println("DQMC_METRIC maximum_direct_product_green_mismatch=", DQMC_DIRECT_GREEN_MAX[])
println("DQMC_METRIC maximum_determinant_ratio_mismatch=", DQMC_RATIO_MAX[])
println("DQMC_METRIC maximum_accepted_update_mismatch=", DQMC_ACCEPTED_UPDATE_MAX[])
println("DQMC_METRIC maximum_sweep_recompute_error=", DQMC_SWEEP_RECOMPUTE_MAX[])

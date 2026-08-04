using AltermagnetDQMCTrial
using LinearAlgebra: BLAS
using Printf: @printf

BLAS.set_num_threads(1)

const L = 4
const BETA = 4.0
const MU = -1.897
const DTAU = 0.1
const COMPARISON_ATOL = 1.0e-10
const THERMODYNAMIC_OBSERVABLES = (
    :density,
    :doping,
    :kinetic_energy_per_cell,
    :shifted_interaction_energy_per_cell,
    :internal_energy_per_cell,
    :grand_energy_per_cell,
)

parameters = DQMCParameters(
    L=L,
    beta=BETA,
    dtau=DTAU,
    U=0.0,
    mu=MU,
    stabilization_period=5,
    burnin_sweeps=0,
    measurement_sweeps=1,
    sweeps_between_measurements=1,
    seed=6701,
)
dimensions = validate_dqmc_parameters(parameters)
fixed_fields = fill(Int8(1), dimensions.Nsp, dimensions.Ltau)
state = initialize_dqmc_state(parameters; hs_fields=fixed_fields)

exact = finite_lattice_free_thermodynamics(L, BETA, MU)
dqmc = measure_thermodynamic_observables(state)

for observable in THERMODYNAMIC_OBSERVABLES
    exact_value = getproperty(exact, observable)
    dqmc_value = getproperty(dqmc, observable)
    absolute_error = abs(dqmc_value - exact_value)
    @printf(
        "%s exact_value=%.16g dqmc_value=%.16g absolute_error=%.16g\n",
        string(observable),
        exact_value,
        dqmc_value,
        absolute_error,
    )
    absolute_error <= COMPARISON_ATOL || error(
        "$observable differs by $absolute_error, exceeding $COMPARISON_ATOL",
    )
end

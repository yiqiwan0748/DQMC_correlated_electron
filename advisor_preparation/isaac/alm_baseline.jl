using AltermagnetDQMCTrial, Test, LinearAlgebra, Random
const TEST_ATOL = 1e-12
const TEST_RTOL = 1e-12
const root = joinpath(@__DIR__, "..", "alm")
for name in (
    "realspace_kinetic_tests.jl", "free_fermion_tests.jl",
    "two_orbital_dqmc_tests.jl", "unequal_time_greens_tests.jl",
    "two_orbital_cluster_ed_tests.jl", "two_orbital_tilted_cluster_ed_tests.jl",
    "prl_101_237004_restarted_lanczos_tests.jl",
)
    @info "Running independent baseline file" name
    include(joinpath(root, "test", name))
end

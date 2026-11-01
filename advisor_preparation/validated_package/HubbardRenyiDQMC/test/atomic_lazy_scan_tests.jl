# Targeted regression for the atomic odd-time-grid sampling failure.
# Run separately from the historical baseline suite; no archived source is modified.
module AtomicLazyScanTests
using Test, LinearAlgebra, Random, Statistics
using HubbardRenyiDQMC
const H=HubbardRenyiDQMC
module ExactScan
include(joinpath(@__DIR__,"../analysis/atomic_odd_sweep_proof.jl"))
end
@testset "Atomic systematic scan and local-lazy correction" begin
    reference=ExactScan.exact_proof()
    original=reference["ordinary_fullcycle_classes"]
    corrected=reference["lazy_local_fullcycle_classes"]
    @test sort([c["size"] for c in original]) == [1,1,30]
    @test all(c["closed"] for c in original)
    @test length(corrected)==1
    @test corrected[1]["size"]==32
    @test corrected[1]["closed"]
    @test reference["ordinary_stationarity_max_error"]<=1e-14
    @test reference["lazy_stationarity_max_error"]<=1e-14
    @test reference["lazy_minimum_diagonal_probability"]>0
    @test reference["correct_closed_class_mass"]≈0.03741523958821903 atol=1e-14
    @test reference["uniform_initialization_closed_mass"]==1/16
    closedseeds=[r["seed"] for r in reference["seed_predictions"] if !isempty(r["frozen_sites_1based"])]
    @test closedseeds==[4001,5003]
    target=reference["exact_target_totals"]
    @test target["energy_centered"]≈-1.830793928648644 atol=1e-12
    @test target["density_total"]≈4.161481159434526 atol=1e-12
    @test target["double_occupancy_total"]≈0.6230420975551021 atol=1e-12
end
@testset "Every active sampler breaks the deterministic free-field two-cycle" begin
    # U=0 makes the physical/replica weights exactly independent of every HS bit.
    # Without per-proposal laziness, a full sweep deterministically complements all
    # bits and visits only two configurations. This fixture exercises each API.
    for method in [:ordinary,:direct_fused,:rank1_fused]
        active=H.initializer(0.,0.,0.,.5,4,.1,81271;n_stab=2)
        passive=H.initializer(0.,0.,0.,.5,4,.1,81272;n_stab=2)
        seen=Set{Tuple{Vararg{Int}}}();push!(seen,Tuple(vec(active.s)))
        acceptance=Float64[]
        for _ in 1:12
            if method==:ordinary
                active,delta,a=H.update_replica!(active,0.)
            elseif method==:direct_fused
                active,delta,a=H.update_replica_fused!(active,passive,0.,[1,2],.5)
            else
                active,delta,a=H.update_replica_rank1_fused!(active,passive,0.,[1,2],.5)
            end
            push!(seen,Tuple(vec(active.s)));push!(acceptance,a)
            @test maximum(abs,active.Gup_00-.5I)<=1e-11
            @test maximum(abs,active.Gdn_00-.5I)<=1e-11
            @test isfinite(delta)&&delta<=1e-10
        end
        @test length(seen)>2
        @test 0<mean(acceptance)<1
    end
end
end

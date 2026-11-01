#!/usr/bin/env julia
# P21 adapter eligibility proof only. These are contemporaneous HS identities,
# not the external ALF comparison and not an independently counted benchmark.
using LinearAlgebra, Random, TOML, SHA, Dates
include(joinpath(@__DIR__,"UserMatrixObserver.jl"))
using .UserMatrixObserver
BLAS.set_num_threads(1)
const BASE=realpath(ENV["P21_USER_BASE"])
const OVERLAY=realpath(ENV["P21_USER_OVERLAY"])
const OUTPUT=ENV["P21_USER_PROOF_OUTPUT"]
ispath(OUTPUT)&&error("Choose a new proof receipt")
const BASESHA="b72bbbf158381b0f776c1fa259d532b069094c47e361d3490c19ea0f43ada4be"
sha(p)=bytes2hex(sha256(read(p)))
sha(joinpath(BASE,"src/dqmc_kernel.jl"))==BASESHA||error("Original kernel changed")
const OVERLAYSHA=sha(joinpath(OVERLAY,"src/dqmc_kernel.jl"))
OVERLAYSHA=="911a61642b23a4465639c3b3e4c63ff9d7cf05d51da288f3224f1775c0aba891"||error("Reviewed observer overlay differs")
module Original
include(joinpath(ENV["P21_USER_BASE"],"src/dqmc_kernel.jl"))
end
module Observed
include(joinpath(ENV["P21_USER_OVERLAY"],"src/dqmc_kernel.jl"))
end

function maxdiff(A,B)
    size(A)==size(B)||error("Buffer dimensions differ")
    maximum(abs.(A.-B))
end
function current_config_matrices(B,tau)
    L=size(first(B),1);I0=Matrix{Float64}(I,L,L)
    total=copy(I0);prefix=copy(I0)
    for (l,b) in enumerate(B)
        # Explicit product is a diagnostic only. Primary observer retains buffers.
        full=Diagonal(b.expmΔτV)*b.expmΔτK
        total=full*total
        l<=tau&&(prefix=full*prefix)
    end
    g00=inv(I0+total)
    (;g00,gt0=prefix*g00,g0t=(g00-I0)/prefix,gtt=prefix*g00/prefix)
end
function check_buffers(event)
    errors=Dict{String,Float64}()
    for spin in (:up,:dn)
        b=getproperty(event,Symbol("B",spin));exact=current_config_matrices(b,event.tau)
        for (suffix,ref) in (("00",exact.g00),("τ0",exact.gt0),("0τ",exact.g0t),("ττ",exact.gtt))
            actual=getproperty(event,Symbol("G",spin,"_",suffix))
            errors[string(spin,"_",suffix)]=maxdiff(actual,ref)
        end
        # Endpoint identities use this callback's G00, never a prior HS sample.
        g=getproperty(event,Symbol("G",spin,"_00"))
        a=getproperty(event,Symbol("G",spin,"_τ0"))
        z=getproperty(event,Symbol("G",spin,"_0τ"))
        if event.tau==0
            errors[string(spin,"_endpoint_gt0")]=maxdiff(a,g)
            errors[string(spin,"_endpoint_g0t")]=maxdiff(z,g-I)
        elseif event.tau==length(b)
            errors[string(spin,"_endpoint_gt0")]=maxdiff(a,I-g)
            errors[string(spin,"_endpoint_g0t")]=maxdiff(z,-g)
        end
    end
    errors
end
function same_state(a,b)
    for name in (:s,:Gup,:Gdn,:Gup_00,:Gup_τ0,:Gup_0τ,:Gdn_00,:Gdn_τ0,:Gdn_0τ)
        isequal(getproperty(a,name),getproperty(b,name))||error("Observer changed "*string(name))
    end
    for name in (:sgndetGup,:sgndetGdn,:logdetGup,:logdetGdn)
        isequal(getproperty(a,name),getproperty(b,name))||error("Observer changed "*string(name))
    end
    for name in (:Bup,:Bdn),l in eachindex(getproperty(a,name))
        A=getproperty(a,name)[l];B=getproperty(b,name)[l]
        isequal(A.expmΔτV,B.expmΔτV)||error("Observer changed B potential")
        isequal(A.expmΔτK,B.expmΔτK)||error("Observer changed B hopping")
    end
    for name in (:fermion_greens_calc_up,:fermion_greens_calc_dn)
        A=getproperty(a,name);B=getproperty(b,name)
        (A.forward,A.l)==(B.forward,B.l)||error("Observer changed calculator direction")
        all(isequal(x,y) for (x,y) in zip(A.B_bar,B.B_bar))||error("Observer changed B_bar")
        for k in eachindex(A.F),field in (:L,:d,:R)
            isequal(getproperty(A.F[k],field),getproperty(B.F[k],field))||error("Observer changed LDR factors")
        end
        isequal(A.G′,B.G′)||error("Observer changed stabilization scratch Green matrix")
    end
    rand(copy(a.rng),UInt64)==rand(copy(b.rng),UInt64)||error("Observer consumed RNG")
end
function main()
    seed=20261022021;cycles=16;tol=5e-10
    a=Original.initializer(1.,2.,0.,1.,4,.1,seed;n_stab=2,boundary_condition=:pbc)
    b=deepcopy(a);events=Any[];worst=0.;pre_worst=0.;post_worst=0.;initial_worst=0.
    rawbin=MatrixBin(4,10)
    observer=function(event)
        observe!(rawbin,event)
        errors=check_buffers(event);e=maximum(values(errors));worst=max(worst,e)
        if event.stage=="pre_stabilization";pre_worst=max(pre_worst,e)
        elseif event.stage=="post_stabilization";post_worst=max(post_worst,e)
        else;initial_worst=max(initial_worst,e);end
        push!(events,Dict("tau_index"=>event.tau,"stage"=>event.stage,"forward"=>event.forward,
            "phase_real"=>real(event.phase),"phase_imag"=>imag(event.phase),"matrix_errors"=>errors,
            "hs_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(copy(event.s)))))))
        nothing
    end
    for cycle in 1:cycles
        a.fermion_greens_calc_up.forward||error("Cycle must begin forward")
        for direction in 1:2
            a,da,aa=Original.update_replica!(a,0.)
            b,db,ab=Observed.update_replica!(b,0.;observer=observer)
            isequal((da,aa),(db,ab))||error("Observer changed acceptance/stability")
            same_state(a,b)
        end
    end
    expected=cycles*(1+2*a.Lτ)
    length(events)==expected||error("Incomplete callback grid")
    all(x["forward"] for x in events)||error("Unexpected reverse event")
    pass=isfinite(worst)&&worst<=tol
    receipt=Dict("schema"=>"p21_user_native_buffer_eligibility_v1","catalogue_id"=>"P21",
        "status"=>pass ? "ADAPTER_PROOF_PASS_NOT_EXTERNAL_BENCHMARK" : "ADAPTER_PROOF_FAIL",
        "created_utc"=>string(now(UTC)),"original_source_sha256"=>BASESHA,"observer_source_sha256"=>OVERLAYSHA,
        "driver_sha256"=>sha(@__FILE__),"julia_version"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),
        "seed"=>seed,"forward_reverse_cycles"=>cycles,"tolerance"=>tol,"max_error"=>worst,
        "pre_stabilization_max_error"=>pre_worst,"post_stabilization_max_error"=>post_worst,
        "initial_max_error"=>initial_worst,"trajectory_and_rng_exactly_matched"=>true,
        "callbacks"=>events,"native_buffer_accumulator"=>raw_record(rawbin),
        "job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"))
    mkpath(dirname(OUTPUT));open(OUTPUT,"w")do io;TOML.print(io,receipt);end
    pass||error("Actual maintained buffers failed contemporaneous-HS identities")
end
main()

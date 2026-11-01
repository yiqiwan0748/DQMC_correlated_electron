module UserMatrixObserver
using LinearAlgebra
export MatrixBin, observe!, raw_record

"""Read-only accumulator of actual maintained forward buffers, before stabilization."""
mutable struct MatrixBin
    L::Int
    M::Int
    nums::Array{ComplexF64,5} # order, tau+1, spin, i, j
    dens::Vector{Float64}    # moving-cut configurations may differ with tau
    counts::Vector{Int}
    endpoint_max_error::Float64
end
MatrixBin(L,M)=MatrixBin(L,M,zeros(ComplexF64,2,M+1,2,L,L),zeros(M+1),zeros(Int,M+1),0.)

function observe!(b::MatrixBin,event)
    event.forward||error("P21 observer receives forward buffers only")
    # Primary measurements exercise the user's rank-one updates before rebuilding.
    event.stage=="initial"&&event.tau==0 || event.stage=="pre_stabilization" || return nothing
    tau=event.tau;0<=tau<=b.M||error("Tau out of range")
    abs(imag(event.phase))<=1e-12||error("This real-Hamiltonian adapter requires real sign")
    sign=real(event.phase);isfinite(sign)&&abs(abs(sign)-1)<=1e-12||error("Invalid phase")
    for (spin,prefix) in enumerate((:up,:dn))
        gt0=getproperty(event,Symbol("G",prefix,"_τ0"))
        g0t=getproperty(event,Symbol("G",prefix,"_0τ"))
        g00=getproperty(event,Symbol("G",prefix,"_00"))
        size(gt0)==(b.L,b.L)&&size(g0t)==(b.L,b.L)||error("Full matrix required")
        all(isfinite,gt0)&&all(isfinite,g0t)||error("Nonfinite maintained buffer")
        b.nums[1,tau+1,spin,:,:].+=sign.*gt0
        b.nums[2,tau+1,spin,:,:].+=sign.*g0t
        if tau==0
            b.endpoint_max_error=max(b.endpoint_max_error,maximum(abs.(gt0.-g00)),maximum(abs.(g0t.-(g00-I))))
        elseif tau==b.M
            b.endpoint_max_error=max(b.endpoint_max_error,maximum(abs.(gt0.-(I-g00))),maximum(abs.(g0t.+g00)))
        end
    end
    b.dens[tau+1]+=sign;b.counts[tau+1]+=1
    nothing
end

function raw_record(b::MatrixBin)
    all(==(first(b.counts)),b.counts)&&first(b.counts)>0||error("Incomplete primary time grid")
    # Explicit loops avoid Julia column-major/Fortran export ambiguity.
    arrays=Dict{String,Any}()
    for (order,name) in enumerate(("GT0","G0T"))
        v=ComplexF64[b.nums[order,t+1,s,i,j] for t in 0:b.M for s in 1:2 for i in 1:b.L for j in 1:b.L]
        arrays[name*"_numerator_real_sums"]=real.(v)
        arrays[name*"_numerator_imag_sums"]=imag.(v)
    end
    merge(arrays,Dict("sign_denominator_sums_by_tau"=>copy(b.dens),"counts_by_tau"=>copy(b.counts),
        "endpoint_contemporaneous_max_error"=>b.endpoint_max_error,
        "ordering"=>"tau-major,spin(up,down),row_i,column_j; all indices zero-based externally",
        "primary_stage"=>"initial(tau0) and pre_stabilization(tau1..M)",
        "measurement_origin"=>"actual user maintained buffers; no dense recomputation"))
end
end

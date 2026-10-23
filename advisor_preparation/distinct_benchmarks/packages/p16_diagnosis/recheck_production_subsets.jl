#!/usr/bin/env julia
using LinearAlgebra,TOML,SHA
const ROOT=ENV["RENYI_ROOT"]
module ED
include(joinpath(Main.ROOT,"analysis/compute_s2_ed.jl"))
end
BLAS.set_num_threads(1)
function main()
    out=ENV["P16_SUBSETS_OUTPUT"];ispath(out)&&error("Output exists")
    space=ED.build_space(4,1.,4.,0.,2,2;boundary="pbc")
    E,psi,it,reported=ED.ground_state_lanczos(space;tol=1e-12,maxiter=1000)
    hv=similar(psi);ED.apply_hamiltonian!(hv,psi,space);res=norm(hv-E*psi);res<=1e-9||error("Ground residual")
    rows=[Dict("sites"=>[i for i in 0:3 if !iszero(bits&(1<<i))],"purity"=>exp(-ED.renyi2_entropy_sites(psi,space,[i for i in 0:3 if !iszero(bits&(1<<i))]))) for bits in 0:15]
    d=Dict("source_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"analysis/compute_s2_ed.jl")))),"driver_sha256"=>bytes2hex(sha256(read(@__FILE__))),"ground_residual"=>res,"energy"=>E,"rows"=>rows,"job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"))
    open(out,"w") do io;TOML.print(io,d);end
end
main()

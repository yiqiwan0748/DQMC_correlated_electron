#!/usr/bin/env julia
using LinearAlgebra,TOML,SHA,Dates
const ROOT=ENV["RENYI_ROOT"]
module ED
include(joinpath(Main.ROOT,"analysis/compute_s2_ed.jl"))
end
include(joinpath(@__DIR__,"regrouping_v2.jl"))
using .P16RegroupingV2
BLAS.set_num_threads(1)
function main()
    output=ENV["P16_RECHECK_OUTPUT"];ispath(output)&&error("Output exists")
    oraclepath=ENV["P16_ORACLE_TOML"];oracle=TOML.parsefile(oraclepath)["P16"]["values"]
    space=ED.build_space(4,1.,4.,0.,2,2;boundary="pbc")
    E,psi,it,reported=ED.ground_state_lanczos(space;tol=1e-12,maxiter=1000)
    hv=similar(psi);ED.apply_hamiltonian!(hv,psi,space);res=norm(hv-E*psi);res<=1e-9||error("Ground residual")
    old=[ED.renyi2_entropy_sites(psi,space,s) for s in [[0],[2],[0,2]]]
    ent=[renyi2_entropy_sites_v2(psi,space,s) for s in [[0],[2],[0,2]]]
    values=vcat(exp.(-ent),[ent[1]+ent[2]-ent[3]]);err=maximum(abs.(values.-oracle))
    rows=Any[]
    for bits in 0:15
        A=[i for i in 0:3 if !iszero(bits&(1<<i))]
        push!(rows,Dict("sites"=>A,"old_purity"=>exp(-ED.renyi2_entropy_sites(psi,space,A)),"corrected_purity"=>exp(-renyi2_entropy_sites_v2(psi,space,A))))
    end
    d=Dict("id"=>"P16","status"=>err<=1e-8 ? "PASS" : "FAIL","scope"=>"Actual frozen user ED ground state plus isolated candidate regrouping_v2 adapter; canonical source unchanged","values"=>values,"oracle_values"=>oracle,"old_values"=>vcat(exp.(-old),[old[1]+old[2]-old[3]]),"maximum_absolute_error"=>err,"all_subset_purities"=>rows,"ground_residual"=>res,"energy"=>E,"job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"),"source_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"analysis/compute_s2_ed.jl")))),"candidate_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"regrouping_v2.jl")))),"driver_sha256"=>bytes2hex(sha256(read(@__FILE__))),"oracle_sha256"=>bytes2hex(sha256(read(oraclepath))))
    open(output,"w") do io;TOML.print(io,d);end
    println("P16 isolated v2 ",d["status"]," maximum error ",err)
end
main()

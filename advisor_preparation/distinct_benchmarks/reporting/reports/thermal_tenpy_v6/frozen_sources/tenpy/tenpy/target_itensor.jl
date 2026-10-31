#!/usr/bin/env julia
# Additive observable/sector adapter around the user's existing Hubbard MPO.
using LinearAlgebra, Random, TOML, SHA, Dates
using ITensors, ITensorMPS
BLAS.set_num_threads(1)

function options(args)
    d=Dict{String,String}(); pilot=false; i=1
    while i<=length(args)
        if args[i]=="--pilot"
            pilot=true; i+=1
        else
            i<length(args) || error("Missing value for $(args[i])")
            d[args[i]]=args[i+1]; i+=2
        end
    end
    return d,pilot
end
const OPTS,PILOT=options(ARGS)
const MANIFEST=abspath(get(OPTS,"--manifest",joinpath(@__DIR__,"manifest.toml")))
const CFG=TOML.parsefile(MANIFEST)
const OUTPUT=abspath(OPTS["--output"])
ispath(OUTPUT) && error("Output exists: $OUTPUT")
const TARGET_SOURCE=haskey(ENV,"RENYI_ROOT") ? joinpath(ENV["RENYI_ROOT"],"analysis/compute_s2_itensor_dmrg.jl") : normpath(joinpath(dirname(MANIFEST),CFG["source_relative"]))
bytes2hex(sha256(read(TARGET_SOURCE)))==CFG["source_sha256"] || error("Target source SHA256 mismatch")
module UserDMRG
include(Main.TARGET_SOURCE)
end

function prefix_sectors(psi,cut,floor)
    orthogonalize!(psi,cut)
    sites=siteinds(psi)
    A=psi[1]
    for j in 2:cut; A*=psi[j]; end
    ordered=(sites[1:cut]...,linkind(psi,cut))
    M=reshape(Array(A,ordered...),4^cut,:)
    tracecheck=sum(abs2,M)
    n=op("Ntot",sites[1]); z=op("Sz",sites[1])
    nq=[round(Int,real(n[i,i])) for i in 1:4]
    sq=[round(Int,2real(z[i,i])) for i in 1:4]
    groups=Dict{Tuple{Int,Int},Vector{Int}}()
    for row in 0:(4^cut-1)
        ds=digits(row;base=4,pad=cut).+1
        q=(sum(nq[d] for d in ds),sum(sq[d] for d in ds))
        push!(get!(groups,q,Int[]),row+1)
    end
    result=Dict{String,Any}[]
    for ((N,S),rows) in sort(collect(groups);by=first)
        p=sort(svdvals(Matrix(M[rows,:])).^2;rev=true)./tracecheck
        push!(result,Dict("N"=>N,"twoSz"=>S,"weight"=>sum(p),
                          "probabilities"=>p[p.>floor]))
    end
    result,tracecheck
end

function parity_strings(psi,sites,L)
    rows=Dict{String,Any}[]
    for i in 1:L,j in i:L
        terms=Any[isodd(j-i+1) ? -1.0 : 1.0]
        for k in i:j; append!(terms,("F",k)); end
        os=OpSum();os+=Tuple(terms)
        product=MPO(os,sites)
        push!(rows,Dict("i"=>i,"j"=>j,"value"=>real(inner(psi',product,psi))))
    end
    rows
end

function solve(c,arm,floor)
    start=time()
    L=c["L"];sites=siteinds("Electron",L;conserve_qns=true)
    H=UserDMRG.hubbard_mpo(sites;L,t=c["t"],U=c["U"],mu=0.0,boundary="obc")
    state=vcat(fill("Up",c["Nup"]),fill("Dn",c["Ndown"]),fill("Emp",L-c["Nup"]-c["Ndown"]))
    Random.seed!(arm["seed"]);shuffle!(state)
    psi0=random_mps(sites,state;linkdims=16)
    ns=arm["nsweeps"];cap=arm["chi_max"]
    schedule=[min(cap,n<=5 ? (16,32,64,128,256)[n] : cap) for n in 1:ns]
    noise=[n<=3 ? (1e-6,1e-7,1e-8)[n] : 0.0 for n in 1:ns]
    observer=DMRGObserver()
    reported,psi=dmrg(H,psi0;nsweeps=ns,maxdim=schedule,cutoff=arm["cutoff"],
                      eigsolve_tol=1e-12,noise,observer,outputlevel=1)
    normcheck=real(inner(psi,psi))
    E=real(inner(psi',H,psi))/normcheck
    variance=real(inner(H,psi,H,psi))/normcheck-E^2
    density=real.(expect(psi,"Ntot"));spin=real.(expect(psi,"Sz"))
    sectors,tracecheck="P11" in c["catalogue_ids"] ? prefix_sectors(psi,c["cut"],floor) : (Dict{String,Any}[],1.0)
    strings="P15" in c["catalogue_ids"] ? parity_strings(psi,sites,L) : Dict{String,Any}[]
    te=truncerrors(observer)
    Dict{String,Any}("energy"=>E,"reported_energy"=>reported,"variance"=>variance,
       "norm"=>normcheck,"N"=>sum(density),"twoSz"=>2sum(spin),"density"=>density,
       "sectors"=>sectors,"strings"=>strings,"rdm_trace"=>tracecheck,
       "max_bond"=>maxlinkdim(psi),"final_truncerr"=>isempty(te) ? -1.0 : last(te),
       "elapsed_seconds"=>time()-start)
end

function main()
    result=Dict{String,Any}("schema"=>CFG["schema"],"backend"=>"ITensorMPS",
        "status"=>"RUNNING","pilot"=>PILOT,"manifest_sha256"=>bytes2hex(sha256(read(MANIFEST))),
        "target_source_sha256"=>bytes2hex(sha256(read(TARGET_SOURCE))),
        "adapter_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "julia_version"=>string(VERSION),"library_versions"=>Dict("ITensorMPS"=>string(pkgversion(ITensorMPS)),"ITensors"=>string(pkgversion(ITensors))),"cases"=>Dict{String,Any}(),
        "attribution"=>"User Hubbard MPO and DMRG workflow; ITensorMPS supplies the optimizer. This file adds doped initialization and observables.")
    start=time()
    try
        for name in ("doped","half")
            PILOT && name!="doped" && continue
            c=CFG["fixtures"][name];arms=Dict{String,Any}()
            result["cases"][name]=Dict("arms"=>arms)
            for armname in ("coarse","refined")
                PILOT && armname!="coarse" && continue
                arms[armname]=solve(c,CFG["arms"][armname],CFG["schmidt_probability_floor"])
            end
        end
        result["status"]=PILOT ? "PILOT_PARTIAL" : "COMPLETED"
    catch e
        result["status"]="EXECUTION_ERROR"
        result["error"]=sprint(showerror,e,catch_backtrace())
    end
    result["elapsed_seconds"]=time()-start
    mkpath(dirname(OUTPUT))
    open(OUTPUT,"w") do io;TOML.print(io,result);end
    result["status"]=="EXECUTION_ERROR" && error(result["error"])
end
main()

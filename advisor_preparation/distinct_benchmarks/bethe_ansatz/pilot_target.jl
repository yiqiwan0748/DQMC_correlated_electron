#!/usr/bin/env julia
# Existing user Hubbard MPO plus ITensorMPS optimizer; additive sector/observable adapter.
using LinearAlgebra, Random, TOML, SHA
using ITensors, ITensorMPS
BLAS.set_num_threads(1)
function options(args)
    length(args)%2==0 || error("Expected --manifest PATH --output PATH")
    Dict(args[i]=>args[i+1] for i in 1:2:length(args))
end
const OPT=options(ARGS)
const CONFIG_PATH=abspath(get(OPT,"--manifest",joinpath(@__DIR__,"pilot.toml")))
const CFG=TOML.parsefile(CONFIG_PATH)
const OUTPUT=abspath(OPT["--output"])
ispath(OUTPUT) && error("Output exists")
CFG["phase"]=="PILOT_ONLY" && !CFG["scientific_completion_allowed"] || error("Pilot protocol only")
string(VERSION)==CFG["julia_version"] || error("Julia version mismatch")
string(pkgversion(ITensors))==CFG["itensors_version"] || error("ITensors version mismatch")
string(pkgversion(ITensorMPS))==CFG["itensormps_version"] || error("ITensorMPS version mismatch")
const SOURCE=haskey(ENV,"RENYI_ROOT") ? joinpath(ENV["RENYI_ROOT"],"analysis/compute_s2_itensor_dmrg.jl") : normpath(joinpath(dirname(CONFIG_PATH),CFG["source_relative"]))
const SOURCE_ROOT=dirname(dirname(SOURCE))
for (file,key) in ((SOURCE,"source_sha256"),(joinpath(SOURCE_ROOT,"Project.toml"),"project_sha256"),(joinpath(SOURCE_ROOT,"Manifest.toml"),"manifest_sha256"))
    bytes2hex(sha256(read(file)))==CFG[key] || error("Source/environment SHA mismatch: $key")
end
module UserDMRG
include(Main.SOURCE)
end

function initial_state(L,sector,seed)
    iseven(L) || error("Even chain required")
    state=[isodd(j) ? "Up" : "Dn" for j in 1:L]
    sector=="minus" && (state[2]="Emp")
    sector=="plus" && (state[2]="UpDn")
    sector in ("half","minus","plus") || error("Unknown sector")
    Random.seed!(seed);shuffle!(state)
    m=L÷2
    Nup=sector=="plus" ? m+1 : m
    Ndown=sector=="minus" ? m-1 : m
    state,Nup,Ndown
end

function solve(L,sector,armname,arm)
    began=time();sector_id=findfirst(==(sector),CFG["sectors"])
    seed=arm["seed_base"]+1000L+10sector_id
    state,Nup,Ndown=initial_state(L,sector,seed)
    sites=siteinds("Electron",L;conserve_qns=true)
    H=UserDMRG.hubbard_mpo(sites;L,t=CFG["t"],U=CFG["U"],mu=CFG["mu_physical"],boundary=CFG["boundary"])
    psi0=random_mps(sites,state;linkdims=16)
    cap=arm["chi_max"];ns=arm["nsweeps"]
    schedule=[min(cap,k<=5 ? (16,32,64,128,256)[k] : cap) for k in 1:ns]
    noise=[k<=3 ? (1e-6,1e-7,1e-8)[k] : 0.0 for k in 1:ns]
    observer=DMRGObserver()
    reported,psi=dmrg(H,psi0;nsweeps=ns,maxdim=schedule,cutoff=arm["cutoff"],
                      eigsolve_tol=1e-12,noise,observer,outputlevel=1)
    normcheck=real(inner(psi,psi));energy=real(inner(psi',H,psi))/normcheck
    variance=real(inner(H,psi,H,psi))/normcheck-energy^2
    up=real.(expect(psi,"Nup"));down=real.(expect(psi,"Ndn"));doubles=real.(expect(psi,"Nupdn"))
    te=truncerrors(observer)
    Dict{String,Any}("L"=>L,"sector"=>sector,"arm"=>armname,"seed"=>seed,
        "Nup_expected"=>Nup,"Ndown_expected"=>Ndown,"Nup"=>sum(up),"Ndown"=>sum(down),
        "energy"=>energy,"reported_energy"=>reported,"energy_per_site"=>energy/L,
        "norm"=>normcheck,"variance"=>variance,"variance_per_site"=>variance/L,
        "double_occupancy_per_site"=>sum(doubles)/L,"site_double_occupancy"=>doubles,
        "site_Nup"=>up,"site_Ndown"=>down,"max_bond"=>maxlinkdim(psi),
        "final_truncerr"=>isempty(te) ? -1.0 : last(te),"sweep_truncerrors"=>te,
        "elapsed_seconds"=>time()-began)
end

function save_new(path,payload)
    ispath(path) && error("Receipt exists: $path")
    mkpath(dirname(path));open(path,"w") do io;TOML.print(io,payload);end
end
function main()
    result=Dict{String,Any}("schema"=>CFG["schema"],"phase"=>"PILOT_ONLY","status"=>"RUNNING",
        "scientific_completion_allowed"=>false,"manifest_sha256"=>bytes2hex(sha256(read(CONFIG_PATH))),
        "source_sha256"=>bytes2hex(sha256(read(SOURCE))),"adapter_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "project_sha256"=>bytes2hex(sha256(read(joinpath(SOURCE_ROOT,"Project.toml")))),
        "environment_manifest_sha256"=>bytes2hex(sha256(read(joinpath(SOURCE_ROOT,"Manifest.toml")))),
        "julia_version"=>string(VERSION),"ITensors_version"=>string(pkgversion(ITensors)),
        "ITensorMPS_version"=>string(pkgversion(ITensorMPS)),"rows"=>Dict{String,Any}[],
        "attribution"=>"Actual existing user hubbard_mpo; ITensorMPS supplies DMRG. Adapter adds N±1 initial states, local double occupancy, and quality/resource receipts. No thermodynamic completion from pilot.")
    began=time()
    try
        for L in CFG["L"],sector in CFG["sectors"],armname in ("coarse","refined")
            row=solve(L,sector,armname,CFG["arms"][armname]);push!(result["rows"],row)
            save_new(joinpath(dirname(OUTPUT),"L$(L)_$(sector)_$(armname).toml"),row)
        end
        result["status"]="PILOT_ONLY_COMPLETE"
    catch e
        result["status"]="EXECUTION_ERROR";result["error"]=sprint(showerror,e,catch_backtrace())
    end
    result["elapsed_seconds"]=time()-began;save_new(OUTPUT,result)
    result["status"]=="EXECUTION_ERROR" && error(result["error"])
end
main()

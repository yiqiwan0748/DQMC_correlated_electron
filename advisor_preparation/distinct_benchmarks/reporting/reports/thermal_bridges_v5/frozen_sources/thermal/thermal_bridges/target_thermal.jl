#!/usr/bin/env julia
# Calls the existing user thermal-ED and prefix-RDM routines; no kernel edits.
using LinearAlgebra, TOML, SHA
BLAS.set_num_threads(1)
function options(args)
    length(args)%2==0 || error("Use --manifest PATH --output PATH")
    Dict(args[i]=>args[i+1] for i in 1:2:length(args))
end
const OPTS=options(ARGS)
const MANIFEST=abspath(get(OPTS,"--manifest",joinpath(@__DIR__,"manifest.toml")))
const CFG=TOML.parsefile(MANIFEST)
string(VERSION)==CFG["julia_version"] || error("Unpinned Julia version")
const OUTPUT=abspath(OPTS["--output"])
ispath(OUTPUT) && error("Output exists: $OUTPUT")
const THERMAL_SOURCE=haskey(ENV,"RENYI_ROOT") ? joinpath(ENV["RENYI_ROOT"],"analysis/audit_u4_l4_thermal_ed.jl") : normpath(joinpath(dirname(MANIFEST),CFG["thermal_source_relative"]))
const ED_SOURCE=joinpath(dirname(THERMAL_SOURCE),"compute_s2_ed.jl")
bytes2hex(sha256(read(THERMAL_SOURCE)))==CFG["thermal_source_sha256"] || error("Thermal source SHA mismatch")
bytes2hex(sha256(read(ED_SOURCE)))==CFG["ed_source_sha256"] || error("ED source SHA mismatch")
include(THERMAL_SOURCE)
const A=U4L4ThermalEDAudit

function compute_bridges(cfg)
    L=cfg["L"];t=cfg["t"];U=cfg["U"];nA=cfg["nA"]
    L==4 && cfg["boundary"]=="pbc" || error("Frozen fixture mismatch")
    K,V,masks=A.full_fock_matrices(L,t,U,0.0;boundary="pbc")
    numbers=[count_ones(u)+count_ones(d) for u in masks for d in masks]
    ups=[count_ones(u) for u in masks for d in masks]
    downs=[count_ones(d) for u in masks for d in masks]
    N=Matrix(Diagonal(Float64.(numbers)))
    Nu=Matrix(Diagonal(Float64.(ups)));Nd=Matrix(Diagonal(Float64.(downs)))
    # Existing source uses U(nup-1/2)(ndn-1/2)-mu_centered*N.
    H=K+Diagonal(V)+U/2*N-U*L/4*I
    H=Matrix(H)
    blocks=Dict{Int,Any}()
    for number in 0:2L
        ids=findall(==(number),numbers);decomp=eigen(Symmetric(H[ids,ids]))
        blocks[number]=(ids,decomp.values,decomp.vectors)
    end
    levels=vcat([blocks[n][2] for n in 0:2L]...)
    charges=vcat([fill(n,length(blocks[n][2])) for n in 0:2L]...)
    diagnostics=Dict{String,Any}("maximum_trace_error"=>0.,"maximum_rdm_trace_error"=>0.,
        "hermiticity_error"=>norm(H-H'),"N_commutator_error"=>norm(H*N-N*H))
    expect(rho,op)=real(tr(rho*op))
    function density(beta,mu)
        kinetic,potential,_=A.full_fock_matrices(L,t,U,mu-U/2;boundary="pbc")
        rho,d=A.exact_thermal_density(kinetic+Diagonal(potential),beta)
        weights=exp.(-beta.*(d.values.-d.values[1]));weights./=sum(weights)
        diagnostics["maximum_trace_error"]=max(diagnostics["maximum_trace_error"],abs(real(tr(rho))-1))
        rho,d,weights
    end
    function purity(rho)
        s,p,r=A.renyi2(rho,masks,L,nA)
        diagnostics["maximum_rdm_trace_error"]=max(diagnostics["maximum_rdm_trace_error"],abs(real(tr(r))-1))
        p,s
    end
    numberweight(rho,n)=real(sum(diag(rho)[numbers.==n]))
    spinweight(rho,up,dn)=real(sum(diag(rho)[(ups.==up).&(downs.==dn)]))
    function moments(rho,mu)
        G=H-mu*N;E=expect(rho,H);Nm=expect(rho,N);Gm=expect(rho,G)
        E2=expect(rho,H*H);N2=expect(rho,N*N);G2=expect(rho,G*G)
        EN=expect(rho,H*N);EG=expect(rho,H*G)
        Dict{String,Any}("Ephysical"=>E,"Nmean"=>Nm,"Gphysical"=>Gm,"E2"=>E2,"N2"=>N2,"G2"=>G2,"EN"=>EN,"EG"=>EG,
          "var_E"=>E2-E^2,"var_N"=>N2-Nm^2,"var_G"=>G2-Gm^2,"cov_EN"=>EN-E*Nm,"cov_EG"=>EG-E*Gm)
    end
    spec=cfg["X02"];mu=spec["mu_physical"];N0=spec["N_ground"]
    _,d,_=density(1.0,mu);psi=d.vectors[:,1];rho0=psi*psi'
    p0,s0=purity(rho0)
    fixed=A.fixed_sector_ground_state(L,t,U,0.0,nA;boundary="pbc")
    ground=Dict{String,Any}("Ephysical"=>expect(rho0,H),"S2"=>s0,"purity"=>p0,"Nmean"=>expect(rho0,N),
      "Nup_mean"=>expect(rho0,Nu),"Ndown_mean"=>expect(rho0,Nd),"global_grand_gap"=>d.values[2]-d.values[1],
      "ground_residual"=>norm((H-mu*N)*psi-(d.values[1]-U*L/4)*psi),
      "mu_minus"=>blocks[N0][2][1]-blocks[N0-1][2][1],"mu_plus"=>blocks[N0+1][2][1]-blocks[N0][2][1],
      "fixed_Nup_Ndown_ground_S2"=>fixed.S2,"fixed_Nup_Ndown_ground_Ephysical"=>fixed.energy,
      "fixed_Nup_Ndown_ground_residual"=>fixed.residual)
    rows=Dict{String,Any}[]
    for beta in spec["beta"]
        rho,_,p=density(beta,mu);pur,s2=purity(rho);epsilon=sum(p[2:end]);B=spec["purity_lipschitz_factor"]*epsilon
        valid=B<p0;bound=valid ? -log1p(-B/p0) : -1.0
        wN=numberweight(rho,N0);delta=s2-s0
        qualifies=valid && bound<=spec["entropy_bias_target"] && max(0.,1-wN)<=spec["particle_sector_leakage_target"] && abs(delta)<=spec["entropy_bias_target"]
        push!(rows,Dict("beta"=>beta,"mu_physical"=>mu,"mu_centered"=>mu-U/2,"S2"=>s2,"purity"=>pur,
          "delta_S2"=>delta,"Nmean"=>expect(rho,N),"weight_N_ground"=>wN,
          "weight_fixed_Nup_Ndown"=>spinweight(rho,spec["Nup_ground"],spec["Ndown_ground"]),
          "ground_weight"=>p[1],"excited_weight"=>epsilon,"entropy_bound_valid"=>valid,"entropy_bias_bound"=>bound,"qualifies"=>qualifies))
    end
    X02=Dict("ground"=>ground,"rows"=>rows)
    spec=cfg["X03"];Nc=spec["N_total"];ids=blocks[Nc][1];rows=Dict{String,Any}[]
    for beta in spec["beta"]
        blockrho,_=A.exact_thermal_density(H[ids,ids],beta)
        canonical=zeros(size(H));canonical[ids,ids]=blockrho
        lo,hi=spec["mu_physical_bracket"]
        function meanN(mu)
            w=exp.(-beta.*((levels.-mu.*charges).-minimum(levels.-mu.*charges)));w./=sum(w)
            m=sum(w.*charges);m,beta*(sum(w.*charges.^2)-m^2)
        end
        first(meanN(lo))<Nc<first(meanN(hi)) || error("Chemical-potential root not bracketed")
        iterations=0
        for iteration in 1:spec["max_bisections"]
            iterations=iteration;mid=(lo+hi)/2
            if first(meanN(mid))<Nc;lo=mid;else;hi=mid;end
            hi-lo<=spec["mu_bracket_width"] && break
        end
        mu=(lo+hi)/2;gce,_,_=density(beta,mu)
        pcA,scA=purity(canonical);pgA,sgA=purity(gce)
        cE=expect(canonical,H)/L;gE=expect(gce,H)/L
        push!(rows,Dict("beta"=>beta,"mu_physical"=>mu,"mu_centered"=>mu-U/2,"Nmean_GCE"=>expect(gce,N),
          "Nmean_canonical"=>expect(canonical,N),"Nvariance_canonical"=>expect(canonical,N*N)-Nc^2,
          "compressibility"=>last(meanN(mu)),"mean_N_error"=>abs(expect(gce,N)-Nc),"bisections"=>iterations,
          "mu_bracket_width"=>hi-lo,"canonical_dimension"=>length(ids),
          "fixed_Nup2_Ndown1_dimension"=>count((ups.==2).&(downs.==1)),
          "canonical_weight_Nup2_Ndown1"=>spinweight(canonical,2,1),
          "Ecanonical_per_site"=>cE,"EGCE_per_site"=>gE,"delta_energy_per_site"=>gE-cE,
          "purity_canonical"=>pcA,"purity_GCE"=>pgA,"delta_purity"=>pgA-pcA,"S2_canonical"=>scA,"S2_GCE"=>sgA))
    end
    X03=Dict("rows"=>rows,"ensemble"=>"fixed TOTAL N; all Sz sectors included")
    spec=cfg["P07"];mu=spec["mu_physical"];rows=Dict{String,Any}[]
    for beta in spec["beta"]
        rho,_,_=density(beta,mu);m=moments(rho,mu);T=1/beta;h=T*spec["finite_difference_relative_T_step"]
        energyat(temp)=expect(first(density(1/temp,mu)),H)
        fd=(energyat(T-2h)-8energyat(T-h)+8energyat(T+h)-energyat(T+2h))/(12h)
        response=beta^2*m["cov_EG"]
        push!(rows,merge(m,Dict("beta"=>beta,"mu_physical"=>mu,"mu_centered"=>mu-U/2,
          "C_mu"=>response,"beta2_var_grand"=>beta^2*m["var_G"],"C_mu_finite_difference"=>fd,
          "finite_difference_error"=>abs(response-fd),
          "covariance_identity_error"=>abs(m["cov_EG"]-(m["var_E"]-mu*m["cov_EN"])))))
    end
    Dict("X02"=>X02,"X03"=>X03,"P07"=>Dict("rows"=>rows),"diagnostics"=>diagnostics)
end

function main()
    start=time();result=Dict{String,Any}("schema"=>CFG["schema"],"status"=>"COMPLETED",
      "manifest_sha256"=>bytes2hex(sha256(read(MANIFEST))),"adapter_sha256"=>bytes2hex(sha256(read(@__FILE__))),
      "thermal_source_sha256"=>bytes2hex(sha256(read(THERMAL_SOURCE))),"ed_source_sha256"=>bytes2hex(sha256(read(ED_SOURCE))),
      "julia_version"=>string(VERSION),"attribution"=>"Existing user full_fock_matrices, exact_thermal_density, renyi2 and fixed_sector_ground_state; additive ensemble/response adapter.")
    merge!(result,compute_bridges(CFG));result["elapsed_seconds"]=time()-start
    mkpath(dirname(OUTPUT));open(OUTPUT,"w") do io;TOML.print(io,result);end
end
main()

#!/usr/bin/env julia
using LinearAlgebra, TOML, SHA, Dates
BLAS.set_num_threads(1)
const HERE=@__DIR__
const MANIFEST=joinpath(HERE,"manifest.toml")
const M=TOML.parsefile(MANIFEST)
const SOURCE=joinpath(HERE,M["source_root"])
string(VERSION)==M["julia_version"]||error("Unpinned Julia version before source load")
for (path,expected) in M["source_hashes"]
    bytes2hex(sha256(read(joinpath(SOURCE,path))))==expected||error("Source hash mismatch before include: $path")
end
include(joinpath(SOURCE,"scripts/two_orbital_open_cluster_ed.jl"))
include(joinpath(SOURCE,"scripts/two_orbital_atom_ed.jl"))
const ED=TwoOrbitalOpenClusterED
const ATOM=TwoOrbitalAtomED
const BAND=ED.ValidatedKinetic
sha(path)=bytes2hex(sha256(read(path)))

function eigengroups(es,vs,tol)
    rows=Any[];start=1
    while start<=length(es)
        stop=start
        while stop<length(es)&&abs(es[stop+1]-es[start])<=tol;stop+=1;end
        p=vs[:,start:stop]*adjoint(vs[:,start:stop])
        push!(rows,Dict("energies"=>es[start:stop],"multiplicity"=>stop-start+1,
              "projector"=>[real(p[i,j]) for i in axes(p,1) for j in axes(p,2)],
              "site_probabilities"=>real.(diag(p)),
              "orbital_weights"=>[real(sum(p[i,i] for i in a:2:8)) for a in 1:2]))
        start=stop+1
    end
    rows
end

function a15()
    p=M["A15"];parameters=BAND.ModelParameters(p["hoppings"]...)
    k=ED.open_cluster_kinetic_matrix(2,2,parameters);rows=Any[]
    for u in p["interactions"],sector in p["spin_sectors"]
        nu,nd=sector;h=ED.sector_hamiltonian(k,u,nu,nd);dense=ED.dense_sector_hamiltonian(h)
        es=eigen(Hermitian(dense));shift=u*(8/4-1/2);physical=es.values.-shift
        # Map production sector basis rows onto physical one-electron spatial-mode order.
        physical_rows=Int[]
        for i in 1:length(h.basis)
            up,down=ED._basis_state(h.basis,i);state=up|down
            count_ones(state)==1||error("Wrong one-electron basis")
            push!(physical_rows,trailing_zeros(state)+1)
        end
        sort(physical_rows)==collect(1:8)||error("Incomplete one-electron basis")
        vectors=zeros(ComplexF64,8,8)
        for (i,j) in enumerate(physical_rows);vectors[j,:]=es.vectors[i,:];end
        native=ED.ground_state_observables(h;method=:dense)
        push!(rows,Dict("U"=>u,"Nup"=>nu,"Ndown"=>nd,"dimension"=>size(dense,1),
            "centered_shift"=>shift,"energies"=>physical,
            "groups"=>eigengroups(physical,vectors,p["degeneracy_tolerance"]),
            "residual"=>maximum(abs,dense*es.vectors-es.vectors*Diagonal(es.values)),
            "orthogonality"=>maximum(abs,adjoint(es.vectors)*es.vectors-I),
            "native_ground_energy_physical"=>native.energy-shift,
            "native_ground_residual"=>native.residual,"native_double_occupancy"=>native.double_occupancy,
            "native_converged"=>native.converged))
    end
    rows
end

function a16(meshes)
    p=M["A16"];parameters=BAND.ModelParameters(p["hoppings"]...);rows=Any[]
    for nk in meshes
        start=time();samples=BAND.sample_band_energies(nk;parameters=parameters)
        raw=BAND.density_of_states(samples;bins=p["bins"],smoothing_width=0.)
        lo,hi=extrema(samples.energies);de=raw.bin_width
        lower=sum(raw.counts[raw.energy.<(-p["edge_window_start"])])/(2nk^2)
        upper=sum(raw.counts[raw.energy.>(p["edge_window_start"])])/(2nk^2)
        curves=Any[]
        for sigma in p["smoothing_widths"]
            dos=BAND.density_of_states(samples;bins=p["bins"],smoothing_width=sigma)
            density=dos.density./2;imax=argmax(density)
            push!(curves,Dict("sigma"=>sigma,"density"=>density,
                "normalization"=>sum(density)*de,"peak"=>density[imax],"peak_energy"=>dos.energy[imax]))
        end
        push!(rows,Dict("Nk"=>nk,"energy"=>raw.energy,"bin_width"=>de,
             "minimum_band_energy"=>lo,"maximum_band_energy"=>hi,"total_counts"=>sum(raw.counts),
             "lower_edge_mass"=>lower,"upper_edge_mass"=>upper,"curves"=>curves,
             "elapsed_seconds"=>time()-start))
    end
    rows
end

function a17()
    p=M["A17"];rows=Any[]
    # Native atom helper has x-up,y-up,x-down,y-down order.
    c=[ATOM._annihilation_operator(i,4) for i in 1:4]
    n=[adjoint(a)*a for a in c];nn=[n[1]+n[3],n[2]+n[4]]
    spin=Any[]
    for a in 1:2
        plus=adjoint(c[a])*c[a+2];minus=adjoint(plus)
        push!(spin,[(plus+minus)/2,(plus-minus)/(2im),(n[a]-n[a+2])/2])
    end
    for u in p["interactions"],mu in p["centered_chemical_potentials"],beta in p["inverse_temperatures"]
        generator=zeros(ComplexF64,16,16);coverage=zeros(Int,16)
        for nu in 0:2,nd in 0:2
            h=ED.sector_hamiltonian(zeros(2,2),u,nu,nd);block=ED.dense_sector_hamiltonian(h)
            indices=Int[]
            for i in 1:length(h.basis)
                up,down=ED._basis_state(h.basis,i);push!(indices,Int(up|(down<<2))+1)
            end
            coverage[indices].+=1
            generator[indices,indices]=block-mu*(nu+nd)*I
        end
        coverage==ones(Int,16)||error("Incomplete grand canonical basis")
        es=eigen(Hermitian(generator));weights=exp.(-beta.*(es.values.-minimum(es.values)));weights./=sum(weights)
        rho=es.vectors*Diagonal(weights)*adjoint(es.vectors)
        mean(op)=real(tr(rho*op))
        densconn=mean(nn[1]*nn[2])-mean(nn[1])*mean(nn[2])
        spinconn=[mean(spin[1][a]*spin[2][b])-mean(spin[1][a])*mean(spin[2][b]) for a in 1:3 for b in 1:3]
        native=ATOM.solve_two_orbital_atom(u,mu,beta,[0.,beta/2,beta])
        push!(rows,Dict("U"=>u,"mu"=>mu,"beta"=>beta,"dimension"=>16,"sector_count"=>9,
             "probabilities"=>real.(diag(rho)),"rho_offdiagonal"=>maximum(abs,rho-Diagonal(diag(rho))),
             "rho_trace"=>real(tr(rho)),"rho_minimum_eigenvalue"=>minimum(eigvals(Hermitian(rho))),
             "residual"=>maximum(abs,generator*es.vectors-es.vectors*Diagonal(es.values)),
             "density"=>mean(nn[1]+nn[2]),"double_occupancy"=>mean(n[1]*n[3]+n[2]*n[4]),
             "density_connected"=>densconn,"spin_connected"=>spinconn,
             "local_spin_squared"=>[mean(s*s) for s in spin[1]],
             "local_density_variance"=>mean(nn[1]*nn[1])-mean(nn[1])^2,
             "S_alm_q0"=>mean((spin[1][3]-spin[2][3])^2),
             "delta_alm"=>mean(spin[1][3]-spin[2][3]),
             "native_density"=>native.density,"native_double_occupancy"=>native.double_occupancy,
             "native_S_alm_q0"=>native.S_alm_q0,"native_delta_alm"=>native.delta_alm))
    end
    rows
end

function main()
    output=ENV["ANALYTIC_TWO_ORBITAL_OUTPUT"];ispath(output)&&error("Refusing overwrite")
    string(VERSION)==M["julia_version"]||error("Unpinned Julia version")
    hashes=Dict(path=>sha(joinpath(SOURCE,path)) for path in keys(M["source_hashes"]))
    hashes==M["source_hashes"]||error("Unpinned source closure")
    mode=get(ENV,"ANALYTIC_TWO_ORBITAL_MODE","full");mode in ["pilot","full"]||error("Invalid mode")
    start=time();r=Dict{String,Any}("release"=>M["release"],"mode"=>mode,
        "manifest_sha256"=>sha(MANIFEST),"driver_sha256"=>sha(@__FILE__),"source_hashes"=>hashes,
        "julia_version"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads(),"julia_threads"=>Threads.nthreads(),
        "host"=>get(ENV,"HOSTNAME","UNKNOWN"),"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","NONE"),
        "created_utc"=>string(now(UTC)))
    # A separate pilot is timing/feasibility only. It cannot complete any catalogue question.
    if mode=="pilot"
        r["A15"]=a15();r["A17"]=a17()
        r["A16"]=a16(M["resources"]["pilot_meshes"])
        r["scientific_completions"]=0
    else
        r["A15"]=a15();r["A16"]=a16(M["A16"]["meshes"]);r["A17"]=a17()
    end
    r["elapsed_seconds"]=time()-start
    open(output,"w") do io;TOML.print(io,r);end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__);main();end

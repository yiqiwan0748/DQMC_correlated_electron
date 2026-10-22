module ALMDMRGBenchmark20260911

using ITensors, ITensorMPS
using LinearAlgebra, Random, Serialization, SHA, TOML, Dates

export main, run_benchmark, preflight, fermionic_local_sign, from_ed, to_ed

const SEEDS = (314159, 271828)
const CAPS = (32, 64, 128, 256)
const NSECTORS = (0, 1)
const ACTION_ATOL = 2e-10
const ENERGY_ATOL = 1e-8
const VARIANCE_ATOL = 1e-8
const SPIN_ATOL = 1e-6
const NUMBER_ATOL = 1e-10
const SOURCE = abspath(@__FILE__)

timestamp() = string(now(UTC)) * "Z"
sha(path) = bytes2hex(sha256(read(path)))

function write_atomic(path, writer)
    ispath(path) && error("Refusing to replace existing artifact: $path")
    tmp = path * ".partial"
    ispath(tmp) && error("Partial artifact already exists: $tmp")
    open(writer, tmp, "w")
    mv(tmp, path)
    return path
end
save_jls(path, x) = write_atomic(path, io -> serialize(io, x))
save_toml(path, x) = write_atomic(path, io -> TOML.print(io, x))

function csvcell(x)
    s = string(x)
    return occursin(r"[,\n\r\"]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end
function append_row(path, names, values)
    fresh = !isfile(path)
    open(path, "a") do io
        fresh && println(io, join(names, ','))
        println(io, join(csvcell.(collect(values)), ','))
        flush(io)
    end
end
function journal(out, event, detail="")
    append_row(joinpath(out, "journal.csv"), ("utc", "event", "detail"),
               (timestamp(), event, detail))
end

"""All-up/all-down Fock basis to interleaved Electron-site basis parity."""
function fermionic_local_sign(up::Integer, down::Integer, m::Int=8)
    inversions = 0
    for i in 0:(m-1)
        ((down >> i) & 1) == 1 && (inversions += count_ones(up >> (i+1)))
    end
    return isodd(inversions) ? -1.0 : 1.0
end

function local_values(up::Integer, down::Integer)
    return ntuple(p -> 1 + Int((up >> (p-1)) & 1) + 2Int((down >> (p-1)) & 1), 8)
end

function mapping(sec)
    labels = NTuple{8,Int}[]
    signs = Float64[]
    lookup = Dict{Tuple{UInt64,UInt64},Int}()
    for up in sec.up_states, down in sec.down_states
        push!(labels, local_values(up, down))
        push!(signs, fermionic_local_sign(up, down))
        lookup[(UInt64(up), UInt64(down))] = length(labels)
    end
    return (; labels, signs, lookup)
end

function initial_labels(Sz::Int, which::Int)
    if Sz == 0
        return which == 1 ? ["Up","Dn","Up","Dn","Up","Dn","Emp","Emp"] :
                            ["UpDn","Emp","Up","Dn","Emp","Up","Dn","Emp"]
    elseif Sz == 1
        return which == 1 ? ["Up","Up","Dn","Up","Up","Dn","Emp","Emp"] :
                            ["UpDn","Emp","Up","Emp","Dn","Up","Up","Emp"]
    end
    error("Only registered Sz=0,1 sectors are allowed")
end

function masks_from_labels(labels)
    up, down = UInt64(0), UInt64(0)
    for (p, label) in enumerate(labels)
        label in ("Up", "UpDn") && (up |= UInt64(1) << (p-1))
        label in ("Dn", "UpDn") && (down |= UInt64(1) << (p-1))
    end
    return up, down
end

"""Lossless fixed-sector ED coefficient tensor to QN-conserving MPS."""
function from_ed(v, sites, sec, map=mapping(sec))
    length(v) == length(map.labels) || error("ED vector length mismatch")
    q = flux(MPS(sites, initial_labels(sec.Sz, 1)))
    tensor = ITensor(Float64, q, sites...)
    for r in eachindex(v)
        pairs = ntuple(p -> sites[p] => map.labels[r][p], 8)
        tensor[pairs...] = map.signs[r] * v[r]
    end
    return MPS(tensor, sites; cutoff=0.0, maxdim=256)
end

"""Contract an MPS and recover ED coefficients with the inverse Fock parity."""
function to_ed(psi, sites, sec, map=mapping(sec))
    tensor = psi[1]
    for p in 2:8
        tensor = tensor * psi[p]
    end
    v = Vector{Float64}(undef, length(map.labels))
    maximag = 0.0
    for r in eachindex(v)
        pairs = ntuple(p -> sites[p] => map.labels[r][p], 8)
        z = tensor[pairs...]
        maximag = max(maximag, abs(imag(z)))
        v[r] = map.signs[r] * real(z)
    end
    maximag <= 1e-12 || error("Unexpected complex coefficient in real benchmark: $maximag")
    return v
end

"""Independent open2x2 hopping assembly from the registered Bloch coefficients."""
function registered_kinetic()
    K = zeros(8, 8)
    t1,t2,t3,t4 = -1.0,-1.75,-0.85,-0.65
    index(x,y,a) = 2(x + 2y) + a
    for y in 0:1, x in 0:1
        for (dx,dy,T) in ((1,0,[-t1 0.0;0.0 -t2]),
                          (0,1,[-t2 0.0;0.0 -t1]),
                          (1,1,[-t3 t4;t4 -t3]),
                          (1,-1,[-t3 -t4;-t4 -t3]))
            xx,yy = x+dx,y+dy
            0 <= xx <= 1 && 0 <= yy <= 1 || continue
            for a in 1:2, b in 1:2
                p,q = index(x,y,a),index(xx,yy,b)
                K[p,q] += T[a,b]
                K[q,p] += T[a,b]
            end
        end
    end
    return K
end

function build_mpo(K, U, sites)
    os = OpSum()
    for p in 1:8
        !iszero(K[p,p]) && add!(os, K[p,p], "Ntot", p)
        os += U, "Nupdn", p
        os += -U/2, "Ntot", p
        os += U/4, "Id", p
    end
    for p in 1:7, q in (p+1):8
        iszero(K[p,q]) && continue
        for (create,destroy) in (("Cdagup","Cup"),("Cdagdn","Cdn"))
            os += K[p,q], create,p,destroy,q
            os += K[q,p], create,q,destroy,p
        end
    end
    return MPO(Float64, os, sites)
end

function validate_reference(ref)
    ref.cells == 4 && ref.N == 6 && ref.U == 20.0 && ref.t4 == -0.65 ||
        error("Reference is outside the registered model")
    size(ref.kinetic) == (8,8) || error("Expected exactly eight orbital modes")
    maximum(abs, ref.kinetic - registered_kinetic()) <= 2e-12 ||
        error("Reference K disagrees with independently assembled open target hopping")
    Set(keys(ref.sectors)) == Set(NSECTORS) || error("Only Sz0 and Sz1 are registered")
    for Sz in NSECTORS
        sec = ref.sectors[Sz]
        sec.Sz == Sz && sec.Nup == 3+Sz && sec.Ndown == 3-Sz ||
            error("Sector particle counts do not match")
        D = length(sec.up_states)*length(sec.down_states)
        size(sec.H) == (D,D) || error("Dense H dimensions do not match masks")
        maximum(abs, sec.H-transpose(sec.H)) <= 2e-12 || error("Reference H is not real symmetric")
        size(sec.eigenvectors) == (D,sec.ground_degeneracy) || error("Ground basis shape mismatch")
        maximum(abs, sec.eigenvectors' * sec.eigenvectors-I) <= 2e-9 ||
            error("Ground basis is not orthonormal")
        maximum(abs, sec.H*sec.eigenvectors-sec.energy*sec.eigenvectors) <= 2e-9 ||
            error("Reference ground eigenvectors fail their residual")
        length(unique(sec.up_states)) == length(sec.up_states) || error("Repeated up masks")
        length(unique(sec.down_states)) == length(sec.down_states) || error("Repeated down masks")
        all(u -> count_ones(u) == sec.Nup, sec.up_states) || error("Wrong up masks")
        all(d -> count_ones(d) == sec.Ndown, sec.down_states) || error("Wrong down masks")
        size(sec.spin_dot) == (8,8) || error("Expected full physical spin-dot matrix")
    end
    return true
end

function negative_hop_fixture(K, sec, spin::Symbol, map)
    for up0 in sec.up_states, down0 in sec.down_states
        up,down = UInt64(up0),UInt64(down0)
        mask = spin == :up ? up : down
        for p in 1:8, q in 1:8
            abs(p-q) > 1 && !iszero(K[p,q]) || continue
            ((mask >> (q-1)) & 1) == 1 && ((mask >> (p-1)) & 1) == 0 || continue
            between = count(k -> ((mask >> (k-1)) & 1) == 1, (min(p,q)+1):(max(p,q)-1))
            isodd(between) || continue
            newmask = xor(mask, (UInt64(1) << (q-1)) | (UInt64(1) << (p-1)))
            dst = spin == :up ? (newmask,down) : (up,newmask)
            srcindex,dstindex = map.lookup[(up,down)],map.lookup[dst]
            return (; srcindex,dstindex,p,q,expected=-K[p,q])
        end
    end
    error("No nonlocal odd-parity hopping fixture found")
end

function preflight(ref, H, sites, out)
    validate_reference(ref)
    rows = NamedTuple[]
    for Sz in NSECTORS
        sec = ref.sectors[Sz]
        map = mapping(sec)
        D = size(sec.H,1)
        fixtures = Pair{String,Vector{Float64}}[]
        rng = Xoshiro(94871 + Sz)
        for j in 1:3
            v = randn(rng,D); normalize!(v)
            all(x -> !iszero(x),v) || error("Random fixture lacks full support")
            push!(fixtures,"full_support_$j"=>v)
        end
        for j in 1:2
            v = zeros(D)
            v[map.lookup[masks_from_labels(initial_labels(Sz,j))]] = 1.0
            push!(fixtures,"product_$j"=>v)
        end
        for spin in (:up,:down)
            fixture = negative_hop_fixture(ref.kinetic,sec,spin,map)
            entryerr = abs(sec.H[fixture.dstindex,fixture.srcindex]-fixture.expected)
            append_row(joinpath(out,"fermion_sign_checks.csv"),
                ("Sz","spin","p","q","src","dst","expected","ED_entry","error"),
                (Sz,spin,fixture.p,fixture.q,fixture.srcindex,fixture.dstindex,
                 fixture.expected,sec.H[fixture.dstindex,fixture.srcindex],entryerr))
            entryerr <= 2e-12 || error("Independent nonlocal fermion sign test failed")
            v=zeros(D);v[fixture.srcindex]=1.0
            push!(fixtures,"negative_nonlocal_$(spin)"=>v)
        end
        for (label,v) in fixtures
            psi = from_ed(v,sites,sec,map)
            back = to_ed(psi,sites,sec,map)
            roundtrip = maximum(abs,back-v)
            # Direct SVD contraction avoids density-matrix squaring in this lossless check.
            hp = apply(H,psi;alg="naive",cutoff=0.0,maxdim=256)
            got = to_ed(hp,sites,sec,map)
            expected = sec.H*v
            action = maximum(abs,got-expected)
            relative = norm(got-expected)/max(1.0,norm(expected))
            leakage = abs(real(inner(hp,hp))-sum(abs2,got))
            passed = roundtrip <= ACTION_ATOL && action <= ACTION_ATOL &&
                     relative <= 2e-11 && leakage <= 1e-8
            row = (;Sz,label,dimension=D,roundtrip,action,relative,leakage,passed)
            push!(rows,row)
            append_row(joinpath(out,"preflight.csv"), keys(row),values(row))
            passed || error("MPO/mapping preflight failed for Sz=$Sz $label: $row")
        end
    end
    save_toml(joinpath(out,"preflight_summary.toml"),Dict(
        "status"=>"PASS","action_fixtures"=>length(rows),"full_support_per_sector"=>3,
        "random_seeds"=>[94871,94872],"max_action_error"=>maximum(r.action for r in rows),
        "max_roundtrip_error"=>maximum(r.roundtrip for r in rows),
        "action_tolerance"=>ACTION_ATOL,"kinetic_entry_tolerance"=>2e-12,
        "algorithm"=>"naive MPO-MPS contraction, cutoff0,maxdim256; full QN tensor mapping"))
    return rows
end

function measure_state(psi, H, sites, sec, map)
    normalize!(psi)
    norm2=real(inner(psi,psi))
    E=real(inner(psi',H,psi))
    H2=real(inner(H,psi,H,psi))
    variance=H2-E^2
    v=to_ed(psi,sites,sec,map)
    vnorm=sum(abs2,v)
    Hv=sec.H*v
    fock_energy=dot(v,Hv)/vnorm
    fock_variance=sum(abs2,Hv-fock_energy*v)/vnorm
    fidelity=sum(abs2,sec.eigenvectors'*v)/vnorm
    Czz=correlation_matrix(psi,"Sz","Sz";ishermitian=false)
    Cpm=correlation_matrix(psi,"S+","S-";ishermitian=false)
    Cmp=correlation_matrix(psi,"S-","S+";ishermitian=false)
    Cdot=Czz+0.5*(Cpm+Cmp)
    imagerr=maximum(abs,imag.(Cdot))
    imagerr<=1e-10 || error("Complex spin scalar in real state")
    Cdot=real.(Cdot);Czz=real.(Czz)
    nup=real.(expect(psi,"Nup"));ndown=real.(expect(psi,"Ndn"))
    double=real.(expect(psi,"Nupdn"));density=nup+ndown
    spinz=real.(expect(psi,"Sz"))
    odd=[isodd(p) ? 1.0 : -1.0 for p in 1:8]
    S2=sum(Cdot)
    spinerr=maximum(abs,Cdot-sec.spin_dot)
    contacterr=maximum(abs,diag(Cdot)-0.75*(density-2double))
    particleerr=max(abs(sum(nup)-sec.Nup),abs(sum(ndown)-sec.Ndown))
    szerr=abs(sum(spinz)-sec.Sz)
    energyerr=abs(E-sec.energy)
    spin2err=abs(S2-sec.total_spin_squared)
    edactionerr=abs(E-fock_energy)
    leakage=abs(norm2-vnorm)
    # A degenerate reference is an ensemble; do not call its selected-vector mismatch bias.
    ensemble_sensitive=sec.ground_degeneracy>1
    energy_pass=energyerr<=ENERGY_ATOL
    variance_pass=(-VARIANCE_ATOL<=variance<=VARIANCE_ATOL && fock_variance<=VARIANCE_ATOL)
    spin_pass=spinerr<=SPIN_ATOL && spin2err<=SPIN_ATOL
    number_pass=particleerr<=NUMBER_ATOL && szerr<=NUMBER_ATOL
    internal_pass=edactionerr<=ENERGY_ATOL && contacterr<=SPIN_ATOL && leakage<=NUMBER_ATOL
    passed=energy_pass && variance_pass && spin_pass && number_pass && internal_pass
    return (;E,H2,variance,fock_energy,fock_variance,fidelity,norm2,leakage,Cdot,Czz,
        density,double,nup,ndown,spinz,S2,spinerr,contacterr,particleerr,szerr,energyerr,
        spin2err,edactionerr,ensemble_sensitive,energy_pass,variance_pass,spin_pass,
        number_pass,internal_pass,passed,Sodd_z=dot(odd,Czz*odd)/4,
        Sodd_scalar=dot(odd,Cdot*odd)/12,Delta=dot(odd,spinz)/4,
        densityerr=maximum(abs,density-sec.density),
        doubleerr=maximum(abs,double-sec.double_occupancy))
end

function max_logged_truncation(path)
    vals=[parse(Float64,m.captures[1]) for m in eachmatch(r"maxerr\s*=\s*([0-9.eE+\-]+)",read(path,String))]
    return isempty(vals) ? NaN : maximum(vals)
end

function save_stage(out,Sz,seed,stage,cap,reported,seconds,maxtrunc,psi,m)
    tag="Sz$(Sz)_seed$(seed)_stage$(stage)_cap$(cap)"
    row=(;Sz,seed,stage,cap,energy=m.E,reported_energy=reported,energy_error=m.energyerr,
        variance=m.variance,fock_variance=m.fock_variance,fidelity=m.fidelity,
        spin_matrix_error=m.spinerr,total_spin_squared=m.S2,spin_squared_error=m.spin2err,
        Sodd_z=m.Sodd_z,Sodd_scalar=m.Sodd_scalar,Delta=m.Delta,
        particle_error=m.particleerr,sz_error=m.szerr,contact_error=m.contacterr,
        density_error=m.densityerr,double_occupancy_error=m.doubleerr,
        energy_action_error=m.edactionerr,sector_leakage=m.leakage,
        maxlinkdim=maxlinkdim(psi),max_truncation=maxtrunc,seconds,
        ensemble_sensitive=m.ensemble_sensitive,energy_pass=m.energy_pass,
        variance_pass=m.variance_pass,spin_pass=m.spin_pass,number_pass=m.number_pass,
        internal_pass=m.internal_pass,passed=m.passed)
    append_row(joinpath(out,"stages.csv"),keys(row),values(row))
    for p in 1:8,q in 1:8
        append_row(joinpath(out,"spin_matrices.csv"),
            ("Sz","seed","stage","cap","p","q","spin_dot","SzSz"),
            (Sz,seed,stage,cap,p,q,m.Cdot[p,q],m.Czz[p,q]))
    end
    for p in 1:8
        append_row(joinpath(out,"local_observables.csv"),
            ("Sz","seed","stage","cap","mode","Nup","Ndown","density","double_occupancy","Sz_local"),
            (Sz,seed,stage,cap,p,m.nup[p],m.ndown[p],m.density[p],m.double[p],m.spinz[p]))
    end
    save_jls(joinpath(out,"checkpoints",tag*".jls"),(;Sz,seed,stage,cap,psi=deepcopy(psi),measurement=m))
    return row
end

function metadata(refpath,mode)
    md=Dict{String,Any}("created_utc"=>timestamp(),"mode"=>mode,"reference_path"=>abspath(refpath),
        "reference_sha256"=>sha(refpath),"runner_sha256"=>sha(SOURCE),"julia_version"=>string(VERSION),
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "ITensors_version"=>string(Base.pkgversion(ITensors)),
        "ITensorMPS_version"=>string(Base.pkgversion(ITensorMPS)),
        "caps"=>collect(CAPS),"seeds"=>collect(SEEDS),"physical_Sz"=>collect(NSECTORS),
        "sweeps_per_cap"=>6,"extra_zero_noise_sweeps"=>6,"cutoff"=>1e-12,
        "eigsolve_tol"=>1e-11,"eigsolve_krylovdim"=>10,"eigsolve_maxiter"=>8,
        "first_cap_noise"=>[1e-6,1e-7,1e-8,0.0,0.0,0.0],
        "later_cap_noise"=>[1e-8,1e-9,0.0,0.0,0.0,0.0],
        "energy_tolerance"=>ENERGY_ATOL,"variance_tolerance"=>VARIANCE_ATOL,
        "spin_entry_tolerance"=>SPIN_ATOL,"particle_Sz_tolerance"=>NUMBER_ATOL,
        "number_convention"=>"Electron sites conserve Nf and 2*physicalSz; no total-S or orbital-population constraint",
        "scalar_convention"=>"spin_dot contains all three components; Sodd_scalar divides sum by3*Ncell",
        "scope"=>"open2x2 U20 N6 t4=-.65; exactly four arms; finite-cluster validation")
    for (name,mod) in (("ITensors",ITensors),("ITensorMPS",ITensorMPS))
        path=pathof(mod)
        md[name*"_entry_source_sha256"]=sha(path)
        md[name*"_entry_source_path"]=path
    end
    for name in ("Project.toml","Manifest.toml","setup.jl")
        path=joinpath(@__DIR__,name)
        isfile(path) && (md[name*"_sha256"]=sha(path))
    end
    return md
end

function run_benchmark(refpath::AbstractString,out::AbstractString;preflight_only::Bool=false)
    ispath(out) && error("Refusing existing scientific output directory: $out")
    isfile(refpath) || error("Reference file does not exist")
    mkpath(out);mkpath(joinpath(out,"checkpoints"));mkpath(joinpath(out,"logs"))
    active_state=nothing
    try
        BLAS.set_num_threads(1)
        Threads.nthreads()==1 || error("Registered benchmark requires one Julia thread")
        Base.pkgversion(ITensors)==v"0.9.30" || error("ITensors pin mismatch")
        Base.pkgversion(ITensorMPS)==v"0.4.1" || error("ITensorMPS pin mismatch")
        mode=preflight_only ? "preflight-only" : "run"
        save_toml(joinpath(out,"metadata.toml"),metadata(refpath,mode))
        cp(SOURCE,joinpath(out,"run_benchmark.jl"))
        for file in ("Project.toml","Manifest.toml")
            p=joinpath(@__DIR__,file);isfile(p) && cp(p,joinpath(out,file))
        end
        journal(out,"reference_load")
        ref=deserialize(refpath)
        validate_reference(ref)
        sites=siteinds("Electron",8;conserve_nf=true,conserve_sz=true)
        H=build_mpo(ref.kinetic,ref.U,sites)
        journal(out,"preflight_start")
        preflight(ref,H,sites,out)
        journal(out,"preflight_pass")
        if preflight_only
            write_atomic(joinpath(out,"status.txt"),io->println(io,"PREFLIGHT_PASS"))
            return 0
        end
        finalstates=Dict{Tuple{Int,Int},Any}()
        finalrows=NamedTuple[]
        finalmeasurements=Dict{Tuple{Int,Int},Any}()
        for Sz in NSECTORS
            sec=ref.sectors[Sz];map=mapping(sec)
            for (initialization,seed) in enumerate(SEEDS)
                labels=initial_labels(Sz,initialization)
                psi=random_mps(Xoshiro(seed),sites,labels;linkdims=4)
                active_state=(;Sz,seed,stage=0,psi)
                save_jls(joinpath(out,"checkpoints","Sz$(Sz)_seed$(seed)_initial.jls"),
                    (;Sz,seed,labels,psi=deepcopy(psi),flux=string(flux(psi))))
                journal(out,"arm_start","Sz=$Sz seed=$seed")
                for stage in 1:5
                    cap=stage<=4 ? CAPS[stage] : 256
                    noise=stage==1 ? [1e-6,1e-7,1e-8,0.0,0.0,0.0] :
                          stage==5 ? zeros(6) : [1e-8,1e-9,0.0,0.0,0.0,0.0]
                    logpath=joinpath(out,"logs","Sz$(Sz)_seed$(seed)_stage$(stage).log")
                    journal(out,"stage_start","Sz=$Sz seed=$seed stage=$stage cap=$cap")
                    t0=time_ns()
                    result=open(logpath,"w") do io
                        redirect_stdout(io) do
                            redirect_stderr(io) do
                                dmrg(H,psi;nsweeps=6,maxdim=fill(cap,6),cutoff=1e-12,
                                    noise,eigsolve_tol=1e-11,eigsolve_krylovdim=10,
                                    eigsolve_maxiter=8,outputlevel=1)
                            end
                        end
                    end
                    reported,psi=result
                    active_state=(;Sz,seed,stage,cap,psi)
                    seconds=(time_ns()-t0)/1e9
                    m=measure_state(psi,H,sites,sec,map)
                    row=save_stage(out,Sz,seed,stage,cap,reported,seconds,
                                   max_logged_truncation(logpath),psi,m)
                    journal(out,"stage_complete","Sz=$Sz seed=$seed stage=$stage passed=$(m.passed)")
                    if stage==5
                        push!(finalrows,row)
                        finalstates[(Sz,seed)]=deepcopy(psi)
                        finalmeasurements[(Sz,seed)]=m
                    end
                end
            end
        end
        length(finalstates)==4 || error("Did not complete exactly four registered arms")
        agreement_pass=true
        for Sz in NSECTORS
            a,b=(finalmeasurements[(Sz,s)] for s in SEEDS)
            dE=abs(a.E-b.E);dC=maximum(abs,a.Cdot-b.Cdot);dS=abs(a.S2-b.S2)
            ok=dE<=2ENERGY_ATOL && dC<=2SPIN_ATOL && dS<=2SPIN_ATOL
            agreement_pass &= ok
            append_row(joinpath(out,"initialization_agreement.csv"),
                ("Sz","energy_difference","spin_matrix_difference","S2_difference","within_twice_reference_tolerance"),
                (Sz,dE,dC,dS,ok))
        end
        save_jls(joinpath(out,"finalstates.jls"),(;sites,states=finalstates,measurements=finalmeasurements))
        degeneracy_sensitive=any(r.ensemble_sensitive for r in finalrows)
        passed=all(r.passed for r in finalrows) && agreement_pass
        status=passed ? "PASS" : degeneracy_sensitive ? "INCONCLUSIVE_ENSEMBLE_SENSITIVE" : "FAIL_REGISTERED_GATES"
        save_toml(joinpath(out,"summary.toml"),Dict(
            "status"=>status,"completed_arms"=>4,"completed_stages"=>20,
            "initialization_agreement_pass"=>agreement_pass,
            "ensemble_sensitive"=>degeneracy_sensitive,
            "max_final_energy_error"=>maximum(r.energy_error for r in finalrows),
            "max_final_variance"=>maximum(r.variance for r in finalrows),
            "max_final_fock_variance"=>maximum(r.fock_variance for r in finalrows),
            "max_final_spin_matrix_error"=>maximum(r.spin_matrix_error for r in finalrows),
            "max_final_S2_error"=>maximum(r.spin_squared_error for r in finalrows),
            "minimum_final_ground_subspace_fidelity"=>minimum(r.fidelity for r in finalrows),
            "finished_utc"=>timestamp(),"claim"=>"finite-cluster method benchmark only"))
        write_atomic(joinpath(out,"status.txt"),io->println(io,status))
        journal(out,"complete",status)
        return passed ? 0 : 2
    catch err
        bt=catch_backtrace()
        open(joinpath(out,"failure.txt"),"w") do io
            println(io,timestamp());showerror(io,err,bt);println(io)
        end
        if active_state!==nothing
            try
                save_jls(joinpath(out,"failure_state.jls"),active_state)
            catch saveerr
                open(joinpath(out,"failure.txt"),"a") do io
                    println(io,"Saving failure state also failed: ",sprint(showerror,saveerr))
                end
            end
        end
        open(joinpath(out,"status.txt"),"w") do io;println(io,"ERROR_PRESERVED");end
        journal(out,"error",sprint(showerror,err))
        rethrow()
    end
end

function main(args=ARGS)
    options=Dict{String,String}();mode=nothing;i=1
    while i<=length(args)
        arg=args[i]
        if arg in ("--run","--preflight-only")
            mode===nothing || error("Choose exactly one execution mode")
            mode=arg;i+=1
        elseif arg in ("--reference","--output")
            i<length(args) || error("Missing value for $arg")
            haskey(options,arg) && error("Duplicate $arg")
            options[arg]=args[i+1];i+=2
        elseif arg=="--help"
            println("run_benchmark.jl (--preflight-only|--run) --reference /path/reference.jls --output /new/output")
            return 0
        else
            error("Unknown argument $arg; no implicit model or run overrides are allowed")
        end
    end
    mode!==nothing && haskey(options,"--reference") && haskey(options,"--output") ||
        error("Explicit mode, --reference and --output are required")
    return run_benchmark(abspath(options["--reference"]),abspath(options["--output"]);
                         preflight_only=mode=="--preflight-only")
end

end # module

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    exit(ALMDMRGBenchmark20260911.main())
end

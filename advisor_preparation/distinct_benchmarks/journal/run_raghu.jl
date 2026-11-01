#!/usr/bin/env julia
using LinearAlgebra, Statistics, TOML, SHA, Dates
BLAS.set_num_threads(1)
const ALM_ROOT=get(ENV,"DISTINCT_ALM_ROOT",normpath(joinpath(@__DIR__,"../frozen_alm")))
module Target
include(joinpath(Main.ALM_ROOT,"src/noninteracting.jl"))
include(joinpath(Main.ALM_ROOT,"src/realspace_kinetic.jl"))
end
const T=(-1.0,1.3,-.85,-.85)
const MU=1.45
const PARAMS=Target.ModelParameters(T...)
const SOURCE="https://arxiv.org/pdf/0804.1113v1"
# Independent evaluation of the author's Eqs. (2)-(4), without target helpers.
function published(kx,ky)
    t1,t2,t3,t4=T
    ex=-2t1*cos(kx)-2t2*cos(ky)-4t3*cos(kx)*cos(ky)
    ey=-2t2*cos(kx)-2t1*cos(ky)-4t3*cos(kx)*cos(ky)
    xy=-4t4*sin(kx)*sin(ky)
    center=(ex+ey)/2;z=(ex-ey)/2;r=sqrt(z*z+xy*xy)
    (;energies=[center-r,center+r],z,x=xy,r)
end
wrap(x)=mod(x+pi,2pi)-pi
distance(a,b)=hypot(wrap(a[1]-b[1]),wrap(a[2]-b[2]))
function components(mask,ks)
    N=length(ks);seen=falses(N,N);out=Any[];queue=Vector{Int}(undef,N*N)
    for start in eachindex(mask)
        mask[start]&&!seen[start]||continue
        first=1;last=1;queue[1]=start;seen[start]=true;sx=0.;cx=0.;sy=0.;cy=0.
        while first<=last
            q=queue[first];first+=1;ix=mod(q-1,N)+1;iy=(q-1)÷N+1
            sx+=sin(ks[ix]);cx+=cos(ks[ix]);sy+=sin(ks[iy]);cy+=cos(ks[iy])
            for (nx,ny) in ((mod1(ix-1,N),iy),(mod1(ix+1,N),iy),(ix,mod1(iy-1,N)),(ix,mod1(iy+1,N)))
                j=nx+(ny-1)*N
                if mask[j]&&!seen[j];last+=1;queue[last]=j;seen[j]=true;end
            end
        end
        push!(out,Dict("cells"=>last,"center"=>[atan(sx,cx),atan(sy,cy)]))
    end
    out
end
function crossing(energy,center,angle,band)
    direction=[cos(angle),sin(angle)]
    f(r)=energy(center[1]+r*direction[1],center[2]+r*direction[2])[band]-MU
    a=0.;fa=f(a)
    for b in range(.01,1.55;length=256)
        fb=f(b)
        if signbit(fa)!=signbit(fb)
            for _ in 1:55
                m=(a+b)/2;fm=f(m)
                signbit(fm)==signbit(fa) ? (a=m;fa=fm) : (b=m)
            end
            return (a+b)/2
        end
        a=b;fa=fb
    end
    error("No registered radial Fermi crossing at center=$center angle=$angle")
end
function cross_realspace()
    L=8;K=Target.realspace_kinetic_matrix(L,PARAMS);err=0.
    for nx in 0:L-1,ny in 0:L-1
        kx=2pi*nx/L;ky=2pi*ny/L
        block=Target.fourier_block(K,L,kx,ky)
        q=published(kx,ky)
        reference=[(q.energies[1]+q.energies[2])/2+q.z q.x;q.x (q.energies[1]+q.energies[2])/2-q.z]
        err=max(err,maximum(abs,block-reference))
    end
    err
end
function j20()
    start=time();rows=Any[];fillings=Float64[];maxenergy=0.
    for N in [128,256,512,1024]
        s=Target.sample_band_energies(N;parameters=PARAMS);lower=reshape(s.energies[1,:],N,N);upper=reshape(s.energies[2,:],N,N)
        holes=components(lower.>MU,s.k_values);electrons=components(upper.<MU,s.k_values)
        reference_occupied=0;err=0.
        for iy in 1:N,ix in 1:N
            e=published(s.k_values[ix],s.k_values[iy]).energies
            reference_occupied+=count(<(MU),e)
            err=max(err,abs(lower[ix,iy]-e[1]),abs(upper[ix,iy]-e[2]))
        end
        filling=Target.zero_temperature_filling(s,MU);reference=2reference_occupied/N^2
        push!(fillings,filling);maxenergy=max(maxenergy,err)
        hc=length(holes)==2&&all(any(distance(c["center"],p)<=2pi/N for c in holes) for p in ([0.,0.],[pi,pi]))
        ec=length(electrons)==2&&all(any(distance(c["center"],p)<=2pi/N for c in electrons) for p in ([pi,0.],[0.,pi]))
        push!(rows,Dict("Nk"=>N,"hole_components"=>holes,"electron_components"=>electrons,"published_topology_matches"=>hc&&ec,"target_filling_per_cell_both_spins"=>filling,"independent_equation_filling"=>reference,"filling_error"=>abs(filling-reference),"maximum_band_error"=>err))
    end
    roots=Any[];maxroot=0.
    targetenergy(kx,ky)=eigvals(Hermitian(Target.bloch_hamiltonian(kx,ky,PARAMS)))
    for (center,band) in [([0.,0.],1),([pi,pi],1),([pi,0.],2),([0.,pi],2)],theta in range(0,2pi;length=17)[1:16]
        rt=crossing(targetenergy,center,theta,band);ro=crossing((x,y)->published(x,y).energies,center,theta,band);maxroot=max(maxroot,abs(rt-ro))
        push!(roots,Dict("center"=>center,"band"=>band,"angle"=>theta,"target_radius"=>rt,"published_equation_radius"=>ro,"error"=>abs(rt-ro)))
    end
    realerr=cross_realspace();convergence=abs(fillings[end]-fillings[end-1])
    ok=all(r["published_topology_matches"]&&r["filling_error"]<=1e-12 for r in rows)&&maxenergy<=1e-10&&maxroot<=1e-9&&realerr<=1e-10&&convergence<=1e-3
    Dict("id"=>"J20","status"=>ok ? "PASS" : "FAIL","scope"=>"Author-preprint-v1 Fig3 four-pocket topology and Eqs2-4; actual production band/grid/realspace code, not interacting DQMC","meshes"=>rows,"fermi_rays"=>roots,"maximum_contour_radius_error"=>maxroot,"realspace_Fourier_block_error"=>realerr,"last_mesh_filling_change"=>convergence,"runtime_seconds"=>time()-start)
end
function j21()
    start=time();rows=Any[];maxerr=0.;maxsum=0.;minsplit=Inf
    for ky in [.27pi,.43pi],kx in range(.08pi,.92pi;length=129)
        p=published(kx,ky);p.r>1e-8||error("Degenerate band cannot have gauge-unique orbital residues")
        es=eigen(Hermitian(Target.bloch_hamiltonian(kx,ky,PARAMS)))
        for band in 1:2
            sign=band==1 ? -1 : 1;reference=[(1+sign*p.z/p.r)/2,(1-sign*p.z/p.r)/2]
            target=abs2.(es.vectors[:,band]);err=maximum(abs,target-reference);maxerr=max(maxerr,err);maxsum=max(maxsum,abs(sum(target)-1));minsplit=min(minsplit,2p.r)
            push!(rows,Dict("kx"=>kx,"ky"=>ky,"band"=>band,"target_orbital_weights"=>target,"published_projector_weights"=>reference,"error"=>err))
        end
    end
    # Opposite sides of kx=ky exchange which orbital dominates each band.
    transfers=Bool[]
    for ky in [.27pi,.43pi],band in 1:2
        left=eigen(Hermitian(Target.bloch_hamiltonian(ky-.15pi,ky,PARAMS))).vectors
        right=eigen(Hermitian(Target.bloch_hamiltonian(ky+.15pi,ky,PARAMS))).vectors
        push!(transfers,(abs2(left[1,band])-.5)*(abs2(right[1,band])-.5)<0)
    end
    ok=maxerr<=1e-8&&maxsum<=1e-12&&all(transfers)
    Dict("id"=>"J21","status"=>ok ? "PASS" : "FAIL","scope"=>"Orbital-resolved band residues derived independently from published two-band projector, not a digitized figure or interacting spectral function","maximum_weight_error"=>maxerr,"maximum_weight_sum_error"=>maxsum,"minimum_band_splitting"=>minsplit,"all_registered_orbital_transfers"=>all(transfers),"rows"=>rows,"runtime_seconds"=>time()-start)
end
function main()
    output=ENV["DISTINCT_JOURNAL_OUTPUT"];isfile(output)&&error("Choose a new output path")
    r=Dict{String,Any}("created_utc"=>string(now(UTC)),"julia_version"=>string(VERSION),"host"=>get(ENV,"HOSTNAME","UNKNOWN"),"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"reference_url"=>SOURCE,"reference_version"=>"0804.1113v1; APS final PDF unavailable during preparation; preprint parameters explicitly pinned","reference_doi"=>"10.1103/PhysRevB.77.220503","parameters"=>Dict("t1"=>T[1],"t2"=>T[2],"t3"=>T[3],"t4"=>T[4],"mu"=>MU),"driver_sha256"=>bytes2hex(sha256(read(@__FILE__))),"source_hashes"=>Dict(p=>bytes2hex(sha256(read(joinpath(ALM_ROOT,p)))) for p in ["src/noninteracting.jl","src/realspace_kinetic.jl"]),"benchmarks"=>Any[])
    for f in [j20,j21]
        row=try f() catch ex;Dict("id"=>f===j20 ? "J20" : "J21","status"=>"UNRESOLVED","reason"=>sprint(showerror,ex,catch_backtrace())) end
        push!(r["benchmarks"],row);println(row["id"]," ",row["status"]);flush(stdout)
        open(output,"w") do io;TOML.print(io,r);end
    end
    r["finished_utc"]=string(now(UTC));open(output,"w") do io;TOML.print(io,r);end
    all(b["status"]=="PASS" for b in r["benchmarks"])||exit(1)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__);main();end

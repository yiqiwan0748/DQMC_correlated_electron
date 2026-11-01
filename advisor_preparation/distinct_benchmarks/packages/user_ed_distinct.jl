#!/usr/bin/env julia
using LinearAlgebra, Random, Statistics, TOML, SHA, Dates
const ROOT=ENV["RENYI_ROOT"]
module ED
include(joinpath(Main.ROOT,"analysis/compute_s2_ed.jl"))
end
BLAS.set_num_threads(1)
const M=TOML.parsefile(joinpath(@__DIR__,"manifest.toml"))
const L=M["L"]
rowflat(A)=[Float64(real(A[i,j])) for i in axes(A,1) for j in axes(A,2)]
function applyops(mask,ops)
    sign=1.
    for (kind,p) in ops
        present=!iszero(mask & (1<<p))
        (kind=='-' ? !present : present) && return mask,0.
        isodd(count_ones(mask & ((1<<p)-1))) && (sign=-sign)
        mask=xor(mask,1<<p)
    end
    mask,sign
end
function action(psi,masks,terms)
    out=zeros(1<<(2L))
    for (k,mask) in enumerate(masks),(coefficient,ops) in terms
        target,sign=applyops(mask,ops);out[target+1]+=coefficient*sign*psi[k]
    end
    out
end
function operator_matrix(masks,terms)
    index=Dict(s=>i for (i,s) in enumerate(masks));A=zeros(length(masks),length(masks))
    for (k,mask) in enumerate(masks),(coefficient,ops) in terms
        target,sign=applyops(mask,ops)
        sign!=0 && haskey(index,target) && (A[index[target],k]+=coefficient*sign)
    end
    A
end
function rdm_adapter(psi,space,sites)
    complement=setdiff(collect(0:L-1),sites);nA=length(sites);nB=L-nA;nd=length(space.dn_masks);X=zeros(1<<(2nA),1<<(2nB))
    for (iu,u) in enumerate(space.up_masks),(id,d) in enumerate(space.dn_masks)
        ua=ED.compressed_bits(u,sites);ub=ED.compressed_bits(u,complement);da=ED.compressed_bits(d,sites);db=ED.compressed_bits(d,complement)
        a=Int(ua | (da<<nA))+1;b=Int(ub | (db<<nB))+1;X[a,b]+=ED.subsystem_sign(ua,ub,da)*psi[(iu-1)*nd+id]
    end
    X*X'
end
function projections(C)
    f=eigen(Symmetric(C));w=f.values;V=f.vectors;out=collect(w);j=1
    while j<=L
        k=j+1
        while k<=L && abs(w[k]-w[j])<1e-7;k+=1;end
        append!(out,rowflat(V[:,j:k-1]*V[:,j:k-1]'));j=k
    end
    out
end
function spin_ladder()
    masks=[m for m in 0:(1<<(2L))-1 if count_ones(m)==L];index=Dict(s=>i for (i,s) in enumerate(masks));H=zeros(length(masks),length(masks))
    for nu in 0:L
        space=ED.build_space(L,M["t"],M["U"],0.,nu,L-nu;boundary=M["boundary"]);states=Int[]
        for u in space.up_masks,d in space.dn_masks;push!(states,Int(u)|(Int(d)<<L));end
        n=ED.dimension(space);x=zeros(n);y=zeros(n)
        for j in 1:n
            fill!(x,0.);x[j]=1;ED.apply_hamiltonian!(y,x,space)
            for i in 1:n;H[index[states[i]],index[states[j]]]=y[i];end
        end
    end
    Sp=operator_matrix(masks,[(1.,[('-',L+i),('+',i)]) for i in 0:L-1]);sz=[(count_ones(m & ((1<<L)-1))-count_ones(m>>L))/2 for m in masks];S2=Diagonal(sz.^2)+(Sp*Sp'+Sp'*Sp)/2
    norm(H*S2-S2*H)<=1e-9 || error("Target fullS2 does not commute with user's ED Hamiltonian")
    f=eigen(Symmetric(H));e=f.values;v=f.vectors;out=Float64[];j=1
    while j<=length(e)
        k=j+1
        while k<=length(e)&&abs(e[k]-e[j])<1e-8;k+=1;end
        if e[j]-e[1]<=M["spin_window_above_ground"]+1e-9
            spins=eigvals(Symmetric(v[:,j:k-1]'*S2*v[:,j:k-1]))
            for spin in spins;append!(out,[mean(e[j:k-1])-e[1],spin]);end
        end
        j=k
    end
    out
end
function main()
    output=ENV["DISTINCT_TARGET_OUTPUT"];isfile(output)&&error("Choose a new target output")
    space=ED.build_space(L,M["t"],M["U"],0.,M["Nup"],M["Ndown"];boundary=M["boundary"])
    E,psi,it,reported=ED.ground_state_lanczos(space;tol=1e-12,maxiter=1000);hv=similar(psi);ED.apply_hamiltonian!(hv,psi,space);res=norm(hv-E*psi);res<=M["ground_residual_limit"]||error("User ED ground residual")
    masks=Int[]
    for u in space.up_masks,d in space.dn_masks;push!(masks,Int(u)|(Int(d)<<L));end
    full=zeros(1<<(2L));for (k,s) in enumerate(masks);full[s+1]=psi[k];end
    expect(terms)=dot(full,action(psi,masks,terms))
    Cu=[expect([(1.,[('-',j),('+',i)])]) for i in 0:L-1,j in 0:L-1];Cd=[expect([(1.,[('-',L+j),('+',L+i)])]) for i in 0:L-1,j in 0:L-1]
    occupations=[Float64((s>>i)&1)+Float64((s>>(L+i))&1) for i in 0:L-1,s in masks];z=[(Float64((s>>i)&1)-Float64((s>>(L+i))&1))/2 for i in 0:L-1,s in masks];p=psi.^2;density=occupations*p;Czz=(z.*reshape(p,1,:))*z';Cnn=(occupations.*reshape(p,1,:))*occupations'-density*density'
    structure(C)=[real(sum(exp(2im*pi*q*(i-j)/L)*C[i+1,j+1] for i in 0:L-1,j in 0:L-1)/L) for q in 0:L-1]
    function pairmatrix(bond)
        states=Vector{Vector{Float64}}()
        for i in 0:L-1
            j=mod(i+1,L)
            terms=bond ? [(1/sqrt(2),[('-',i),('-',L+j)]),(-1/sqrt(2),[('-',L+i),('-',j)])] : [(1.,[('-',i),('-',L+i)])]
            push!(states,action(psi,masks,terms))
        end
        [dot(x,y) for x in states,y in states]
    end
    result=Dict{String,Any}("metadata"=>Dict("created_utc"=>string(now(UTC)),"job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"),"source_root"=>ROOT,"user_solver_source_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"analysis/compute_s2_ed.jl")))),"adapter_sha256"=>bytes2hex(sha256(read(@__FILE__))),"manifest_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"manifest.toml")))),"ground_energy_uncentered"=>E,"ground_residual"=>res,"ground_iterations"=>it,"julia_version"=>string(VERSION),"blas_threads"=>BLAS.get_num_threads()))
    for c in M["cases"]
        id=c["id"];c["ready"]||continue
        try
            values=if id=="P01"
                probs=max.(eigvals(Symmetric(rdm_adapter(psi,space,[0,1]))),0.);traces=[sum(probs.^q) for q in 2:4];vcat(traces,[log(traces[q-1])/(1-q) for q in 2:4])
            elseif id=="P02";spin_ladder()
            elseif id=="P03";vcat(rowflat(Cu),projections(Cu),rowflat(Cd),projections(Cd))
            elseif id=="P04";sq=structure(Czz);vcat(rowflat(Czz),sq,[Float64(argmax(sq)-1)])
            elseif id=="P05";vcat(rowflat(Cnn),structure(Cnn))
            elseif id=="P06";vcat(rowflat(pairmatrix(false)),rowflat(pairmatrix(true)))
            elseif id=="P16"
                ent=[ED.renyi2_entropy_sites(psi,space,s) for s in [[0],[2],[0,2]]];vcat(exp.(-ent),[ent[1]+ent[2]-ent[3]])
            else;error("Unsupported adapter")
            end
            result[id]=Dict("status"=>"TARGET_READY","values"=>values)
        catch ex
            result[id]=Dict("status"=>"UNRESOLVED","reason"=>sprint(showerror,ex,catch_backtrace()))
        end
        open(output,"w") do io;TOML.print(io,result);end
    end
    println("Wrote ",output)
end
main()

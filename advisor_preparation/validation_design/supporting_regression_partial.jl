#!/usr/bin/env julia
# Unfinished supporting regression draft: 25 registered kernel checks; no execution driver.
# These internal checks are excluded from the 100 external scientific comparisons.
using Test, LinearAlgebra, Random, TOML, Dates, SHA
const ROOT=get(ENV,"RENYI_ROOT",normpath(joinpath(@__DIR__,"../source_snapshot/HubbardRenyiDQMC")))
using HubbardRenyiDQMC
const H=HubbardRenyiDQMC
module ED
include(joinpath(Main.ROOT,"analysis/compute_s2_ed.jl"))
end
module DM
include(joinpath(Main.ROOT,"analysis/compute_s2_itensor_dmrg.jl"))
end
module SUM
include(joinpath(Main.ROOT,"analysis/summarize_grid.jl"))
end
module NODE
include(joinpath(Main.ROOT,"run/run_lambda_node.jl"))
end
module LAUNCH
include(joinpath(Main.ROOT,"run/run_grid_local.jl"))
end
module THERM
include(joinpath(Main.ROOT,"analysis/audit_u4_l4_thermal_ed.jl"))
end
using ITensors, ITensorMPS
BLAS.set_num_threads(1)
const CASES=Any[]
const CACHE=Dict{Any,Any}()
function case(f,id,title,fixture,oracle,acceptance;tier="T0")
    push!(CASES,(;f,id,title,fixture,oracle,acceptance,tier))
end
cached(f,k)=get!(f,CACHE,k)
# Independent Jordan-Wigner full Fock oracle: bits 0:L-1 are up, L:2L-1 down.
function annihilation(L,p)
    d=1<<(2L); a=zeros(d,d)
    for s in 0:d-1
        if !iszero(s&(1<<p))
            sign=isodd(count_ones(s&((1<<p)-1))) ? -1. : 1.
            a[(s ⊻ (1<<p))+1,s+1]=sign
        end
    end
    a
end
function oracle(L=4,U=2.;t=1.,mu=0.,boundary="pbc",centered=false)
    cached((:oracle,L,U,t,mu,boundary,centered)) do
        cs=[annihilation(L,p) for p in 0:2L-1]; ns=[c'c for c in cs]
        d=1<<(2L); K=zeros(d,d); V=zeros(d,d)
        bonds=[(j,j+1,-t) for j in 1:L-1]
        boundary!="obc" && push!(bonds,(L,1,boundary=="pbc" ? -t : t))
        for spin in 0:1,(i,j,v) in bonds
            a,b=i+spin*L,j+spin*L
            K .+= v.*(cs[a]'cs[b]+cs[b]'cs[a])
        end
        for i in 1:L
            V .+= U.*(ns[i]*ns[i+L])
            centered && (V .+= (-U/2).*(ns[i]+ns[i+L]) .+ (U/4).*Matrix{Float64}(I,d,d))
        end
        N=sum(ns); Hamiltonian=K+V-mu*N
        (;L,U,t,mu,boundary,centered,cs,ns,K,V,N,H=Hamiltonian)
    end
end
function sector(L=4,nu=2,nd=2)
    up=ED.fixed_particle_masks(L,nu); dn=ED.fixed_particle_masks(L,nd)
    idx=[Int(u)|(Int(d)<<L) for u in up for d in dn].+1
    (;up,dn,idx)
end
function edfixture(U=2.;boundary="apbc",L=4,nu=2,nd=2)
    cached((:ed,U,boundary,L,nu,nd)) do
        o=oracle(L,U;boundary); sec=sector(L,nu,nd)
        space=ED.build_space(L,1.,U,0.,nu,nd;boundary)
        Hb=o.H[sec.idx,sec.idx]; ev=eigen(Symmetric(Hb)); E,psi,it,res=ED.ground_state_lanczos(space;tol=1e-11,maxiter=200)
        (;o,sec,space,Hb,ev,E,psi,it,res)
    end
end
function full_state(f)
    v=zeros(1<<(2f.space.L)); v[f.sec.idx]=f.psi; v
end
# Signed regrouping of occupied creation operators to interleaved A then interleaved B.
function signed_labels(s,L,A)
    B=setdiff(collect(1:L),A)
    order=vcat([[j-1,j-1+L] for j in A]...,[[j-1,j-1+L] for j in B]...)
    occupied_order=[p for p in order if !iszero(s&(1<<p))]
    inversions=sum(occupied_order[i]>occupied_order[j] for i in eachindex(occupied_order) for j in i+1:length(occupied_order);init=0)
    label(sites)=1+sum((((s>>(j-1))&1)+2*((s>>(j-1+L))&1))*4^(k-1) for (k,j) in enumerate(sites);init=0)
    label(A),label(B),isodd(inversions) ? -1. : 1.
end
function coefficient_matrix(v,L,A)
    M=zeros(eltype(v),4^length(A),4^(L-length(A)))
    for s in 0:length(v)-1
        a,b,sgn=signed_labels(s,L,A); M[a,b]=sgn*v[s+1]
    end
    M
end
function pure_entropy(v,L,A)
    M=coefficient_matrix(v,L,A); rho=M*M'; rho ./= tr(rho)
    -log(real(tr(rho*rho))),rho
end
function thermal_density(o,beta)
    ev=eigen(Symmetric(o.H)); w=exp.(-beta.*(ev.values.-minimum(ev.values))); w./=sum(w)
    ev.vectors*Diagonal(w)*ev.vectors'
end
function reduced_density(rho,L,A)
    da,db=4^length(A),4^(L-length(A)); out=zeros(eltype(rho),da,da)
    labels=[signed_labels(s,L,A) for s in 0:size(rho,1)-1]
    for s in eachindex(labels),t in eachindex(labels)
        a,b,sgn=labels[s]; aa,bb,sgnt=labels[t]
        b==bb && (out[a,aa]+=sgn*sgnt*rho[s,t])
    end
    out
end
function to_full(psi,sites,L)
    T=prod(psi); v=zeros(1<<(2L))
    for s in 0:length(v)-1
        invs=sum(((s>>(i-1+L))&1)*((s>>(j-1))&1) for i in 1:L for j in i+1:L;init=0)
        labels=[1+((s>>(j-1))&1)+2*((s>>(j-1+L))&1) for j in 1:L]
        v[s+1]=(isodd(invs) ? -1. : 1.)*real(T[(sites[j]=>labels[j] for j in 1:L)...])
    end
    v
end
function from_full(v,sites,L)
    T=ITensor(Float64,sites...)
    for s in 0:length(v)-1
        iszero(v[s+1]) && continue
        invs=sum(((s>>(i-1+L))&1)*((s>>(j-1))&1) for i in 1:L for j in i+1:L;init=0)
        labels=[1+((s>>(j-1))&1)+2*((s>>(j-1+L))&1) for j in 1:L]
        T[(sites[j]=>labels[j] for j in 1:L)...]=(isodd(invs) ? -1. : 1.)*v[s+1]
    end
    MPS(T,sites;cutoff=0.,maxdim=256)
end
function dmfixture(U=2.;boundary="apbc",seed=1123,maxdim=64,nsweeps=12,cutoff=1e-12)
    cached((:dm,U,boundary,seed,maxdim,nsweeps,cutoff)) do
        sites=siteinds("Electron",4;conserve_qns=false)
        mpo=DM.hubbard_mpo(sites;L=4,t=1.,U,mu=0.,boundary)
        Random.seed!(seed); init=random_mps(sites,["Up","Dn","Up","Dn"];linkdims=8)
        E,psi=dmrg(mpo,init;nsweeps,maxdim,cutoff,outputlevel=0)
        v=to_full(psi,sites,4); f=edfixture(U;boundary)
        (;sites,mpo,E,psi,v,f)
    end
end
function direct_green(rep,spin=:up)
    Bs=spin==:up ? rep.Bup : rep.Bdn; P=Matrix{Float64}(I,rep.L,rep.L)
    for B in Bs
        P=Diagonal(B.expmΔτV)*B.expmΔτK*P
    end
    inv(I+P)
end
replica(;L=4,U=2.,beta=.2,dtau=.1,seed=101,n_stab=1,bc=:pbc)=H.initializer(1.,U,0.,beta,L,dtau,seed;n_stab,boundary_condition=bc)
function tracker_fixture()
    rng=Xoshiro(712); N=4; C=[1,3]; G1=.4Matrix{Float64}(I,N,N)+.01randn(rng,N,N); G2=.6Matrix{Float64}(I,N,N)+.01randn(rng,N,N)
    tr=H.GroverRank1Tracker(Float64,2); H.set_passive_grover_block!(tr,G2,C); H.rebuild_active_grover!(tr,G1,C)
    gt=.1randn(rng,N,N); g0=.1randn(rng,N,N)
    (;tr,G1,G2,C,gt,g0)
end
function tracker_trial(f)
    p=(.1/1.2).*f.g0[f.C,2]; q=vec(f.gt[2,f.C]); f.G1[f.C,f.C]+p*q'
end
const META=(Ninc=1,S2_ED=nothing,root_seeds=[7],N_bins=2,run_name="fixture")
const ROW0=["fixture,0,0.0,7,1,1.0,1.0,0.1,0.80,","fixture,0,0.0,7,2,1.0,1.0,0.1,0.82,"]
const ROW1=["fixture,1,1.0,7,1,1.0,2.0,1.0,,1.20","fixture,1,1.0,7,2,1.0,2.0,1.0,,1.22"]
function writegrid(dir;rows0=ROW0,rows1=ROW1,header=SUM.BINS_HEADER)
    for (k,rows) in [(0,rows0),(1,rows1)]
        rows===nothing && continue
        open(joinpath(dir,"fixture_k$(k)_seed7_bins.csv"),"w") do io
            println(io,header); foreach(row->println(io,row),rows)
        end
    end
end
# 25 DQMC and optimized-kernel tests.
case("DQMC-01","Boundary input validation","pbc/apbc/invalid tokens","declared boundary enum","invalid tokens throw") do
    @test H.normalize_boundary_condition("pbc")==:pbc
    @test H.normalize_boundary_condition("apbc")==:apbc
    @test_throws ArgumentError H.normalize_boundary_condition("periodic")
end
case("DQMC-02","Weighted Hubbard-Stratonovich sum","U2.3,dtau.1, four occupations","explicit onsite exponential","atol=2e-15,rtol=2e-15") do
    U=2.3;dt=.1;a=acosh(exp(dt*U/2))
    for nu in 0:1,nd in 0:1
        lhs=exp(-dt*U*(nu-.5)*(nd-.5)); rhs=.5exp(-dt*U/4)*sum(exp(a*s*(nu-nd)) for s in [-1,1])
        @test isapprox(lhs,rhs;atol=2e-15,rtol=2e-15)
    end
end
case("DQMC-03","PBC free dispersion","L8,t1,PBC","-2cos(2pi n/L)","eigenvalue atol=1e-12") do
    K=H.hubbard_chain_kinetic_matrix(1.,8,:pbc)
    @test sort(eigvals(Symmetric(K)))≈sort([-2cos(2pi*n/8) for n in 0:7]) atol=1e-12
end
case("DQMC-04","APBC free dispersion","L8,t1,APBC","-2cos(2pi(n+1/2)/L)","eigenvalue atol=1e-12") do
    K=H.hubbard_chain_kinetic_matrix(1.,8,:apbc)
    @test sort(eigvals(Symmetric(K)))≈sort([-2cos(2pi*(n+.5)/8) for n in 0:7]) atol=1e-12
end
case("DQMC-05","Boundary seam changes only two entries","L8,PBC/APBC","explicit seam hopping sign","exact matrix pattern") do
    p=H.hubbard_chain_kinetic_matrix(1.,8,:pbc);a=H.hubbard_chain_kinetic_matrix(1.,8,:apbc)
    @test a==a'; @test count(!iszero,a-p)==2; @test a[1,8]==a[8,1]==1.
end
case("DQMC-06","Bipartite particle-hole kinetic symmetry","L8,PBC/APBC","staggered sublattice unitary","exact matrix identity") do
    D=Diagonal([(-1.)^j for j in 1:8]);for bc in [:pbc,:apbc]
        K=H.hubbard_chain_kinetic_matrix(1.,8,bc); @test D*K*D == -K
    end
end
case("DQMC-07","Free equal-time Green function","L4,U0,beta.2,dtau.1","independent Fermi kernel inv(I+exp(-beta K))","atol=1e-11") do
    r=replica(U=0.);K=H.hubbard_chain_kinetic_matrix(1.,4)
    @test r.Gup≈inv(I+exp(-.2K)) atol=1e-11
    @test r.Gdn≈r.Gup atol=1e-12
end
case("DQMC-08","Fixed-field Green full product","L4,U2,beta.2,two slices","direct dense propagator multiplication","atol=1e-10") do
    r=replica();@test r.Gup≈direct_green(r) atol=1e-10; @test r.Gdn≈direct_green(r,:dn) atol=1e-10
end
case("DQMC-09","Fermion determinant bookkeeping","same fixed-field replica","logabsdet of direct Green","logdet atol=1e-10, exact sign") do
    r=replica();lg,sg=logabsdet(direct_green(r)); @test r.logdetGup≈lg atol=1e-10; @test r.sgndetGup==sg
end
case("DQMC-10","HS field dimensions and reproducibility","L4,beta.2,dtau.1,seed101","integer grid and repeated explicit seed","exact equality") do
    a=replica();b=replica();@test size(a.s)==(4,2);@test all(x->x in [-1,1],a.s);@test a.s==b.s
end
case("DQMC-11","Local G00 rank-one formula","random real 4x4 arrays","explicit outer product","atol=1e-13") do
    rng=Xoshiro(11);g=randn(rng,4,4);a=randn(rng,4,4);b=randn(rng,4,4);expected=g+.1/1.2*b[:,2]*a[2,:]';actual=copy(g)
    H.local_update_G00!(actual,a,b,1.2,.1,2,zeros(4),zeros(4));@test actual≈expected atol=1e-13
end
case("DQMC-12","Local Gtau0 rank-one formula","random real 4x4 arrays","explicit outer product from fixed Gtautau","atol=1e-13") do
    rng=Xoshiro(12);g=randn(rng,4,4);a=randn(rng,4,4);expected=g+.1/1.2*a[:,2]*g[2,:]';actual=copy(g)
    H.local_update_Gτ0!(actual,a,1.2,.1,2,zeros(4),zeros(4));@test actual≈expected atol=1e-13
end
case("DQMC-13","Local G0tau contact subtraction","random real 4x4 arrays","explicit outer product with Gtautau-I row","atol=1e-13") do
    rng=Xoshiro(13);g=randn(rng,4,4);a=randn(rng,4,4);row=copy(a[2,:]);row[2]-=1;expected=g+.1/1.2*g[:,2]*row';actual=copy(g)
    H.local_update_G0τ!(actual,a,1.2,.1,2,zeros(4),zeros(4));@test actual≈expected atol=1e-13
end
case("DQMC-14","Stabilization grouping with identical HS fields","L4,U2,beta.2,period1/2","same seed and dense product oracle","atol=1e-10") do
    a=replica(n_stab=1);b=replica(n_stab=2);@test a.s==b.s;@test a.Gup≈b.Gup atol=1e-10;@test b.Gup≈direct_green(b) atol=1e-10
end
case("DQMC-15","Free Green time-step independence","L4,U0,beta.2,dtau.1/.05","commuting free propagators","atol=1e-11") do
    a=replica(U=0.);b=replica(U=0.,dtau=.05);@test a.Gup≈b.Gup atol=1e-11
end
case("DQMC-16","Zero-lambda fused update equals ordinary update","same-seed active/passive replicas","unmodified base Markov kernel at lambda0","same accepted fields, G atol1e-9") do
    a=replica();b=deepcopy(a);p=replica(seed=102);a,da,aa=H.update_replica!(a,0.);b,db,ab=H.update_replica_fused!(b,p,0.,[1,2],0.)
    @test a.s==b.s;@test aa==ab;@test a.Gup≈b.Gup atol=1e-9
end
case("DQMC-17","Optimized and direct fused trajectories","L4,U2,beta.2,lambda.5,4 rounds","full Grover determinant implementation with same RNG","same fields/signs/acceptance, G atol1e-9") do
    a=replica();p=replica(seed=102);b=deepcopy(a);q=deepcopy(p)
    for _ in 1:4
        a,_,aa=H.update_replica_fused!(a,p,0.,[1,2],.5);b,_,ab=H.update_replica_rank1_fused!(b,q,0.,[1,2],.5)
        @test a.s==b.s;@test aa==ab;@test a.Gup≈b.Gup atol=1e-9
        p,_,_=H.update_replica_fused!(p,a,0.,[1,2],.5);q,_,_=H.update_replica_rank1_fused!(q,b,0.,[1,2],.5)
        @test p.s==q.s
    end
end
case("DQMC-18","Grover determinant lemma","nC2,seed18,20 random draws","ratio of independent full determinants","atol1e-10,rtol1e-9") do
    rng=Xoshiro(18)
    for _ in 1:20
        a=randn(rng,2,2);b=randn(rng,2,2);p=randn(rng,2);q=randn(rng,2);g=a*b+(I-a)*(I-b)
        abs(det(g))<1e-6 && continue
        gn=(a+p*q')*b+(I-a-p*q')*(I-b);r=(2b-I)'q
        @test det(gn)/det(g)≈1+dot(r,g\p) atol=1e-10 rtol=1e-9
    end
end
case("DQMC-19","Tracker proposed log determinant","seed712,nC2,site2","explicit trial matrix logabsdet","atol1e-10") do
    f=tracker_fixture();v=H.prepare_grover_rank1_proposal!(f.tr,f.g0,f.gt,1.2,.1,2,f.C);a=tracker_trial(f);g=a*f.G2[f.C,f.C]+(I-a)*(I-f.G2[f.C,f.C]);@test v≈logabsdet(g)[1] atol=1e-10
end
case("DQMC-20","Accepted tracker matrix and inverse","same proposal fixture","full rebuilt Grover matrix and inverse","atol1e-9") do
    f=tracker_fixture();H.prepare_grover_rank1_proposal!(f.tr,f.g0,f.gt,1.2,.1,2,f.C);H.accept_grover_rank1_proposal!(f.tr);a=tracker_trial(f);g=a*f.G2[f.C,f.C]+(I-a)*(I-f.G2[f.C,f.C]);@test f.tr.grover≈g atol=1e-9;@test f.tr.inverse_grover*g≈Matrix{Float64}(I,2,2) atol=1e-9
end
case("DQMC-21","Unaccepted proposal preserves accepted state","same proposal fixture","bytewise copied accepted arrays","exact equality") do
    f=tracker_fixture();g=copy(f.tr.grover);a=copy(f.tr.G1C);iv=copy(f.tr.inverse_grover);H.prepare_grover_rank1_proposal!(f.tr,f.g0,f.gt,1.2,.1,2,f.C);@test f.tr.grover==g;@test f.tr.G1C==a;@test f.tr.inverse_grover==iv
end
case("DQMC-22","Explicit tracker rebuild after accepted update","same fixture, accepted then rebuild","dense product using accepted G1C","atol1e-10") do
    f=tracker_fixture();H.prepare_grover_rank1_proposal!(f.tr,f.g0,f.gt,1.2,.1,2,f.C);H.accept_grover_rank1_proposal!(f.tr);f.G1[f.C,f.C]=tracker_trial(f);H.rebuild_active_grover!(f.tr,f.G1,f.C);@test f.tr.accepted_since_rebuild==0;@test f.tr.G1C≈f.G1[f.C,f.C] atol=1e-10
end
case("DQMC-23","Ill-conditioned inverse shortcut disabled","nC2,Grover eigenvalue near0","explicit logabsdet and inverse conditioning","unusable inverse or explicit rebuild") do
    f=tracker_fixture();a=.5Matrix{Float64}(I,2,2);b=copy(a);a[1,1]=1.;b[1,1]=1e-16
    H.set_passive_grover_block!(f.tr,b,[1,2]);H.rebuild_active_grover!(f.tr,a,[1,2]);@test !f.tr.usable_inverse
end
case("DQMC-24","Fused log-ratio special cases","lambda0, finite logs and zero determinants","power weight exp(lambda delta logX)","exact special values, atol1e-14") do
    @test H.fused_ratio_from_logs(-Inf,-Inf,0.)==1.;@test H.fused_ratio_from_logs(-Inf,0.,.5)==0.;@test H.fused_ratio_from_logs(0.,-Inf,.5)==Inf;@test H.fused_ratio_from_logs(log(4.),0.,.5)≈2. atol=1e-14
end
case("DQMC-25","Short sweep Green matrices versus full products","L4,U2,beta.2,5 complete update rounds","recomputed dense propagator product","atol1e-8") do
    r=replica()
    for _ in 1:5
        r,_,_=H.update_replica!(r,0.)
        @test r.Gup_00≈direct_green(r) atol=1e-8
        @test r.Gdn_00≈direct_green(r,:dn) atol=1e-8
    end
end

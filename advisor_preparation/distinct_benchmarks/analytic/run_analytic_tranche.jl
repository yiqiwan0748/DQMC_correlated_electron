#!/usr/bin/env julia
# Additive scientific adapters. No production or archived source is modified.
using LinearAlgebra, SparseArrays, Statistics, TOML, Dates, SHA, Sockets, Printf
BLAS.set_num_threads(1)
const ROOT=normpath(get(ENV,"RENYI_ROOT",joinpath(@__DIR__,"..","..","validated_package","HubbardRenyiDQMC")))
module UserED
include(joinpath(Main.ROOT,"analysis","compute_s2_ed.jl"))
end
module UserThermal
include(joinpath(Main.ROOT,"analysis","audit_u4_l4_thermal_ed.jl"))
end
const TE=UserThermal.U4L4ThermalEDAudit
const THERMAL_CACHE=Dict{Any,Any}()
const OPS_CACHE=Dict{Int,Any}()
const IMPLEMENTED=["A"*lpad(i,2,'0') for i in vcat(collect(1:14),[22,23])]

flat(x::Number)=[Float64(x)]
flat(x)=Float64.(vec(x))
function check!(fixture,name,observed,reference;tol=1e-10,metric="maximum_absolute")
    a,b=flat(observed),flat(reference)
    length(a)==length(b)||error("Dimension mismatch: $name")
    err=all(isfinite,a)&&all(isfinite,b) ? (metric=="L1" ? sum(abs,a-b) : maximum(abs,a-b;init=0.)) : Inf
    push!(fixture["checks"],Dict("observable"=>name,"observed"=>a,"analytical_reference"=>b,
        "metric"=>metric,"error"=>err,"tolerance"=>tol,"status"=>err<=tol ? "PASS" : "FAIL"))
end
function fixture!(result,parameters)
    f=Dict{String,Any}("parameters"=>parameters,"checks"=>Any[])
    push!(result["fixtures"],f);f
end
function dense_hamiltonian(space)
    n=UserED.dimension(space);H=zeros(n,n);x=zeros(n);y=zeros(n)
    for j in 1:n
        fill!(x,0.);x[j]=1.;UserED.apply_hamiltonian!(y,x,space);H[:,j]=y
    end
    maximum(abs,H-H')<=1e-12||error("User Hamiltonian action is not Hermitian")
    H
end
function ground(L,t,U,Nup,Ndn;boundary="pbc")
    space=UserED.build_space(L,t,U,0.,Nup,Ndn;boundary)
    E,psi,it,res=UserED.ground_state_lanczos(space;maxiter=max(300,UserED.dimension(space)),tol=1e-12)
    H=dense_hamiltonian(space);d=eigen(Symmetric(H))
    abs(E-first(d.values))<=1e-10||error("Lanczos state is not the dense certified minimum")
    norm(H*psi-E*psi)<=1e-10||error("Eigenvector residual exceeds certification threshold")
    (;space,E,psi,H,d,it,res)
end
index_full(u,d,L)=Int(u)*(1<<L)+Int(d)+1
function embed(psi,space)
    v=zeros(1<<(2space.L));nd=length(space.dn_masks)
    for (iu,u) in enumerate(space.up_masks),(id,d) in enumerate(space.dn_masks)
        v[index_full(u,d,space.L)]=psi[(iu-1)*nd+id]
    end
    v
end
function particle_vectors(L)
    masks=collect(UInt(0):UInt((1<<L)-1))
    N=Float64[count_ones(u)+count_ones(d) for u in masks for d in masks]
    D=Float64[count_ones(u&d) for u in masks for d in masks]
    N,D
end
function operators(L)
    get!(OPS_CACHE,L) do
        dim=1<<(2L);maskL=UInt((1<<L)-1);c=Any[]
        for mode in 0:2L-1
            mat=spzeros(Float64,dim,dim)
            for u in UInt(0):maskL,d in UInt(0):maskL
                bits=u|(d<<L);bit=UInt(1)<<mode
                iszero(bits&bit)&&continue
                dest=xor(bits,bit);sgn=isodd(count_ones(bits&(bit-1))) ? -1. : 1.
                mat[index_full(dest&maskL,dest>>L,L),index_full(u,d,L)]=sgn
            end
            push!(c,mat)
        end
        N,D=particle_vectors(L);sz=spzeros(Float64,dim,dim);plus=spzeros(Float64,dim,dim)
        for i in 1:L
            sz+=(c[i]'*c[i]-c[L+i]'*c[L+i])/2
            plus+=c[i]'*c[L+i]
        end
        S2=sz*sz+(plus*plus'+plus'*plus)/2
        (;c,N,D,S2,pair=[c[L+i]*c[i] for i in 1:L])
    end
end
function thermal(L,t,U,mu,beta,boundary="pbc")
    get!(THERMAL_CACHE,(L,t,U,mu,beta,boundary)) do
        K,V,masks=TE.full_fock_matrices(L,t,U,mu;boundary)
        H=K+Diagonal(V);rho,d=TE.exact_thermal_density(H,beta)
        p=exp.(-beta.*(d.values.-first(d.values)));logZ=-beta*first(d.values)+log(sum(p));p./=sum(p)
        (;L,t,U,mu,beta,boundary,K,V,masks,H,rho,d,p,logZ)
    end
end
function source_propagator(a,tau)
    cs=[a.d.vectors'*c*a.d.vectors for c in operators(a.L).c[1:a.L]]
    E=a.d.values.-first(a.d.values)
    W=exp.(-(a.beta-tau).*E).*transpose(exp.(-tau.*E))
    W./=sum(exp.(-a.beta.*E))
    [sum(W.*cs[i].*conj.(cs[j])) for i in 1:a.L,j in 1:a.L]
end
function free_modes(L,t,mu,beta,boundary)
    phi=boundary=="apbc" ? pi : 0.
    k=[(2pi*m+phi)/L for m in 0:L-1];eps=-2t.*cos.(k);xi=eps.-mu
    f=1 ./ (1 .+ exp.(beta.*xi));C=[sum(f.*exp.(im.*k.*(j-i)))/L for i in 0:L-1,j in 0:L-1]
    (;k,eps,xi,f,C)
end
function atomic(U,mu,beta)
    w=[exp(-beta*U/4),exp(beta*(U/4+mu)),exp(beta*(U/4+mu)),exp(-beta*(U/4-2mu))];p=w/sum(w)
    N=p[2]+p[3]+2p[4];D=p[4]
    (;p,N,D,E=U*(D-N/2+1/4),logZ=log(sum(w)))
end
function source_onebody(g)
    s=g.space;nd=length(s.dn_masks);index=UserED.mask_index(s.up_masks);C=zeros(s.L,s.L)
    for i in 0:s.L-1,j in 0:s.L-1,(iu,u) in enumerate(s.up_masks),(id,d) in enumerate(s.dn_masks)
        src=(iu-1)*nd+id
        if i==j
            C[i+1,j+1]+=abs2(g.psi[src])*(UserED.occupied(u,i) ? 1. : 0.)
        elseif UserED.occupied(u,j)&&!UserED.occupied(u,i)
            dest=xor(xor(u,UInt(1)<<j),UInt(1)<<i)
            C[i+1,j+1]+=g.psi[(index[dest]-1)*nd+id]*g.psi[src]*UserED.fermion_hop_sign(u,j,i)
        end
    end
    C
end

function atomic_eos!(r)
    for U in [0.,2.,6.],mu in [-.7,0.,.6],beta in [.25,2.]
        f=fixture!(r,Dict("L"=>1,"t"=>0.,"U"=>U,"mu_centered"=>mu,"beta"=>beta))
        a=thermal(1,0.,U,mu,beta);ref=atomic(U,mu,beta);N,D=particle_vectors(1)
        check!(f,"N",dot(diag(a.rho),N),ref.N);check!(f,"D",dot(diag(a.rho),D),ref.D)
        check!(f,"centered_internal_energy",tr(a.rho*a.H)+mu*dot(diag(a.rho),N),ref.E)
        check!(f,"logZ",a.logZ,ref.logZ)
    end
end
function atomic_poles!(r)
    for U in [1.,4.],mu in [-.3,.7],beta in [.5,2.]
        f=fixture!(r,Dict("L"=>1,"U"=>U,"mu_centered"=>mu,"beta"=>beta))
        a=thermal(1,0.,U,mu,beta);ref=atomic(U,mu,beta);nbar=ref.N/2
        poles=[-U/2-mu,U/2-mu];residues=[1-nbar,nbar]
        C=a.d.vectors'*operators(1).c[1]*a.d.vectors;observed=zeros(2)
        for i in eachindex(a.p),j in eachindex(a.p)
            weight=(a.p[i]+a.p[j])*abs2(C[i,j]);weight<1e-14&&continue
            pole=a.d.values[j]-a.d.values[i];b=argmin(abs.(poles.-pole))
            check!(f,"transition_pole_$(i)_$(j)",pole,poles[b]);observed[b]+=weight
        end
        check!(f,"two_pole_residues",observed,residues)
        for q in [0.,1/8,1/3,1/2,7/8,1.]
            tau=q*beta;G=sum(residues.*exp.(-tau.*poles)./(1 .+ exp.(-beta.*poles)))
            check!(f,"Gplus_tau_fraction_$q",source_propagator(a,tau)[1,1],G)
        end
    end
end
function dimer!(r,kind)
    for U in [0.,1.,4.,12.],t in [.7,1.3]
        f=fixture!(r,Dict("physical_one_bond_t"=>t,"source_L2_ring_t"=>t/2,"U"=>U,"mu_uncentered"=>0.))
        g=ground(2,t/2,U,1,1);R=sqrt(U^2+16t^2);E2=(U-R)/2
        check!(f,"singlet_ground_energy",g.E,E2;tol=1e-9)
        if kind==3
            trip=ground(2,t/2,U,2,0);ops=operators(2);v=embed(g.psi,g.space);vt=embed(trip.psi,trip.space)
            check!(f,"singlet_S2",dot(v,ops.S2*v),0.);check!(f,"triplet_S2",dot(vt,ops.S2*vt),2.)
            check!(f,"spin_gap",trip.E-g.E,(R-U)/2;tol=1e-9)
        elseif kind==4
            one=ground(2,t/2,U,1,0);three=ground(2,t/2,U,2,1)
            check!(f,"N1_energy",one.E,-abs(t);tol=1e-9);check!(f,"N3_energy",three.E,U-abs(t);tol=1e-9)
            check!(f,"charge_gap",three.E+one.E-2g.E,R-2abs(t);tol=1e-9)
        else
            v=embed(g.psi,g.space);rhoA=TE.reduced_operator(v*v',collect(UInt(0):UInt(3)),2,1)
            D=(1-U/R)/2;p=sort([D/2,D/2,(1-D)/2,(1-D)/2])
            check!(f,"one_site_RDM_eigenvalues",eigvals(Symmetric(rhoA)),p)
            check!(f,"user_S2_functional",UserED.renyi2_entropy(g.psi,g.space,1),-log(sum(abs2,p)))
        end
    end
end
function free_thermo!(r)
    for bc in ["pbc","apbc"],beta in [.25,1.,3.],mu in [-.4,0.,.6]
        f=fixture!(r,Dict("L"=>4,"t"=>1.2,"U"=>0.,"beta"=>beta,"mu_centered"=>mu,"boundary"=>bc))
        a=thermal(4,1.2,0.,mu,beta,bc);b=free_modes(4,1.2,mu,beta,bc);N,D=particle_vectors(4)
        check!(f,"density_per_site",dot(diag(a.rho),N)/4,2sum(b.f)/4)
        check!(f,"kinetic_energy_total",tr(a.rho*a.K),2dot(b.eps,b.f))
        check!(f,"grand_log_partition",a.logZ,2sum(log1p.(exp.(-beta.*b.xi))))
        check!(f,"thermal_entropy_total",-sum(p>0 ? p*log(p) : 0. for p in a.p),-2sum(b.f.*log.(b.f).+(1 .-b.f).*log.(1 .-b.f)))
    end
end
function free_propagator!(r)
    for bc in ["pbc","apbc"],beta in [.5,2.],mu in [-.3,.5]
        f=fixture!(r,Dict("L"=>4,"U"=>0.,"beta"=>beta,"mu_centered"=>mu,"boundary"=>bc,"tau_endpoint"=>"0plus and beta-minus"))
        a=thermal(4,1.,0.,mu,beta,bc);b=free_modes(4,1.,mu,beta,bc)
        for q in [0.,1/8,1/3,1/2,7/8,1.]
            tau=q*beta;G=source_propagator(a,tau)
            target=[sum(exp.(-tau.*b.xi).*exp.(im.*b.k.*(i-j))./(1 .+ exp.(-beta.*b.xi)))/4 for i in 0:3,j in 0:3]
            check!(f,"Gplus_real_tau_$q",real.(G),real.(target);tol=1e-9)
            check!(f,"Gplus_imag_tau_$q",imag.(G),imag.(target);tol=1e-9)
        end
    end
end
function lindhard!(r)
    L=4;beta=1.3;mu=.2;bc="apbc";a=thermal(L,1.,0.,mu,beta,bc);b=free_modes(L,1.,mu,beta,bc)
    for q in [0.,pi/2,pi],m in [0,1,2]
        f=fixture!(r,Dict("L"=>L,"beta"=>beta,"mu_centered"=>mu,"boundary"=>bc,"q"=>q,"Matsubara_m"=>m));omega=2pi*m/beta
        nq=ComplexF64[sum(exp(-im*q*i)*(((u>>i)&1)+((d>>i)&1)) for i in 0:L-1) for u in 0:(1<<L)-1 for d in 0:(1<<L)-1]
        O=a.d.vectors'*Diagonal(nq)*a.d.vectors;chi=0im
        for i in eachindex(a.p),j in eachindex(a.p)
            gap=a.d.values[j]-a.d.values[i]
            factor=abs(gap)<1e-10 ? (m==0 ? beta*a.p[i] : 0.) : (a.p[i]-a.p[j])/(gap-im*omega)
            chi+=factor*abs2(O[i,j])/L
        end
        m==0&&(chi-=beta*abs2(dot(diag(a.rho),nq))/L)
        target=0im
        for k in eachindex(b.k)
            ep=-2cos(b.k[k]+q);fp=1/(1+exp(beta*(ep-mu)));gap=ep-b.eps[k]
            target+=2/L*(abs(gap)<1e-10 ? (m==0 ? beta*b.f[k]*(1-b.f[k]) : 0.) : (b.f[k]-fp)/(gap-im*omega))
        end
        check!(f,"connected_chi_real",real(chi),real(target);tol=1e-9);check!(f,"connected_chi_imag",imag(chi),imag(target);tol=1e-9)
    end
end
function spectral_poles!(r)
    a=thermal(4,1.,0.,.2,1.3,"apbc");b=free_modes(4,1.,.2,1.3,"apbc")
    for k in b.k
        f=fixture!(r,Dict("L"=>4,"beta"=>1.3,"mu_centered"=>.2,"boundary"=>"apbc","k"=>k,"pole_support_weight_threshold"=>1e-13))
        ck=sum(exp(-im*k*i)*operators(4).c[i+1]/2 for i in 0:3);C=a.d.vectors'*ck*a.d.vectors;pole=-2cos(k)-.2
        total=0.;maxerr=0.;discarded=0.
        for i in eachindex(a.p),j in eachindex(a.p)
            w=(a.p[i]+a.p[j])*abs2(C[i,j]);total+=w
            if w<1e-13;discarded+=w;continue;end
            maxerr=max(maxerr,abs(a.d.values[j]-a.d.values[i]-pole))
        end
        check!(f,"spectral_pole_maximum_error",maxerr,0.;tol=1e-9);check!(f,"spectral_weight",total,1.;tol=1e-9)
        check!(f,"discarded_spectral_support_weight",discarded,0.;tol=1e-10)
    end
end
function compressibility!(r)
    for bc in ["pbc","apbc"],beta in [.5,2.],mu in [-.3,.5]
        f=fixture!(r,Dict("L"=>4,"U"=>0.,"beta"=>beta,"mu_centered"=>mu,"boundary"=>bc));a=thermal(4,1.,0.,mu,beta,bc);b=free_modes(4,1.,mu,beta,bc);N,D=particle_vectors(4)
        meanN=dot(diag(a.rho),N);varN=dot(diag(a.rho),N.^2)-meanN^2
        check!(f,"dn_dmu_via_thermal_number_fluctuation",beta*varN/4,2beta*sum(b.f.*(1 .-b.f))/4;tol=1e-9)
    end
end
function gaussian_levels!(r)
    for bc in ["pbc","apbc"],beta in [.5,3.],mu in [-.3,.5],nA in [1,2,3]
        f=fixture!(r,Dict("L"=>4,"U"=>0.,"beta"=>beta,"mu_centered"=>mu,"boundary"=>bc,"nA"=>nA));a=thermal(4,1.,0.,mu,beta,bc);b=free_modes(4,1.,mu,beta,bc)
        nu=eigvals(Hermitian(b.C[1:nA,1:nA]));p=[1.]
        for v in vcat(nu,nu);p=vcat(p.*v,p.*(1-v));end
        rhoA=TE.reduced_operator(a.rho,a.masks,4,nA)
        check!(f,"complete_RDM_spectrum_L1",eigvals(Symmetric(rhoA)),sort(p);tol=1e-9,metric="L1")
    end
end
choose(n,k)=0<=k<=n ? binomial(n,k) : 0
function canonical_infinite!(r)
    for (nu,nd) in [(1,1),(2,2),(3,1)],U in [0.,4.],nA in [1,2,3]
        f=fixture!(r,Dict("L"=>4,"U"=>U,"Nup"=>nu,"Ndown"=>nd,"beta"=>0.,"nA"=>nA))
        space=UserED.build_space(4,1.,U,.4,nu,nd);H=dense_hamiltonian(space);rho,d=TE.exact_thermal_density(H,0.)
        inds=[index_full(u,d,4) for u in space.up_masks for d in space.dn_masks];full=zeros(256,256);full[inds,inds]=rho
        rhoA=TE.reduced_operator(full,collect(UInt(0):UInt(15)),4,nA)
        purity=prod(sum(choose(nA,j)*(choose(4-nA,N-j)/choose(4,N))^2 for j in 0:nA) for N in [nu,nd])
        check!(f,"canonical_RDM_purity",tr(rhoA*rhoA),purity)
        check!(f,"canonical_S2",-log(real(tr(rhoA*rhoA))),-log(purity))
    end
end
function grand_infinite!(r)
    for U in [0.,4.],nA in [1,2,3,4]
        f=fixture!(r,Dict("L"=>4,"U"=>U,"mu_centered"=>.7,"beta"=>0.,"nA"=>nA));a=thermal(4,1.,U,.7,0.)
        S2,purity,rhoA=TE.renyi2(a.rho,a.masks,4,nA)
        check!(f,"GCE_S2",S2,nA*log(4.));check!(f,"GCE_purity",purity,4.0^(-nA))
    end
end
function polarized!(r)
    for (L,N,bc) in [(4,1,"pbc"),(4,2,"apbc"),(6,3,"pbc"),(6,2,"apbc")],U in [0.,3.,10.]
        f=fixture!(r,Dict("L"=>L,"Nup"=>N,"Ndown"=>0,"U"=>U,"boundary"=>bc));g=ground(L,1.,U,N,0;boundary=bc);b=free_modes(L,1.,0.,1.,bc);order=sortperm(b.eps);occ=order[1:N]
        b.eps[order[N+1]]-b.eps[order[N]]>1e-8||error("Polarized fixture is not a closed shell")
        C=[sum(exp.(im.*b.k[occ].*(j-i)))/L for i in 0:L-1,j in 0:L-1]
        check!(f,"spinless_free_ground_energy",g.E,sum(b.eps[occ]);tol=1e-9)
        check!(f,"onebody_real",source_onebody(g),real.(C);tol=1e-9);check!(f,"onebody_imag",zeros(L,L),imag.(C);tol=1e-9)
        D=sum(abs2(g.psi[(iu-1)*length(g.space.dn_masks)+id])*count_ones(u&d) for (iu,u) in enumerate(g.space.up_masks),(id,d) in enumerate(g.space.dn_masks))
        check!(f,"double_occupancy_total",D,0.;tol=1e-12)
    end
end
function eta_pairs!(r)
    L=4;ops=operators(L)
    for U in [1.,4.],m in [1,2,3]
        f=fixture!(r,Dict("L"=>L,"U"=>U,"eta_pairs"=>m,"boundary"=>"pbc","ensemble"=>"specified excited eta tower, not ground state"));space=UserED.build_space(L,1.,U,0.,m,m);psi=zeros(UserED.dimension(space));nd=length(space.dn_masks)
        for (iu,u) in enumerate(space.up_masks),(id,d) in enumerate(space.dn_masks)
            u==d||continue
            phase=sum(i for i in 0:L-1 if UserED.occupied(u,i));psi[(iu-1)*nd+id]=(isodd(phase) ? -1. : 1.)/sqrt(binomial(L,m))
        end
        H=dense_hamiltonian(space);E=m*U;check!(f,"eta_eigenvector_residual",norm(H*psi-E*psi),0.;tol=1e-9)
        d=eigen(Symmetric(H));support=findall(e->abs(e-E)<1e-9,d.values)
        check!(f,"exact_excited_eigenspace_overlap",sum(abs2,d.vectors[:,support]'*psi),1.;tol=1e-9)
        v=embed(psi,space);pv=[p*v for p in ops.pair];P=[dot(pv[i],pv[j]) for i in 1:L,j in 1:L]
        target=[i==j ? m/L : (-1.0)^(i+j)*m*(L-m)/(L*(L-1)) for i in 1:L,j in 1:L]
        check!(f,"staggered_pair_density_matrix",P,target;tol=1e-9)
    end
end
function lieb_spin!(r)
    L=3;ops=operators(L);N,D=particle_vectors(L);inds=findall(==(Float64(L)),N)
    for U in [.5,2.,8.]
        f=fixture!(r,Dict("L"=>L,"U"=>U,"N"=>L,"boundary"=>"open via average of actual PBC/APBC seam matrices","sublattice_sizes"=>[2,1]))
        kp,vp,masks=TE.full_fock_matrices(L,1.,U,0.;boundary="pbc");ka,va,_=TE.full_fock_matrices(L,1.,U,0.;boundary="apbc")
        H=(kp+ka)/2+Diagonal(vp.+U.*N./2 .-U*L/4);Hs=H[inds,inds];S=Matrix(ops.S2[inds,inds]);d=eigen(Symmetric(Hs));ground=findall(e->abs(e-first(d.values))<1e-9,d.values)
        check!(f,"ground_multiplet_dimension",length(ground),2.;tol=0.)
        check!(f,"spin_commutator_norm",norm(Hs*S-S*Hs),0.;tol=1e-9)
        for g in ground
            v=d.vectors[:,g];check!(f,"S_squared_in_ground_component_$g",dot(v,S*v),.75;tol=1e-9)
            check!(f,"S_squared_variance_component_$g",norm(S*v-.75v),0.;tol=1e-9)
        end
    end
end
const HANDLERS=Dict{String,Function}("A01"=>atomic_eos!,"A02"=>atomic_poles!,"A03"=>(r->dimer!(r,3)),"A04"=>(r->dimer!(r,4)),"A05"=>(r->dimer!(r,5)),"A06"=>free_thermo!,"A07"=>free_propagator!,"A08"=>lindhard!,"A09"=>spectral_poles!,"A10"=>compressibility!,"A11"=>gaussian_levels!,"A12"=>canonical_infinite!,"A13"=>grand_infinite!,"A14"=>polarized!,"A22"=>eta_pairs!,"A23"=>lieb_spin!)
function main()
    output=get(ENV,"ANALYTIC_OUTPUT",joinpath(@__DIR__,"analytic_results.toml"));isfile(output)&&error("Choose a new ANALYTIC_OUTPUT; historical receipts are preserved")
    selected=haskey(ENV,"ANALYTIC_IDS") ? split(ENV["ANALYTIC_IDS"],',') : IMPLEMENTED
    all(id->id in IMPLEMENTED,selected)||error("ANALYTIC_IDS contains unsupported IDs")
    manifest=TOML.parsefile(joinpath(@__DIR__,"manifest.toml"));definitions=Dict(c["id"]=>c for c in manifest["benchmarks"])
    files=["analysis/compute_s2_ed.jl","analysis/audit_u4_l4_thermal_ed.jl"]
    report=Dict{String,Any}("created_utc"=>string(now(UTC)),"host"=>gethostname(),"slurm_job_id"=>get(ENV,"SLURM_JOB_ID","LOCAL"),"julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"source_root"=>ROOT,"methods_exercised"=>["user ED Hamiltonian action and Lanczos","user thermal ED density and reduced-operator routines"],"methods_not_certified"=>["DQMC stochastic/unequal-time APIs","DMRG APIs"],"selected_ids"=>selected,"source_sha256"=>Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in files),"tranche_sha256"=>Dict(p=>bytes2hex(sha256(read(joinpath(@__DIR__,p)))) for p in ["run_analytic_tranche.jl","manifest.toml","PROTOCOL.md"]),"results"=>Any[])
    for n in 1:25
        id="A"*lpad(n,2,'0');definition=definitions[id]
        result=Dict{String,Any}("id"=>id,"title"=>definition["title"],"status"=>"NOT_RUN","fixtures"=>Any[],"source_urls"=>definition["source_urls"])
        if id in selected
            started=time()
            try
                HANDLERS[id](result)
                checks=[c for f in result["fixtures"] for c in f["checks"]]
                result["status"]=isempty(checks) ? "UNRESOLVED" : all(c["status"]=="PASS" for c in checks) ? "PASS" : "FAIL"
                result["check_count"]=length(checks);result["maximum_absolute_error"]=maximum([c["error"] for c in checks if c["metric"]=="maximum_absolute"];init=0.)
            catch ex
                result["status"]="UNRESOLVED";result["reason"]=sprint(showerror,ex,catch_backtrace())
            end
            result["elapsed_seconds"]=time()-started
        else
            result["reason"]=get(definition,"not_run_reason","Not selected for this run")
        end
        push!(report["results"],result);open(output,"w") do io;TOML.print(io,report);end
        println(id," ",result["status"]," fixtures=",length(result["fixtures"]));flush(stdout)
    end
    report["finished_utc"]=string(now(UTC));report["distinct_ideas_passed"]=count(r->r["status"]=="PASS",report["results"])
    open(output,"w") do io;TOML.print(io,report);end
    any(r->r["id"] in selected && r["status"]!="PASS",report["results"])&&exit(1)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end

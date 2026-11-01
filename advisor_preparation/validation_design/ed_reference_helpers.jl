# Pure standard-library ED helpers copied verbatim from run_external_validation.jl.
# Rebuild/check against the canonical consumer when its selected methods change.
function sector_expectation(space,psi,c)
    ob=c["observable"];L=space.L;nd=length(space.dn_masks);ans=0.
    if ob=="S2"
        return ED.renyi2_entropy(psi,space,c["nA"])
    end
    for (iu,u) in enumerate(space.up_masks),(id,d) in enumerate(space.dn_masks)
        k=(iu-1)*nd+id;p=abs2(psi[k]);nu=count_ones(u);nn=count_ones(d);D=count_ones(u&d)
        if ob=="double_occupancy_total"
            ans+=p*D
        elseif ob=="density_total"
            ans+=p*(nu+nn)
        elseif ob=="spin_z_correlation"
            i,j=c["sites"];si=Int((u>>i)&1)-Int((d>>i)&1);sj=Int((u>>j)&1)-Int((d>>j)&1);ans+=p*si*sj/4
        elseif ob=="onebody_up"
            i,j=c["sites"]
            if i==j
                ans+=p*Int((u>>i)&1)
            elseif ((u>>j)&1)==1 && ((u>>i)&1)==0
                mask=u ⊻ (UInt(1)<<i) ⊻ (UInt(1)<<j)
                iu2=findfirst(==(mask),space.up_masks)
                sign=ED.fermion_hop_sign(u,j,i)
                ans+=psi[(iu2-1)*nd+id]*psi[k]*sign
            end
        else
            error("Unsupported sector observable $ob")
        end
    end
    ans
end
function user_ed(c)
    f=cached((:ed,c["L"],c["t"],c["U"],c["boundary"],c["Nup"],c["Ndown"])) do
        # The target solver is uncentered and fixed N. Its mu is set to zero.
        space=ED.build_space(c["L"],c["t"],c["U"],0.,c["Nup"],c["Ndown"];boundary=c["boundary"])
        E,psi,it,reported=ED.ground_state_lanczos(space;tol=1e-11,maxiter=1000)
        hv=similar(psi);ED.apply_hamiltonian!(hv,psi,space);res=norm(hv-E*psi)
        res<=1e-9 || error("Explicit ED eigenpair residual $res exceeds 1e-9")
        (;space,E,psi,it,res)
    end
    v=c["observable"]=="energy_uncentered" ? f.E : sector_expectation(f.space,f.psi,c)
    (;value=v,residual=f.res,iterations=f.it)
end
function fock_observable(rho,masks,L,c)
    ob=c["observable"];d=length(masks);ans=0.
    if ob=="S2"
        return first(THERM.U4L4ThermalEDAudit.renyi2(rho,masks,L,c["nA"]))
    end
    for (iu,u) in enumerate(masks),(id,v) in enumerate(masks)
        k=(iu-1)*d+id;p=real(rho[k,k])
        if ob=="double_occupancy_total"
            ans+=p*count_ones(u&v)
        elseif ob=="density_total"
            ans+=p*(count_ones(u)+count_ones(v))
        elseif ob=="spin_z_correlation"
            i,j=c["sites"];a=Int((u>>i)&1)-Int((v>>i)&1);b=Int((u>>j)&1)-Int((v>>j)&1);ans+=p*a*b/4
        elseif ob=="onebody_up"
            i,j=c["sites"]
            if i==j
                ans+=p*Int((u>>i)&1)
            elseif ((u>>j)&1)==1 && ((u>>i)&1)==0
                u2=u ⊻ (UInt(1)<<i) ⊻ (UInt(1)<<j);iu2=Int(u2)+1
                ans+=real(rho[k,(iu2-1)*d+id])*ED.fermion_hop_sign(u,j,i)
            end
        else
            error("Unsupported thermal observable $ob")
        end
    end
    ans
end
function user_thermal(c)
    f=cached((:thermal,c["L"],c["t"],c["U"],c["mu_centered"],c["beta"],get(c,"dtau",.1),c["boundary"])) do
        a=THERM.U4L4ThermalEDAudit
        K,V,masks=a.full_fock_matrices(c["L"],c["t"],c["U"],c["mu_centered"];boundary=c["boundary"])
        if c["ensemble"]=="grand_canonical_trotter"
            rho,_=a.asymmetric_trotter_density(K,V,c["beta"],c["dtau"])
        else
            rho,_=a.exact_thermal_density(K+Diagonal(V),c["beta"])
        end
        (;K,V,masks,rho)
    end
    ob=c["observable"]
    if ob in ["energy_centered","energy_uncentered"]
        grand=real(tr((f.K+Diagonal(f.V))*f.rho))
        n=fock_observable(f.rho,f.masks,c["L"],merge(c,Dict("observable"=>"density_total")))
        v=grand+c["mu_centered"]*n
        ob=="energy_uncentered" && (v+=c["U"]*n/2-c["U"]*c["L"]/4)
    else
        v=fock_observable(f.rho,f.masks,c["L"],c)
    end
    (;value=v,trace_error=abs(real(tr(f.rho))-1))
end
function analytic_atom(c)
    U=c["U"];mu=c["mu_centered"];beta=c["beta"]
    e=[U/4,-U/4-mu,-U/4-mu,U/4-2mu];w=exp.(-beta.*(e.-minimum(e)));p=w./sum(w)
    N=c["L"]*(p[2]+p[3]+2p[4]);D=c["L"]*p[4]
    ob=c["observable"]
    if ob=="energy_centered"
        return U*(D-N/2+c["L"]/4) # internal centered energy; no -mu*N
    elseif ob=="density_total"
        return N
    elseif ob=="double_occupancy_total"
        return D
    elseif ob=="S2"
        return -c["nA"]*log(sum(p.^2))
    end
    error("Unknown atom observable")
end

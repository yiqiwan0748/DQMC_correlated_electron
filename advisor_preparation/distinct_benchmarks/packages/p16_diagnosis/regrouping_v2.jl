# Additive candidate correction; frozen user's ED source is unchanged.
module P16RegroupingV2
using LinearAlgebra
export fermionic_regrouping_sign, renyi2_entropy_sites_v2
occupied(mask,site)=!iszero(mask & (UInt(1)<<site))
function compressed(mask,sites)
    out=UInt(0)
    for (k,i) in enumerate(sites)
        occupied(mask,i) && (out|=UInt(1)<<(k-1))
    end
    out
end
function fermionic_regrouping_sign(up,dn,A,B)
    # Original modes: up(0:L-1),down(0:L-1).
    # Grouped modes: up(A),down(A),up(B),down(B), each set sorted.
    crossings=count(occupied(up,i)&&occupied(up,j)&&i<j for i in B for j in A)
    crossings+=count(occupied(dn,i)&&occupied(dn,j)&&i<j for i in B for j in A)
    crossings+=count(i->occupied(up,i),B)*count(i->occupied(dn,i),A)
    isodd(crossings) ? -1. : 1.
end
function renyi2_entropy_sites_v2(psi,space,sitesA)
    A=sort(Int.(collect(sitesA)))
    length(unique(A))==length(A) || error("Duplicate subsystem sites")
    all(i->0<=i<space.L,A) || error("Subsystem site out of range")
    B=setdiff(collect(0:space.L-1),A);na=length(A);nb=length(B)
    X=zeros(eltype(psi),1<<(2na),1<<(2nb));nd=length(space.dn_masks)
    for (iu,up) in enumerate(space.up_masks),(id,dn) in enumerate(space.dn_masks)
        a=Int(compressed(up,A)|(compressed(dn,A)<<na))+1
        b=Int(compressed(up,B)|(compressed(dn,B)<<nb))+1
        X[a,b]+=fermionic_regrouping_sign(up,dn,A,B)*psi[(iu-1)*nd+id]
    end
    rho=X*X'
    -log(real(sum(abs2,rho)))
end
end

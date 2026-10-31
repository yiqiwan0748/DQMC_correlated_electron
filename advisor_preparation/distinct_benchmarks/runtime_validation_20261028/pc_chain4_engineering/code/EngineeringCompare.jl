# NEW independent Base-only exact finite-data engineering replay. No CLI,
# imports of target builders, eigensolvers, ground truth, or catalogue credit.
module EngineeringCompare
export run_engineering
const Q=Rational{BigInt}
const Z=Complex{Q}
const MAX_BITS=65536
const C_COEFF=Q(1,BigInt(10)^11)
const C_HERM=Q(1,BigInt(10)^12)
const C_ACTION=Q(2,BigInt(10)^10)
const C_OBSERVER=Q(1,BigInt(10)^15)
const C_OBSERVABLE=Q(1,BigInt(10)^9)
const C_ENTROPY=Q(1,BigInt(10)^8)
const C_NORM=Q(1,BigInt(10)^10)

function bounded(x::Q)
    max(ndigits(abs(numerator(x));base=2),ndigits(denominator(x);base=2))<=MAX_BITS ||
        error("UNRESOLVED_EXACT_INTEGER_GROWTH_BOUND")
    return x
end
bounded(x::Z)=Z(bounded(real(x)),bounded(imag(x)))
function exact(x::Float64)
    isfinite(x) || error("Nonfinite binary64 input")
    u=reinterpret(UInt64,x);negative=(u>>63)==1;e=Int((u>>52)&0x7ff);m=BigInt(u&0x000fffffffffffff)
    e==0 || (m+=BigInt(1)<<52)
    shift=e==0 ? -1074 : e-1023-52
    negative && (m=-m)
    return bounded(shift>=0 ? Q(m<<shift,1) : Q(m,BigInt(1)<<(-shift)))
end
exact(x::Integer)=Q(BigInt(x),1)
qrecord(x::Q)=Dict("numerator_decimal"=>string(numerator(x)),"denominator_decimal"=>string(denominator(x)))
zrecord(x::Z)=Dict("real"=>qrecord(real(x)),"imag"=>qrecord(imag(x)))
function packed_array(p)
    p["order"]=="F" && p["eltype"] in ("Float64","ComplexF64") || error("Raw array encoding differs")
    all(v->v isa Integer && !(v isa Bool),p["shape"]) || error("Noninteger raw shape")
    shape=Int.(p["shape"]);1<=length(shape)<=4 && all(v->1<=v<=256,shape) && prod(shape)==length(p["real"])==length(p["imag"]) || error("Raw array shape differs")
    all(v->typeof(v)===Float64,p["real"]) && all(v->typeof(v)===Float64,p["imag"]) || error("Binary64 raw array required")
    data=Z[Z(exact(a),exact(b)) for (a,b) in zip(p["real"],p["imag"])]
    p["eltype"]=="Float64" && any(!iszero(imag(x)) for x in data) && error("Real tag hides imaginary data")
    return reshape(data,Tuple(shape))
end
function packed_scalar(p)
    p["eltype"] in ("Float64","ComplexF64") || error("Scalar precision differs")
    typeof(p["real"])===Float64 && typeof(p["imag"])===Float64 || error("Binary64 scalar components required")
    value=Z(exact(p["real"]),exact(p["imag"]))
    p["eltype"]=="Float64" && !iszero(imag(value)) && error("Real scalar tag hides imaginary data")
    value
end

# Canonical fermion order is up0..up3,down0..down3. Every sign is
# computed from occupations before that elementary creation/annihilation.
function fermion_word(mask::Int,word)
    phase=1;state=mask
    for (kind,mode) in reverse(word)
        0<=mode<8 && kind in (:create,:annihilate) || error("Bad independent fermion word")
        occupied=((state>>mode)&1)==1
        (kind===:create ? occupied : !occupied) && return (state,0)
        isodd(count_ones(state&((1<<mode)-1))) && (phase=-phase)
        state=xor(state,1<<mode)
    end
    return state,phase
end
const MASKS=[m for m in 0:15 if count_ones(m)==2]
const PAIRS=[(u,d) for u in MASKS for d in MASKS]
const STATES=[u|(d<<4) for (u,d) in PAIRS]
const LOOKUP=Dict(m=>i for (i,m) in enumerate(STATES))
function local_digits(u,d)
    [Int((u>>(p-1))&1)+2Int((d>>(p-1))&1) for p in 1:4]
end
function occupation_phase(u,d)
    parity=sum(Int((d>>i)&1)*count_ones(u>>(i+1)) for i in 0:3)
    return isodd(parity) ? -1 : 1
end
function chain_hamiltonian()
    H=fill(Q(0),36,36)
    for (j,(u,d)) in enumerate(PAIRS)
        H[j,j]+=4count_ones(u&d)
        for species in (0,4),p in 0:2,(a,b) in ((p,p+1),(p+1,p))
            target,sign=fermion_word(STATES[j],[(:create,species+a),(:annihilate,species+b)])
            sign==0 || (H[LOOKUP[target],j]-=sign)
        end
    end
    H
end
norm2(v)=bounded(sum(abs2(x) for x in v))
function local_matrix(name)
    A=fill(Q(0),4,4)
    if name in ("Id","Nup","Ndn","Ntot","Nupdn","Sz","F")
        for m in 0:3
            up=m&1;down=(m>>1)&1
            A[m+1,m+1]=name=="Id" ? Q(1) : name=="Nup" ? Q(up) : name=="Ndn" ? Q(down) :
                name=="Ntot" ? Q(up+down) : name=="Nupdn" ? Q(up*down) : name=="Sz" ? Q(up-down,2) :
                Q(isodd(up+down) ? -1 : 1)
        end
    else
        words=Dict("Cdagup"=>[(:create,0)],"Cup"=>[(:annihilate,0)],
            "Cdagdn"=>[(:create,1)],"Cdn"=>[(:annihilate,1)],
            "S+"=>[(:create,0),(:annihilate,1)],"S-"=>[(:create,1),(:annihilate,0)])
        haskey(words,name) || error("Unregistered local operator")
        for m in 0:3
            target,sign=fermion_word(m,words[name]);sign==0 || (A[target+1,m+1]=sign)
        end
    end
    A
end
function full_grouped(amplitudes)
    result=fill(Z(0),256)
    for index in 0:255
        digits=[(index>>(2p))&3 for p in 0:3]
        up=sum((digits[p]&1)<<(p-1) for p in 1:4)
        down=sum(((digits[p]>>1)&1)<<(p-1) for p in 1:4)
        result[(up|(down<<4))+1]=occupation_phase(up,down)*amplitudes[index+1]
    end
    result
end
function full_word_action(v,word)
    length(v)==256 || error("All256 Fock amplitudes required")
    result=fill(Z(0),256)
    for mask in 0:255
        iszero(v[mask+1]) && continue
        target,sign=fermion_word(mask,word)
        sign==0 || (result[target+1]+=sign*v[mask+1])
    end
    result
end
function full_h_action(v)
    result=fill(Z(0),256)
    for mask in 0:255
        result[mask+1]=4count_ones((mask&15)&(mask>>4))*v[mask+1]
    end
    for species in (0,4),p in 0:2,(i,j) in ((p,p+1),(p+1,p))
        result .-=full_word_action(v,[(:create,species+i),(:annihilate,species+j)])
    end
    result
end
function full_density_values(name,p)
    values=Q[]
    for mask in 0:255
        up=(mask>>(p-1))&1;down=(mask>>(4+p-1))&1
        push!(values,name=="Nup" ? Q(up) : name=="Ndn" ? Q(down) : name=="Ntot" ? Q(up+down) :
            name=="Nupdn" ? Q(up*down) : name=="Sz" ? Q(up-down,2) : error("Unknown full density"))
    end
    values
end
function full_spin_expectations(v)
    zz=fill(Z(0),4,4);pm=copy(zz);mp=copy(zz)
    for p in 1:4,q in 1:4
        z=full_density_values("Sz",p).*full_density_values("Sz",q)
        av=full_word_action(v,[(:create,p-1),(:annihilate,4+p-1),(:create,4+q-1),(:annihilate,q-1)])
        bv=full_word_action(v,[(:create,4+p-1),(:annihilate,p-1),(:create,q-1),(:annihilate,4+q-1)])
        zz[p,q]=sum(abs2(v[j])*z[j] for j in eachindex(v))
        pm[p,q]=sum(conj(v[j])*av[j] for j in eachindex(v))
        mp[p,q]=sum(conj(v[j])*bv[j] for j in eachindex(v))
    end
    return zz,pm,mp,zz+(pm+mp)/2
end

function full_tensor_amplitudes(state)
    sites=state["sites"];length(sites)==4 || error("Exactly four native MPS tensors required")
    arrays=[packed_array(v["tensor"]) for v in sites]
    for p in 1:4
        ndims(arrays[p])==3 && size(arrays[p],3)==4 || error("MPS local axes differ")
        sites[p]["axes"]==["left","right","physical"] || error("MPS axes metadata differ")
        p==1 ? size(arrays[p],1)==1 || error("Left boundary differs") :
            size(arrays[p-1],2)==size(arrays[p],1) || error("MPS link dimensions differ")
        p==4 && size(arrays[p],2)!=1 && error("Right boundary differs")
    end
    full=fill(Z(0),256)
    for index in 0:255
        vector=[Z(1)]
        for p in 1:4
            physical=((index>>(2(p-1)))&3)+1;A=arrays[p];out=fill(Z(0),size(A,2))
            for l in eachindex(vector)
                iszero(vector[l]) && continue
                for r in eachindex(out)
                    iszero(A[l,r,physical]) || (out[r]+=vector[l]*A[l,r,physical])
                end
            end
            vector=out
        end
        full[index+1]=bounded(only(vector))
    end
    full
end
function mpo_selected_matrix(packed)
    packed["axes"]==["left","right","bra","ket"] || error("MPO axes differ")
    arrays=[packed_array(v) for v in packed["sites"]];length(arrays)==4 || error("Four MPO tensors required")
    for p in 1:4
        ndims(arrays[p])==4 && size(arrays[p],3)==size(arrays[p],4)==4 || error("MPO physical axes differ")
        p==1 ? size(arrays[p],1)==1 || error("MPO left boundary") :
            size(arrays[p-1],2)==size(arrays[p],1) || error("MPO link dimensions")
        p==4 && size(arrays[p],2)!=1 && error("MPO right boundary")
    end
    H=fill(Z(0),36,36)
    for (i,(ui,di)) in enumerate(PAIRS),(j,(uj,dj)) in enumerate(PAIRS)
        bra=local_digits(ui,di);ket=local_digits(uj,dj);vector=[Z(1)]
        for p in 1:4
            A=arrays[p];out=fill(Z(0),size(A,2))
            for l in eachindex(vector)
                iszero(vector[l]) && continue
                for r in eachindex(out)
                    a=A[l,r,bra[p]+1,ket[p]+1];iszero(a) || (out[r]+=vector[l]*a)
                end
            end
            vector=out
        end
        H[i,j]=bounded(occupation_phase(ui,di)*occupation_phase(uj,dj)*only(vector))
    end
    H
end

# For q in[1/2,2], z=(q-1)/(q+1) has |z|<=1/3. The16-term
# atanh series tail is at most 2|z|^33/(33(1-z^2)); every operation
# below is exact rational arithmetic. No libm or MPFR log is assumed.
function unit_log(q::Q)
    Q(1,2)<=q<=2 || error("Log reduction domain")
    z=(q-1)/(q+1);term=z;value=Q(0)
    for k in 0:15
        value=bounded(value+2term/(2k+1));term=bounded(term*z*z)
    end
    tail=bounded(2*abs(z)^33/(33(1-z*z)))
    return (lo=value-tail,hi=value+tail)
end
function negative_log(q::Q)
    q>0 || error("Nonpositive exact purity")
    y=q;power=0
    while y<Q(1,2);y*=2;power-=1;abs(power)<=16 || error("Operational log scaling cap");end
    while y>1;y/=2;power+=1;abs(power)<=16 || error("Operational log scaling cap");end
    a=unit_log(y);b=unit_log(Q(2))
    lo=a.lo+(power>=0 ? power*b.lo : power*b.hi)
    hi=a.hi+(power>=0 ? power*b.hi : power*b.lo)
    return (lo=bounded(-hi),hi=bounded(-lo),power=power)
end
function exact_purity(amplitudes,cut)
    cut in (1,2,3) && length(amplitudes)==256 || error("Purity dimension/cut")
    rows=4^cut;cols=4^(4-cut);A=reshape(amplitudes,rows,cols);norm=norm2(amplitudes);norm>0 || error("Zero state")
    purity=Q(0)
    if rows<=cols
        for a in 1:rows,b in 1:rows
            value=sum(A[a,j]*conj(A[b,j]) for j in 1:cols);purity+=abs2(value)
        end
    else
        for a in 1:cols,b in 1:cols
            value=sum(conj(A[i,a])*A[i,b] for i in 1:rows);purity+=abs2(value)
        end
    end
    p=bounded(purity/norm^2)
    Q(1,min(rows,cols))<=p<=1 || error("Exact pure-state RDM rank/purity inequality failed")
    return p
end

mutable struct Audit
    checks::Vector{Any}
    on_record::Function
end
function check!(a::Audit,name,passed;details=Dict{String,Any}())
    record=Dict("name"=>name,"passed"=>Bool(passed),"details"=>details)
    push!(a.checks,record)
    a.on_record("ENGINEERING_CHECK",deepcopy(record))===true || error("Comparison checkpoint failed")
    passed || error("ENGINEERING_REPLAY_REJECTED: "*name)
    nothing
end
function close_data!(a,name,observed,expected,tolerance::Q)
    same=size(observed)==size(expected)
    check!(a,name*" dimensions",same)
    residual=Q(0);worst=0
    for (i,(x,y)) in enumerate(zip(observed,expected))
        delta=abs2(x-y)
        if delta>residual;residual=delta;worst=i;end
    end
    check!(a,name,residual<=tolerance^2;details=Dict("maximum_squared_residual"=>qrecord(residual),
        "absolute_gate"=>qrecord(tolerance),"worst_column_major_index"=>worst,
        "evaluation"=>"Exact rationals on recorded binary64 data; no floating replay uncertainty"))
    return residual
end
close_scalar!(a,name,x,y,t)=close_data!(a,name,[Z(x)],[Z(y)],t)
function state_data!(a,label,state;unit_norm=true)
    check!(a,label*" raw state conventions",state["M_spatial_modes"]==4 && state["spin_orbitals"]==8 &&
        state["complete_local_array_shape"]==[4,4,4,4] && state["external_norm_scalar"]==1.0 &&
        state["normalized_or_phase_fitted"]===false && state["local_physical_order"]==["Emp","Up","Dn","UpDn"])
    check!(a,label*" actual QN physical indices retained",length(state["sites"])==4 &&
        all(s["physical_index_hasqns"]===true && s["left_index_hasqns"] isa Bool &&
            s["right_index_hasqns"] isa Bool for s in state["sites"]))
    amplitudes=vec(packed_array(state["complete_local_amplitudes"]))
    check!(a,label*" all256 amplitudes and supported real path",length(amplitudes)==256 && all(iszero(imag(x)) for x in amplitudes))
    reconstructed=full_tensor_amplitudes(state)
    close_data!(a,label*" exact contraction of actual MPS tensors",amplitudes,reconstructed,C_NORM)
    association_error=norm2(amplitudes.-reconstructed)
    check!(a,label*" complete tensor-to-export L2 association",association_error<=C_NORM^2;
        details=Dict("squared_L2_error"=>qrecord(association_error),"absolute_L2_gate"=>qrecord(C_NORM)))
    rows=state["rows"]
    check!(a,label*" complete36-row ordered sector census",length(rows)==36 && [(r["up_mask"],r["down_mask"]) for r in rows]==PAIRS)
    grouped=Z[]
    for (r,(u,d)) in zip(rows,PAIRS)
        digits=local_digits(u,d);phase=occupation_phase(u,d);index=1+sum(digits[p]*4^(p-1) for p in 1:4)
        check!(a,label*" fixed occupation/parity "*string(u)*":"*string(d),r["local_labels"]==digits.+1 && r["grouped_spin_sign"]==phase)
        value=Z(exact(r["real"]),exact(r["imag"]));push!(grouped,value)
        close_scalar!(a,label*" actual local-to-grouped coefficient "*string(u)*":"*string(d),value,phase*amplitudes[index],Q(0))
    end
    fullnorm=norm2(amplitudes);selectednorm=norm2(grouped)
    check!(a,label*" positive exact raw state norm",fullnorm>0)
    outside=Q(0);tensor_outside=Q(0)
    for index in 0:255
        digits=[(index>>(2p))&3 for p in 0:3]
        up=sum((digits[p]&1)<<(p-1) for p in 1:4);down=sum(((digits[p]>>1)&1)<<(p-1) for p in 1:4)
        if count_ones(up)!=2 || count_ones(down)!=2
            outside+=abs2(amplitudes[index+1]);tensor_outside+=abs2(reconstructed[index+1])
        end
    end
    # Original protocol left norm/sector units implicit. Before any numerical
    # data, conservatively require vector L2<=C_NORM via squared comparisons.
    check!(a,label*" full256 outside-sector L2 gate (sqrt-free)",outside<=C_NORM^2 && tensor_outside<=C_NORM^2;
        details=Dict("exported_point_outside_sector_norm_squared"=>qrecord(outside),
            "exact_tensor_outside_sector_norm_squared"=>qrecord(tensor_outside),"full_norm_squared"=>qrecord(fullnorm),
            "outside_sector_absolute_L2_gate"=>qrecord(C_NORM),
            "outside_sector_squared_L2_gate"=>qrecord(C_NORM^2),
            "gate_scope"=>"Conservative pre-data engineering strengthening; original protocol units were implicit"))
    close_scalar!(a,label*" exported complete norm",exact(state["complete_local_norm2"]),fullnorm,C_NORM)
    close_scalar!(a,label*" exported selected norm",exact(state["selected_grouped_norm2"]),selectednorm,C_NORM)
    close_scalar!(a,label*" selected/completed norm agreement",selectednorm,fullnorm,C_NORM)
    reconstructed_norm=norm2(reconstructed)
    unit_norm && close_scalar!(a,label*" unit norm",reconstructed_norm,Q(1),C_NORM)
    a.on_record("COMPARISON_STATE",Dict("label"=>label,"status"=>"EXACT_RAW_STATE_CONTRACTION_AND_SECTOR_REPLAY_COMPLETE",
        "full_norm_squared"=>qrecord(fullnorm),"exported_point_outside_sector_norm_squared"=>qrecord(outside),
        "exact_tensor_outside_sector_norm_squared"=>qrecord(tensor_outside)))===true || error("State comparison checkpoint failed")
    return (amplitudes=reconstructed,exported_amplitudes=amplitudes,grouped=grouped,
        full_grouped=full_grouped(reconstructed),point_full_grouped=full_grouped(amplitudes),norm=reconstructed_norm)
end
function entropy_check!(a,label,amplitudes,cut,observed;schmidt=nothing,normcheck=nothing)
    purity=exact_purity(amplitudes,cut);ball=negative_log(purity)
    check!(a,label*" cut"*string(cut)*" finite log enclosure",ball.hi-ball.lo<=C_HERM;
        details=Dict("purity"=>qrecord(purity),"S2_lower"=>qrecord(ball.lo),"S2_upper"=>qrecord(ball.hi),
            "atanh_terms"=>16,"powers_of_two_reduction"=>ball.power,"rounding_model"=>"Exact rational series and absolute tail bound"))
    value=exact(observed)
    check!(a,label*" cut"*string(cut)*" native entropy",max(abs(value-ball.lo),abs(value-ball.hi))<=C_ENTROPY;
        details=Dict("native_binary64_exact"=>qrecord(value),"absolute_gate"=>qrecord(C_ENTROPY)))
    if normcheck!==nothing
        close_scalar!(a,label*" native entropy raw norm cut"*string(cut),exact(normcheck),norm2(amplitudes),C_NORM)
    end
    if schmidt!==nothing
        singular=exact.(schmidt["singular_values"])
        check!(a,label*" nonnegative raw Schmidt values cut"*string(cut),!isempty(singular) && all(v->v>=0,singular) && schmidt["active_MPS_mutated"]===false)
        snorm=sum(x^2 for x in singular);spurity=sum(x^4 for x in singular)
        close_scalar!(a,label*" Schmidt norm cut"*string(cut),snorm,norm2(amplitudes),C_NORM)
        close_scalar!(a,label*" raw SVD norm point cut"*string(cut),exact(schmidt["raw_norm2"]),snorm,C_NORM)
        close_scalar!(a,label*" raw SVD purity point cut"*string(cut),exact(schmidt["raw_purity_numerator"]),spurity,C_NORM)
        close_scalar!(a,label*" raw SVD explicit denominator cut"*string(cut),exact(schmidt["normalization_denominator"]),snorm^2,C_NORM)
        snorm>0 && spurity>0 || error("Nonpositive exact Schmidt norm/purity")
        p=spurity/snorm^2
        close_scalar!(a,label*" normalized SVD purity cut"*string(cut),p,purity,C_NORM)
        close_scalar!(a,label*" explicit SVD purity point cut"*string(cut),exact(schmidt["explicit_ratio_purity"]),p,C_NORM)
        sball=negative_log(p);svalue=exact(schmidt["S2_from_explicit_ratio"])
        check!(a,label*" raw SVD natural log point cut"*string(cut),max(abs(svalue-sball.lo),abs(svalue-sball.hi))<=C_ENTROPY)
    end
    a.on_record("COMPARISON_CUT",Dict("label"=>label,"prefix_cut"=>cut,"status"=>"EXACT_PURITY_AND_FINITE_LOG_ENCLOSURE_COMPLETE",
        "purity"=>qrecord(purity),"S2_lower"=>qrecord(ball.lo),"S2_upper"=>qrecord(ball.hi)))===true || error("Cut comparison checkpoint failed")
    return purity
end
function observer_checks!(a,raw)
    trace=raw["full_bond_trace"];expected=[(s,h,b) for s in 1:12 for h in 1:2 for b in (h==1 ? (1:3) : (3:-1:1))]
    check!(a,"Actual72 ordered DMRG measure callbacks",length(trace)==72 && [(r["sweep"],r["half_sweep"],r["bond"]) for r in trace]==expected && raw["actual_bond_callback_records"]==trace)
    stops=raw["checkdone_records"]
    check!(a,"All12 actual original checkdone callbacks without early stop",length(stops)==12 && [v["sweep"] for v in stops]==collect(1:12) && all(v["inner_requested_stop"]===false for v in stops))
    noise=[1e-6,1e-7,1e-8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for (i,r) in enumerate(trace)
        check!(a,"Actual bond controls "*string(i),r["actual_sweeps_maxdim"]==2 &&
            r["actual_sweeps_mindim"]==1 && exact(r["actual_sweeps_cutoff"])==exact(1e-4) &&
            exact(r["actual_sweeps_noise"])==exact(noise[r["sweep"]]) &&
            1<=r["retained_linkdim"]<=2 && exact(r["truncerror"])>=0 &&
            r["sweep_is_done"]==((r["half_sweep"],r["bond"])==(2,1)))
    end
    half2=[maximum(exact(r["truncerror"]) for r in trace if r["sweep"]==s && r["half_sweep"]==2) for s in 1:12]
    both=[maximum(exact(r["truncerror"]) for r in trace if r["sweep"]==s) for s in 1:12]
    for stage in raw["stages"]
        s=stage["stage_sweep"];obs=stage["observer"]
        close_data!(a,"Original half_sweep2 maxima atstage"*string(s),exact.(obs["native_half_sweep2_maxima"]),half2[1:s],C_OBSERVER)
        close_data!(a,"Independently replayed half_sweep2 maxima atstage"*string(s),exact.(obs["replayed_half_sweep2_maxima"]),half2[1:s],C_OBSERVER)
        close_data!(a,"Separately named both-direction maxima atstage"*string(s),exact.(obs["NEW_adapter_both_direction_maxima"]),both[1:s],C_OBSERVER)
        close_scalar!(a,"Native final half_sweep2 weight stage"*string(s),exact(obs["native_final_half_sweep2_maximum"]),half2[s],C_OBSERVER)
        check!(a,"Source-native observer discrepancy stage"*string(s),exact(obs["native_trace_max_abs_discrepancy"])<=C_OBSERVER)
        window=obs["zero_noise_native_weight_window"]
        if s==3
            check!(a,"Empty stage3 zero-noise window remains null/notapplicable",window==Dict("status"=>"NOT_APPLICABLE_NO_ZERO_NOISE_SWEEPS","value"=>nothing,"count"=>0,"first_sweep"=>nothing,"last_sweep"=>nothing))
        else
            check!(a,"Nonempty zero-noise window metadata stage"*string(s),window["status"]=="OBSERVED_ZERO_NOISE_SWEEPS_4_TO_STAGE" && window["count"]==s-3 && window["first_sweep"]==4 && window["last_sweep"]==s)
            close_scalar!(a,"Nonempty zero-noise maximum stage"*string(s),exact(window["value"]),maximum(half2[4:s]),C_OBSERVER)
        end
    end
    energies=raw["native_energies_each_completed_sweep"]
    check!(a,"All12 original optimizer energies retained",length(energies)==12)
    for s in 1:12
        last_row=only([r for r in trace if r["sweep"]==s && r["sweep_is_done"]])
        close_scalar!(a,"Original completed-sweep optimizer energy "*string(s),packed_scalar(energies[s]),packed_scalar(last_row["reported_local_optimizer_energy"]),C_OBSERVER)
    end
end

function observable_checks!(a,label,stage,data)
    v=data.full_grouped;n=norm2(v);hv=full_h_action(v)
    energy=sum(conj(v[i])*hv[i] for i in eachindex(v))
    close_scalar!(a,label*" raw MPS norm contraction",exact(stage["raw_norm2"]),n,C_NORM)
    E=stage["energy"]
    close_scalar!(a,label*" raw Hamiltonian numerator",packed_scalar(E["raw_numerator"]),energy,C_OBSERVABLE)
    close_scalar!(a,label*" explicit energy denominator",exact(E["norm2_denominator"]),n,C_NORM)
    close_scalar!(a,label*" normalized energy",packed_scalar(E["explicit_ratio"]),energy/n,C_OBSERVABLE)
    point=data.point_full_grouped;pn=norm2(point);phv=full_h_action(point)
    close_scalar!(a,label*" same exported-point state energy",packed_scalar(E["explicit_ratio"]),sum(conj(point[i])*phv[i] for i in eachindex(point))/pn,C_OBSERVABLE)
    for name in ("Nup","Ndn","Ntot","Nupdn","Sz")
        reference=Z[]
        for p in 1:4
            diagonal=full_density_values(name,p)
            push!(reference,Z(sum(abs2(v[i])*diagonal[i] for i in eachindex(v))))
        end
        local_record=stage["local_observables"][name]
        close_data!(a,label*" raw local "*name,packed_array(local_record["raw_numerator"]),reference,C_OBSERVABLE)
        close_scalar!(a,label*" local denominator "*name,exact(local_record["norm2_denominator"]),n,C_NORM)
        close_data!(a,label*" normalized local "*name,packed_array(local_record["explicit_ratio"]),reference./n,C_OBSERVABLE)
    end
    zz,pm,mp,spin=full_spin_expectations(v)
    close_data!(a,label*" raw full physical spin-dot",packed_array(stage["spin_dot"]["raw_numerator"]),spin,C_OBSERVABLE)
    close_scalar!(a,label*" spin denominator",exact(stage["spin_dot"]["norm2_denominator"]),n,C_NORM)
    close_data!(a,label*" normalized full physical spin-dot",packed_array(stage["spin_dot"]["explicit_ratio"]),spin./n,C_OBSERVABLE)
    native=stage["native_correlation_API"]
    for (key,reference) in (("SzSz",zz),("SpSm",pm),("SmSp",mp),("spin_dot",spin))
        close_data!(a,label*" actual correlation_matrix normalized "*key,packed_array(native[key]),reference./n,C_OBSERVABLE)
    end
    entropy=stage["S2"]
    check!(a,label*" all three registered prefix cuts",[x["prefix_cut"] for x in entropy]==[1,2,3])
    for e in entropy
        b=e["prefix_cut"]
        entropy_check!(a,label*" exact-tensor",data.amplitudes,b,e["native_S2"];
            schmidt=e["NEW_raw_SVD_record"],normcheck=e["native_normcheck"])
        entropy_check!(a,label*" exact-exported-point",data.exported_amplitudes,b,e["native_S2"])
    end
end
function witness_word(term)
    length(term)==5 && exact(term[1])==1 || error("Unit witness term convention differs")
    create,p,destroy,q=term[2:5]
    p isa Integer && q isa Integer && 1<=p<=4 && 1<=q<=4 || error("Witness physical site differs")
    if (create,destroy)==("Cdagup","Cup")
        return [(:create,p-1),(:annihilate,q-1)],p,q,"up"
    elseif (create,destroy)==("Cdagdn","Cdn")
        return [(:create,4+p-1),(:annihilate,4+q-1)],p,q,"down"
    elseif (create,destroy)==("S+","S-")
        return [(:create,p-1),(:annihilate,4+p-1),(:create,4+q-1),(:annihilate,q-1)],p,q,"spin"
    end
    error("Unknown fixed engineering witness operator")
end
function apply_trace_checks!(a,label,action)
    trace=action["truncation_trace"]
    check!(a,label*" source-supported actual apply trace",action["alg"]=="naive" && exact(action["cutoff"])==0 &&
        action["maxdim"]==256 && action["explicit_normalize_called"]===false &&
        [(r["link_from"],r["link_to"]) for r in trace]==[(4,3),(3,2),(2,1)] &&
        all(exact(r["truncation_error"])>=0 for r in trace))
    dims=action["retained_link_dimensions"]
    check!(a,label*" actual retained apply link dimensions",length(dims)==3 && all(v->v isa Integer && 1<=v<=256,dims) &&
        dims==[s["tensor"]["shape"][2] for s in action["state"]["sites"][1:3]])
end
function source_regressions!(a)
    cases=[(UInt64(0),Q(0)),(UInt64(0x8000000000000000),Q(0)),
        (UInt64(0x3ff0000000000000),Q(1)),(UInt64(0xbff0000000000000),Q(-1)),
        (UInt64(0x3ff8000000000000),Q(3,2)),(UInt64(1),Q(1,BigInt(1)<<1074)),
        (UInt64(0x0010000000000000),Q(1,BigInt(1)<<1022)),
        (UInt64(0x7fefffffffffffff),Q(((BigInt(1)<<53)-1)<<971,1))]
    check!(a,"Actual binary64-to-dyadic helper finite/subnormal/max/zero regressions",
        all(exact(reinterpret(Float64,bits))==wanted for (bits,wanted) in cases))
    rejected=Bool[]
    for bits in (UInt64(0x7ff0000000000000),UInt64(0x7ff8000000000000))
        try;exact(reinterpret(Float64,bits));push!(rejected,false);catch;push!(rejected,true);end
    end
    check!(a,"Actual dyadic helper rejects infinity and NaN",all(rejected))
    car=true;cc=true
    for mask in 0:255,i in 0:7,j in 0:7
        terms=Dict{Int,Int}()
        for word in ([(:annihilate,i),(:create,j)],[(:create,j),(:annihilate,i)])
            target,sign=fermion_word(mask,word);terms[target]=get(terms,target,0)+sign
        end
        filter!(p->last(p)!=0,terms)
        car &= terms==(i==j ? Dict(mask=>1) : Dict{Int,Int}())
        terms=Dict{Int,Int}()
        for word in ([(:create,i),(:create,j)],[(:create,j),(:create,i)])
            target,sign=fermion_word(mask,word);terms[target]=get(terms,target,0)+sign
        end
        filter!(p->last(p)!=0,terms);cc &= isempty(terms)
    end
    check!(a,"Actual independent fermion helper satisfies all256x8x8 canonical mixed anticommutators",car)
    check!(a,"Actual independent fermion helper satisfies all256x8x8 creation anticommutators",cc)
    phases=true
    for up in 0:15,down in 0:15
        local_order=Int[]
        for p in 0:3
            ((up>>p)&1)==1 && push!(local_order,p)
            ((down>>p)&1)==1 && push!(local_order,4+p)
        end
        inversions=0
        for i in eachindex(local_order),j in (i+1):length(local_order)
            local_order[i]>local_order[j] && (inversions+=1)
        end
        phases &= occupation_phase(up,down)==(isodd(inversions) ? -1 : 1)
    end
    check!(a,"All256 occupation signs equal independently counted creation-order inversions",phases)
    b=negative_log(Q(1,16));check!(a,"Actual rational16-term log enclosure has adequate pure-state precision",unit_log(Q(1))==(lo=Q(0),hi=Q(0)) && b.lo>0 && b.hi-b.lo<C_HERM)
    overflow_rejected=false
    try;bounded(Q(BigInt(1)<<MAX_BITS,1));catch;overflow_rejected=true;end
    check!(a,"Actual exact-integer growth guard refuses oversize arithmetic",overflow_rejected)
end

function run_engineering(raw;on_record::Function)
    audit=Audit(Any[],on_record)
    source_regressions!(audit)
    check!(audit,"Engineering schema and zero catalogue credit",raw["schema"]=="X08_NATIVE_ARM_ENGINEERING_RAW_V1" && raw["scientific_credit"]==0 && raw["catalogue_id"]=="X08")
    model=raw["model"];arm=raw["arm"]
    check!(audit,"Fixed registered model",model==Dict("id"=>"chain4_U4","M"=>4,"N"=>4,"Nup"=>2,"Ndown"=>2,"t"=>1.0,"U"=>4.0,"mu"=>0.0,"boundary"=>"obc","cuts"=>[1,2,3]))
    check!(audit,"Fixed registered engineering arm",arm["chi_cap"]==2 && exact(arm["cutoff"])==exact(1e-4) && arm["start"]=="A" && arm["seed"]==698071 &&
        arm["initial_labels"]==["Up","Dn","Up","Dn"] && arm["initial_linkdims"]==2 && arm["nsweeps"]==12 && arm["snapshot_sweeps"]==[3,6,12] &&
        arm["maxdim"]==fill(2,12) && arm["noise"]==[1e-6,1e-7,1e-8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] &&
        exact(arm["eigsolve_tol"])==exact(1e-12) && arm["eigsolve_krylovdim"]==20 && arm["eigsolve_maxiter"]==20 && arm["fresh_start"]===true && arm["warm_start_other_arm"]===false)
    check!(audit,"All registered stages retained",[s["stage_sweep"] for s in raw["stages"]]==[3,6,12])
    K=fill(Z(0),4,4)
    for i in 1:3;K[i,i+1]=K[i+1,i]=-1;end
    close_data!(audit,"Independent nearest-neighbor kinetic coefficients",packed_array(raw["native_K"]),K,C_COEFF)
    operators=("Id","Nup","Ndn","Ntot","Nupdn","Sz","S+","S-","Cdagup","Cup","Cdagdn","Cdn","F")
    check!(audit,"Four actual local operator tables",length(raw["native_local_operators"])==4)
    for (p,table) in enumerate(raw["native_local_operators"])
        check!(audit,"Exact13 local operator names site"*string(p),Set(keys(table))==Set(operators))
        for name in operators
            close_data!(audit,"Independent two-mode Fock operator "*name*" site"*string(p),packed_array(table[name]),local_matrix(name),C_COEFF)
        end
    end
    H=chain_hamiltonian();MPO=mpo_selected_matrix(raw["native_MPO"])
    close_data!(audit,"All36x36 selected-sector native MPO coefficients",MPO,H,C_COEFF)
    close_data!(audit,"Exact native MPO Hermiticity on selected sector",MPO,conj.(permutedims(MPO)),C_HERM)
    initial=state_data!(audit,"initial",raw["initial_state"])
    observer_checks!(audit,raw)
    witnesses=raw["NEW_mapping_operator_witnesses"]
    check!(audit,"All independent sign diagnostic names",length(witnesses)==3 && Set(w["name"] for w in witnesses)==Set(["nonlocal_odd_parity_up","nonlocal_odd_parity_down","transverse_spin_exchange"]))
    for w in witnesses
        label="witness "*w["name"];input=state_data!(audit,label*" input",w["input"]);action=state_data!(audit,label*" action",w["action"];unit_norm=false)
        apply_trace_checks!(audit,label,w["actual_apply_diagnostics"])
        check!(audit,label*" action diagnostic state join",w["actual_apply_diagnostics"]["state"]==w["action"])
        word,p,q,species=witness_word(w["operator_term"])
        nonzero=findall(!iszero,input.full_grouped)
        check!(audit,label*" unit grouped product input",length(nonzero)==1 && input.full_grouped[only(nonzero)]==1)
        expected=full_word_action(input.full_grouped,word)
        close_data!(audit,label*" independent fermionic action",action.full_grouped,expected,C_ACTION)
        if species!="spin"
            mask=only(nonzero)-1;occupied=species=="up" ? (mask&15) : (mask>>4)
            middle=sum((occupied>>(k-1))&1 for k in (min(p,q)+1):(max(p,q)-1))
            check!(audit,label*" declared odd nonlocal parity",abs(p-q)>1 && isodd(middle) &&
                exact(w["expected_unit_operator_grouped_entry"])==-1 &&
                w["support_pair_has_nonzero_Hamiltonian_hopping"]===false && exact(w["native_H_coefficient_at_pair"])==0)
            target=w["expected_grouped_target_up_mask"]|(w["expected_grouped_target_down_mask"]<<4)
            check!(audit,label*" declared target unit coefficient",expected[target+1]==-1 && count(!iszero,expected)==1)
        end
    end
    stage_data=Dict{Int,Any}()
    for stage in raw["stages"]
        s=stage["stage_sweep"];label="sweep"*string(s);data=state_data!(audit,label,stage["state"]);stage_data[s]=data
        check!(audit,label*" snapshot copy measurement policy",stage["snapshot_was_deepcopy"]===true && stage["optimizer_MPS_mutated_by_measurement"]===false)
        observable_checks!(audit,label,stage,data)
        action=stage["native_MPO_action"]
        apply_trace_checks!(audit,label,action)
        check!(audit,label*" declared actual apply controls",action["alg"]=="naive" && exact(action["cutoff"])==0 && action["maxdim"]==256 && action["explicit_normalize_called"]===false)
        got=state_data!(audit,label*" actual MPO action",action["state"];unit_norm=false);wanted=full_h_action(data.full_grouped)
        close_data!(audit,label*" full256 independent Hamiltonian state action",got.full_grouped,wanted,C_ACTION)
        residual=norm2(got.full_grouped.-wanted);denominator=norm2(wanted)
        if denominator==0
            check!(audit,label*" zero reference-action relative policy",residual==0;details=Dict("relative_status"=>"NOT_APPLICABLE_ZERO_REFERENCE_NORM","squared_residual"=>qrecord(residual)))
        else
            check!(audit,label*" independent relative action2norm",residual<=C_ACTION^2*denominator;
                details=Dict("squared_residual"=>qrecord(residual),"reference_action_norm_squared"=>qrecord(denominator),"relative_gate"=>qrecord(C_ACTION)))
        end
    end
    close_scalar!(audit,"Returned optimizer energy joins final original observer value",
        packed_scalar(raw["final_reported_optimizer_energy"]),packed_scalar(last(raw["native_energies_each_completed_sweep"])),C_OBSERVER)
    close_scalar!(audit,"Returned final MPS raw norm joins exact final-stage tensors",
        packed_scalar(raw["final_raw_norm2"]),stage_data[12].norm,C_NORM)
    replays=raw["checkpoint_replays"]
    check!(audit,"All3 binary checkpoint replays retained",[v["stage_sweep"] for v in replays]==[3,6,12])
    for r in replays
        s=r["stage_sweep"];data=state_data!(audit,"checkpoint replay"*string(s),r["state"])
        check!(audit,"Checkpoint descriptor present"*string(s),r["checkpoint_bytes"]>0 && length(r["checkpoint_sha256"])==64)
        close_data!(audit,"Checkpoint actual MPS state association"*string(s),data.amplitudes,stage_data[s].amplitudes,C_NORM)
    end
    scaled=raw["NEW_explicit_norm_scaling_witness"]
    check!(audit,"Explicit two-times tensor norm diagnostic",scaled["tensor_multiplier"]==2.0 && scaled["source_snapshot_sweep"]==3)
    scaled_data=state_data!(audit,"scaled norm diagnostic",scaled["state"];unit_norm=false)
    close_data!(audit,"Actual scaled MPS represents twice the source snapshot",scaled_data.amplitudes,2 .*stage_data[3].amplitudes,C_NORM)
    close_scalar!(audit,"Actual scaled raw norm is four times",scaled_data.norm,4*stage_data[3].norm,C_NORM)
    zz,pm,mp,spin=full_spin_expectations(scaled_data.full_grouped)
    for (key,wanted) in (("Czz",zz),("Cpm",pm),("Cmp",mp),("normalized_spin_dot",spin))
        close_data!(audit,"Source-expected normalized correlation API on nonunit norm "*key,packed_array(scaled[key]),wanted./scaled_data.norm,C_OBSERVABLE)
    end
    check!(audit,"All scaled diagnostic prefix entropies",[e["prefix_cut"] for e in scaled["S2"]]==[1,2,3])
    for e in scaled["S2"]
        entropy_check!(audit,"scaled norm diagnostic",scaled_data.amplitudes,e["prefix_cut"],e["native_S2"];normcheck=e["native_normcheck"])
    end
    result=Dict("status"=>"PASS_EXACT_FINITE_DATA_ENGINEERING_REPLAY_NO_X08_CREDIT","scientific_credit"=>0,
        "catalogue_status"=>"NOT_RUN","check_count"=>length(audit.checks),"failed_checks"=>0,
        "reference"=>"Independently authored exact eight-mode fermion algebra for fixed chain4 Hamiltonian; no native K/MPO used as truth input",
        "arithmetic"=>"Recorded binary64 converted to exact dyadic rationals; all tensor, operator, norm, expectation and purity arithmetic exact; natural logs enclosed by rational16-term atanh series",
        "operational_max_state_norm_purity_and_log_rational_bits"=>MAX_BITS,
        "scope"=>"Only this excluded small numerical engineering arm; no global fixed-N ground-state reference, solver accuracy, discarded-weight calibration, full6arm pilot, or X08 completion")
    on_record("COMPARISON_COMPLETE",deepcopy(result))===true || error("Final comparison checkpoint failed")
    result
end
end # module EngineeringCompare

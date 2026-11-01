#!/usr/bin/env python3
"""Generate additive P21 overlays in a NEW directory; never modify input trees.

ALF-derived source retains its original GPL notices. This observer modification
is not an official ALF release. No simulation, compilation or submission here.
"""
import argparse
import difflib
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PIN = "f4c500cc5d7bcf089400c63e4e75382c558a4ed3"
USER_SHA = "b72bbbf158381b0f776c1fa259d532b069094c47e361d3490c19ea0f43ada4be"

def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError("Patch anchor not unique: " + old[:100])
    return text.replace(old, new, 1)

def alf_hubbard(text):
    text = once(text, "          Allocate ( Obs_scal(4) )", """          ! P21 additive observer; original ALF sampling/propagation is unchanged.
          ! This modified source is not an official ALF version.
          If (Ltau == 1) then
             If (Projector .or. Symm .or. Checkerboard) error stop 'P21 representation mismatch'
             If (.not.Mz .or. .not.Continuous) error stop 'P21 requires exact continuous Mz HS'
             If (Ndim /= 4 .or. N_FL /= 2 .or. N_SUN /= 1 .or. Ltrot /= 10) &
                  & error stop 'P21 fixture dimensions mismatch'
             If (abs(Dtau-0.1d0)>1.d-14 .or. abs(Beta-1.d0)>1.d-14 .or. &
                  & abs(Ham_U-2.d0)>1.d-14 .or. abs(Ham_T-1.d0)>1.d-14 .or. &
                  & abs(Ham_chem)>1.d-14 .or. abs(ham_h0)>1.d-14 .or. &
                  & abs(Phi_X)>1.d-14 .or. abs(Phi_Y)>1.d-14 .or. N_Phi/=0) &
                  & error stop 'P21 interacting fixture mismatch'
             Allocate ( Obs_scal(7) )
          else
             Allocate ( Obs_scal(4) )
          endif""")
    text = once(text, 'N = 1;   Filename = "Ener"', '''N = 1;   Filename = "Ener"
              case (5)
                 N = Ndim*Ndim*N_FL*(Ltrot+1); Filename = "P21_GT0"
              case (6)
                 N = Ndim*Ndim*N_FL*(Ltrot+1); Filename = "P21_G0T"
              case (7)
                 ! Metadata counts per tau, printed with N=1, sign=1.
                 N = Ltrot+1; Filename = "P21_Counts"''')
    start = text.index("        subroutine Obser(GR,")
    before, body = text[:start], text[start:]
    body = once(body, "Do I = 1,Size(Obs_scal,1)", "Do I = 1,4  ! P21 vectors are counted by ObserT, not equal-time Obser")
    text = before + body
    start = text.index("        Subroutine ObserT(")
    before, body = text[:start], text[start:]
    body = once(body, "Integer :: IMJ, I, J, I1, J1, no_I, no_J", "Integer :: IMJ, I, J, I1, J1, no_I, no_J, P21_nf, P21_idx")
    body = once(body, "! Standard two-point correlations", """! Native full-matrix buffers; no recomputation and no channel averaging.
           ! Flattening: tau-major, flavor, row i, column j, all zero-based externally.
           If (size(Obs_scal) /= 7) error stop 'P21 observables not allocated'
           If (NT<0 .or. NT>Ltrot) error stop 'P21 tau out of range'
           If (NT == 0) then
              Do I=5,6
                 Obs_scal(I)%N = Obs_scal(I)%N+1
                 Obs_scal(I)%Ave_sign = Obs_scal(I)%Ave_sign+real(ZS,kind(0.d0))
              enddo
           endif
           ! Counts are raw counts, not sign weighted; their print divisor is one.
           Obs_scal(7)%N=1
           Obs_scal(7)%Ave_sign=1.d0
           Obs_scal(7)%Obs_vec(NT+1)=Obs_scal(7)%Obs_vec(NT+1)+1.d0
           Do P21_nf=1,N_FL
              Do I=1,Ndim
                 Do J=1,Ndim
                    P21_idx=1+(((NT*N_FL+P21_nf-1)*Ndim+I-1)*Ndim+J-1)
                    Obs_scal(5)%Obs_vec(P21_idx)=Obs_scal(5)%Obs_vec(P21_idx)+GT0(I,J,P21_nf)*ZP*ZS
                    Obs_scal(6)%Obs_vec(P21_idx)=Obs_scal(6)%Obs_vec(P21_idx)+G0T(I,J,P21_nf)*ZP*ZS
                 enddo
              enddo
           enddo

           ! Standard two-point correlations""")
    return before + body

def alf_main(text):
    text = once(text, "        Use Tau_m_mod", "        Use Tau_m_mod\n        Use Observables, only: Print_bin_Vec")
    text = once(text, "        Integer :: Seed_in", "        Integer :: Seed_in\n        Character(len=16) :: P21_mode\n        Logical :: P21_proof")
    text = once(text, "        Call ham%Alloc_obs(Ltau)", """        Call ham%Alloc_obs(Ltau)
        ! P21 fixed-configuration proof only; ordinary runs take the unchanged path.
        P21_mode=''
        Call get_environment_variable('P21_PROOF_ONLY',P21_mode)
        P21_proof=trim(P21_mode)=='1'
        If (P21_proof) then
           If (Projector .or. Symm .or. reconstruction_needed .or. Ltau/=1) &
                & error stop 'P21 fixed-configuration mode mismatch'
           If (Ndim/=4 .or. N_FL/=2 .or. N_SUN/=1 .or. Ltrot/=10) &
                & error stop 'P21 proof dimensions mismatch'
           ! Deterministic, nonuniform fields; overwrite before native UDV construction.
           Do NT=1,Ltrot
              Do N=1,size(nsigma%f,1)
                 nsigma%f(N,NT)=cmplx(real(mod(3*N+2*NT,7)-3,kind(0.d0))/5.d0,0.d0,kind(0.d0))
              enddo
           enddo
        endif""")
    text = once(text, "        DO  NBC = 1, NBIN", """        If (P21_proof) then
           Call P21_export_propagators()
           Call ham%Init_obs(Ltau)
           ! This is native Tau_m on one fixed configuration, before any MC sweeps.
           Call Tau_m(udvst,GR,PHASE,NSTM,NWRAP,STAB_NT,LOBS_ST,LOBS_EN)
           Do I=5,7
              Call Print_bin_Vec(Obs_scal(I),Group_Comm)
           enddo
           open(unit=5,file='RUNNING',status='old')
           close(5,status='delete')
           write(*,*) 'P21_FIXED_CONFIGURATION_EXPORTED_NOT_SCIENTIFIC_BENCHMARK'
           stop
        endif

        DO  NBC = 1, NBIN""")
    text = once(text, "      end Program Main", """      contains

        Subroutine P21_export_propagators()
          ! Additive proof export. Native PROPR/PROPRM1 and Tau_m remain unmodified.
          Implicit none
          complex(kind=kind(0.d0)) :: P21_B(Ndim,Ndim,N_FL),P21_inv(Ndim,Ndim,N_FL)
          complex(kind=kind(0.d0)) :: P21_K(Ndim,Ndim,N_FL)
          integer :: P21_u,P21_i,P21_j,P21_f,P21_t,P21_n
          open(newunit=P21_u,file='P21_native_propagators.tsv',status='new')
          write(P21_u,'(A)') '# kind tau flavor row col real imag; indices zero-based'
          Do P21_t=1,Ltrot
             P21_B=cmplx(0.d0,0.d0,kind(0.d0))
             P21_inv=P21_B
             P21_K=P21_B
             Do P21_f=1,N_FL
                Do P21_i=1,Ndim
                   P21_B(P21_i,P21_i,P21_f)=1.d0
                   P21_inv(P21_i,P21_i,P21_f)=1.d0
                   P21_K(P21_i,P21_i,P21_f)=1.d0
                enddo
                Call Hop_mod_mmthr(P21_K(:,:,P21_f),P21_f,P21_t)
             enddo
             Call PROPR(P21_B,P21_t)
             Call PROPRM1(P21_inv,P21_t)
             Do P21_f=1,N_FL
                Do P21_i=1,Ndim
                   Do P21_j=1,Ndim
                      write(P21_u,'(A,4(1X,I0),2(1X,ES26.17E3))') &
                           & 'B',P21_t,P21_f-1,P21_i-1,P21_j-1,P21_B(P21_i,P21_j,P21_f)
                      write(P21_u,'(A,4(1X,I0),2(1X,ES26.17E3))') &
                           & 'Binv',P21_t,P21_f-1,P21_i-1,P21_j-1,P21_inv(P21_i,P21_j,P21_f)
                      write(P21_u,'(A,4(1X,I0),2(1X,ES26.17E3))') &
                           & 'Kexp',P21_t,P21_f-1,P21_i-1,P21_j-1,P21_K(P21_i,P21_j,P21_f)
                   enddo
                enddo
                Do P21_n=1,size(Op_V,1)
                   If (Op_V(P21_n,P21_f)%N/=1 .or. Op_V(P21_n,P21_f)%type/=3) &
                        & error stop 'P21 proof requires onsite continuous fields'
                   write(P21_u,'(A,4(1X,I0),2(1X,ES26.17E3))') &
                        & 'field',P21_t,P21_f-1,Op_V(P21_n,P21_f)%P(1)-1,0,nsigma%f(P21_n,P21_t)
                   write(P21_u,'(A,4(1X,I0),2(1X,ES26.17E3))') &
                        & 'coupling',P21_t,P21_f-1,Op_V(P21_n,P21_f)%P(1)-1,0,Op_V(P21_n,P21_f)%g
                enddo
             enddo
          enddo
          Do P21_f=1,N_FL
             Do P21_i=1,Ndim
                Do P21_j=1,Ndim
                   write(P21_u,'(A,4(1X,I0),2(1X,ES26.17E3))') &
                        & 'G00',0,P21_f-1,P21_i-1,P21_j-1,GR(P21_i,P21_j,P21_f)
                enddo
             enddo
          enddo
          close(P21_u)
        end Subroutine P21_export_propagators

      end Program Main""")
    return text

def user_kernel(text):
    text = once(text, "    s, μ, α, Δτ, δG, rng\n) where {T<:Number}", "    s, μ, α, Δτ, δG, rng; observer = nothing\n) where {T<:Number}")
    start = text.index("function local_update_replica!(")
    end = text.index("function local_update_replica_fused!(", start)
    before, body, after = text[:start], text[start:end], text[end:]
    payload = """(; tau=P21_TAU, stage=P21_STAGE, forward=fermion_greens_calc_up.forward,
            phase=sgndetGup*sgndetGdn, Gup_00,Gup_τ0,Gup_0τ,Gup_ττ=Gup,
            Gdn_00,Gdn_τ0,Gdn_0τ,Gdn_ττ=Gdn, Bup,Bdn,s)"""
    callback = lambda tau, stage: "        if observer !== nothing && fermion_greens_calc_up.forward\n            observer(" + payload.replace("P21_TAU",tau).replace("P21_STAGE",'"'+stage+'"') + ")\n        end\n"
    body = once(body, "    for l in fermion_greens_calc_up", callback("0","initial") + "\n    for l in fermion_greens_calc_up")
    body = once(body, "        acceptance_rate += acceptance_rate_l / Lτ", callback("l","pre_stabilization") + "\n        acceptance_rate += acceptance_rate_l / Lτ")
    body = once(body, "        δG = max(δG, δGup, δGdn)", callback("l","post_stabilization") + "\n        δG = max(δG, δGup, δGdn)")
    text = before + body + after
    text = once(text, "function update_replica!(replica, δG)", "function update_replica!(replica, δG; observer = nothing)")
    start=text.index("function update_replica!(replica,")
    end=text.index("function update_replica_fused!",start)
    before,body,after=text[:start],text[start:end],text[end:]
    body=once(body, "        δG, replica.rng\n", "        δG, replica.rng; observer = observer\n")
    return before + body + after

def digest(b):
    return hashlib.sha256(b).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--alf-root',type=Path,required=True)
    p.add_argument('--user-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists(): raise ValueError('Output must be new')
    refs=json.loads((HERE/'upstream_files.json').read_text())
    expected={r['path']:r['sha256'] for r in refs['files']}
    operations=[('alf','Prog/main.F90',alf_main),('alf','Prog/Hamiltonians/Hamiltonian_Hubbard_smod.F90',alf_hubbard),('user','src/dqmc_kernel.jl',user_kernel)]
    prepared=[]
    for origin,rel,fn in operations:
        source=(a.alf_root if origin=='alf' else a.user_root)/rel
        raw=source.read_bytes()
        want=expected[rel] if origin=='alf' else USER_SHA
        if digest(raw)!=want: raise ValueError('Frozen source mismatch: '+str(source))
        patched=fn(raw.decode()).encode()
        prepared.append((origin,rel,raw,patched))
    a.output.mkdir(parents=True)
    rows=[]
    for origin,rel,raw,patched in prepared:
        target=a.output/origin/rel;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(patched)
        patch=''.join(difflib.unified_diff(raw.decode().splitlines(True),patched.decode().splitlines(True),fromfile='a/'+rel,tofile='b/'+rel))
        patchfile=a.output/(origin+'_'+Path(rel).name+'.patch');patchfile.write_text(patch)
        rows.append({'origin':origin,'path':rel,'original_sha256':digest(raw),'observer_sha256':digest(patched),'patch':patchfile.name,'patch_sha256':digest(patch.encode())})
    receipt={'schema':'p21_observer_overlay_v1','status':'GENERATED_NOT_COMPILED_NOT_SCIENTIFICALLY_RUN','alf_commit':PIN,'generator_sha256':digest(Path(__file__).read_bytes()),'files':rows,'sampling_kernel_changes':False,'user_original_preserved':str(a.user_root.resolve()),'alf_original_preserved':str(a.alf_root.resolve())}
    (a.output/'overlay_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()

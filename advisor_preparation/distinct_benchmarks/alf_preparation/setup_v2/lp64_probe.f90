! Setup/link compatibility probe only; no ALF or physical simulation.
program lp64_probe
  use iso_fortran_env, only: int32, real64, error_unit
  implicit none
  integer(int32) :: n(2), lda(2), lwork(2), inc(2), info(2)
  integer(int32), parameter :: sentinel = 123456789_int32
  real(real64) :: x(2), y(2), ra(2,2), rb(2,2), rc(2,2), w(2), rwork(16)
  complex(real64) :: za(2,2), zb(2,2), zc(2,2), hz(2,2), original(2,2), work(100)
  real(real64), external :: ddot
  external :: dgemm, zgemm, zheev
  if (storage_size(0) /= 32 .or. storage_size(0_int32) /= 32) error stop 'Fortran default integer is not32bits'
  ! Padded, zero-high-word inputs are safe even if an unexpected ILP64 library
  ! reads8bytes. The INFO canary detects an8-byte LAPACK output write.
  n=[2_int32,0_int32]; lda=[2_int32,0_int32]
  lwork=[100_int32,0_int32]; inc=[1_int32,0_int32]
  x=[3.0_real64,1415.0_real64]; y=[42.0_real64,0.0_real64]
  if (abs(ddot(n(1),x,inc(1),y,inc(1))-126.0_real64)>1e-12_real64) error stop 'DDOT failed'
  ra=reshape([1.0_real64,2.0_real64,3.0_real64,4.0_real64],[2,2]); rb=0; rb(1,1)=1; rb(2,2)=1; rc=0
  call dgemm('N','N',n(1),n(1),n(1),1.0_real64,ra,lda(1),rb,lda(1),0.0_real64,rc,lda(1))
  if (maxval(abs(rc-ra))>1e-12_real64) error stop 'DGEMM failed'
  za=cmplx(ra,2.0_real64*ra,kind=real64);zb=cmplx(rb,0.0_real64,kind=real64);zc=0
  call zgemm('N','N',n(1),n(1),n(1),cmplx(1.0_real64,0.0_real64,kind=real64), &
       za,lda(1),zb,lda(1),cmplx(0.0_real64,0.0_real64,kind=real64),zc,lda(1))
  if (maxval(abs(zc-za))>1e-12_real64) error stop 'ZGEMM failed'
  original=0;original(1,2)=cmplx(0.0_real64,-1.0_real64,kind=real64);original(2,1)=conjg(original(1,2));hz=original
  info=[-777_int32,sentinel]
  call zheev('V','U',n(1),hz,lda(1),w,work,lwork(1),rwork,info(1))
  if (info(1)/=0) error stop 'ZHEEV INFO failed'
  if (info(2)/=sentinel) error stop 'ILP64 or invalid integer output width detected'
  if (maxval(abs(w-[-1.0_real64,1.0_real64]))>1e-12_real64) error stop 'ZHEEV eigenvalues failed'
  if (maxval(abs(matmul(original,hz)-hz*spread(w,1,2)))>1e-12_real64) error stop 'ZHEEV residual failed'
  if (maxval(abs(matmul(conjg(transpose(hz)),hz)-cmplx(rb,0.0_real64,kind=real64)))>1e-12_real64) error stop 'ZHEEV orthogonality failed'
  print '(A)', 'LP64_BLAS_LAPACK_PROBE_PASS'
  print '(A,I0)', 'DEFAULT_INTEGER_BITS=', storage_size(0)
  print '(A,ES24.16)', 'EIGENPAIR_RESIDUAL=',maxval(abs(matmul(original,hz)-hz*spread(w,1,2)))
end program lp64_probe

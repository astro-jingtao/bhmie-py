! ===========================================================================
! bhmiepy Fortran sources exposed to Python through f2py (module
! _bhmiepy_ext). This file is new code written for bhmiepy (MIT license,
! see the LICENSE file at the root of the repository) -- it is not derived
! from upstream bhmie code.
!
! Design constraints (keep them when editing):
!  - bare subroutines only: no Fortran modules, no "use", no interface
!    blocks, no named kind parameters. This file is the ONLY file parsed
!    by f2py, and those constructs are what f2py's parser handles badly
!    for externally-linked code. External symbols (bhmie_core, bhmie)
!    are resolved at link time from bhmie.f90 / bhmie_f77.f.
!  - plain "double precision" / "double complex" declarations, matching
!    the dp kind used by _fortran/types.f90.
!  - pure ASCII only (f2py reads sources as ASCII).
! ===========================================================================

subroutine backend_info(prec_digits, range_digits)

  ! Report the precision of the real kind used by the extension, so the
  ! Python side can verify the toolchain built double-precision code.

  implicit none

  integer, intent(out) :: prec_digits, range_digits

  prec_digits = precision(1.0d0)
  range_digits = range(1.0d0)

end subroutine backend_info

subroutine bhmie_vec(n, x, m, nang, qext, qsca, qback, g, s1, s2)

  ! Vectorized Mie computation over n (x, m) pairs. Thin loop around the
  ! scalar bhmie_core (which calls the upstream bhmie routine). All input
  ! validation (nang >= 2, series-order limit, Im(m) >= 0) happens on the
  ! Python side BEFORE this is called: the upstream routine responds to
  ! bad input with Fortran "stop", which would kill the whole Python
  ! process.

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: nang
  double precision, intent(in) :: x(n)
  double complex, intent(in) :: m(n)
  double precision, intent(out) :: qext(n), qsca(n), qback(n), g(n)
  double complex, intent(out) :: s1(2*nang-1, n), s2(2*nang-1, n)

  integer :: i

  do i = 1, n
     call bhmie_core(x(i), m(i), nang, s1(:, i), s2(:, i), &
          & qext(i), qsca(i), qback(i), g(i))
  end do

end subroutine bhmie_vec

subroutine bhmie_vec_ang(n, x, m, nang, angles, qext, qsca, qback, g, s1, s2)

  ! Vectorized Mie computation over n (x, m) pairs on a CUSTOM angle grid:
  ! angles(1:nang) must cover 0 to pi/2 (angles(1) == 0 and angles(nang)
  ! == pi/2 exactly, enforced by the upstream routine), and the amplitudes
  ! are returned on the mirrored 0..pi grid (2*nang-1 points). Same
  ! Python-side-validation rule as bhmie_vec applies.

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: nang
  double precision, intent(in) :: x(n)
  double complex, intent(in) :: m(n)
  double precision, intent(in) :: angles(nang)
  double precision, intent(out) :: qext(n), qsca(n), qback(n), g(n)
  double complex, intent(out) :: s1(2*nang-1, n), s2(2*nang-1, n)

  integer :: i

  do i = 1, n
     call bhmie_core_ang(x(i), m(i), nang, angles, s1(:, i), s2(:, i), &
          & qext(i), qsca(i), qback(i), g(i))
  end do

end subroutine bhmie_vec_ang

subroutine bhmie_dust_accum(n, nwav, nang, x, m, weights, cross, angles, &
     & cext, csca, cback, gsca, s11, s12, s33, s34)

  ! One call of the dust pipeline's inner accumulation (port of the loop
  ! body of upstream bhmie_caller.f90 compute_dust_properties, driven
  ! through the package's optimized core): for each of the n size bins,
  ! run the core on every wavelength of the x(n, nwav) / m(n, nwav)
  ! grids and accumulate the weighted running totals IN FORTRAN. This
  ! avoids materializing the huge per-point s1/s2 amplitude arrays on
  ! the Python side, which dominated the pipeline's non-core runtime.
  !
  ! weights(n) hold each bin's number weight (already including the
  ! component's number abundance); cross(n) the geometric cross-sections
  ! in cm^2. Cross-section totals use weights*cross, the scattering-
  ! matrix totals use weights alone -- exactly the upstream convention.
  ! The accumulators are intent(inout) so consecutive components (and
  ! bins with weight <= 0, which the Python side excludes beforehand,
  ! preserving the upstream 'if (weight_number > 0)' skip) add into the
  ! same arrays. Bins accumulate sequentially, like upstream's own loop.
  ! Same Python-side-validation rule as bhmie_vec applies.

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: nwav
  integer, intent(in) :: nang
  double precision, intent(in) :: x(n, nwav)
  double complex, intent(in) :: m(n, nwav)
  double precision, intent(in) :: weights(n)
  double precision, intent(in) :: cross(n)
  double precision, intent(in) :: angles(nang)
  double precision, intent(inout) :: cext(nwav), csca(nwav), cback(nwav), gsca(nwav)
  double precision, intent(inout) :: s11(nwav, 2*nang-1), s12(nwav, 2*nang-1), &
       & s33(nwav, 2*nang-1), s34(nwav, 2*nang-1)

  integer :: i, iw, j, nang2
  double precision :: qext, qsca, qback, g, w, wc, a1sqj, a2sqj
  double complex, allocatable :: s1(:), s2(:)

  nang2 = 2*nang - 1
  allocate(s1(nang2), s2(nang2))

  do i = 1, n
     w = weights(i)
     wc = w * cross(i)
     do iw = 1, nwav
        call bhmie_core_ang(x(i, iw), m(i, iw), nang, angles, s1, s2, &
             & qext, qsca, qback, g)
        cext(iw) = cext(iw) + qext * wc
        csca(iw) = csca(iw) + qsca * wc
        cback(iw) = cback(iw) + qback * wc
        gsca(iw) = gsca(iw) + g * qsca * wc
        do j = 1, nang2
           a1sqj = abs(s1(j)) * abs(s1(j))
           a2sqj = abs(s2(j)) * abs(s2(j))
           s11(iw, j) = s11(iw, j) + w * 0.5d0 * (a1sqj + a2sqj)
           s12(iw, j) = s12(iw, j) + w * 0.5d0 * (-a1sqj + a2sqj)
           s33(iw, j) = s33(iw, j) + w * dble(s2(j) * dconjg(s1(j)))
           s34(iw, j) = s34(iw, j) + w * aimag(s2(j) * dconjg(s1(j)))
        end do
     end do
  end do

  deallocate(s1, s2)

end subroutine bhmie_dust_accum

subroutine bhmie_vec_upstream(n, x, m, nang, qext, qsca, qback, g, s1, s2)

  ! Vectorized Mie computation over n (x, m) pairs through the PRISTINE
  ! upstream routine (bhmie_upstream.f90: verbatim body, including the
  ! original allocate(d(nmxx)) per call). Kept for benchmarking and
  ! auditing the package's optimized copy; not used by the public API.
  ! Same Python-side-validation rule as bhmie_vec applies.

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: nang
  double precision, intent(in) :: x(n)
  double complex, intent(in) :: m(n)
  double precision, intent(out) :: qext(n), qsca(n), qback(n), g(n)
  double complex, intent(out) :: s1(2*nang-1, n), s2(2*nang-1, n)

  integer :: i

  do i = 1, n
     call bhmie_upstream_core(x(i), m(i), nang, s1(:, i), s2(:, i), &
          & qext(i), qsca(i), qback(i), g(i))
  end do

end subroutine bhmie_vec_upstream

subroutine bhmie_f77_ref(x, refrel, nang, qext, qsca, qback, gsca, s1, s2, ierr)

  ! Scalar reference computation through the ORIGINAL fixed-form F77
  ! BHMIE (bhmie_f77.f, Bohren & Huffman 1983 App. A / Draine). Its
  ! Fortran interface is single precision, so values are converted on the
  ! way in and out. Used by the test suite as an independent oracle for
  ! the double-precision port; kept in the extension so it can also be
  ! used interactively to cross-check results.
  !
  ! ierr = 0 on success; 1 if nang is outside [2, 1000] (MXNANG); 2 if the
  ! Mie series order exceeds the F77 NMXX = 150000. The F77 routine would
  ! otherwise STOP and kill the whole Python process, so callers MUST
  ! check ierr.

  implicit none

  integer, intent(in) :: nang
  double precision, intent(in) :: x
  double complex, intent(in) :: refrel
  double precision, intent(out) :: qext, qsca, qback, gsca
  double complex, intent(out) :: s1(2*nang-1), s2(2*nang-1)
  integer, intent(out) :: ierr

  real :: xr, qextr, qscar, qbackr, gscar
  complex :: refrelr, s1r(2*nang-1), s2r(2*nang-1)
  double precision :: xstop, nmx_d
  integer :: j, nmx

  ierr = 0

  if(nang < 2 .or. nang > 1000) then
     ierr = 1
     qext = 0.d0
     qsca = 0.d0
     qback = 0.d0
     gsca = 0.d0
     s1 = (0.d0, 0.d0)
     s2 = (0.d0, 0.d0)
     return
  end if

  ! mirror of the F77 series-order check (single-precision routine,
  ! NMXX = 150000 there)
  xstop = x + 4.d0*x**0.3333d0 + 2.d0
  nmx_d = max(xstop, abs(refrel)*x) + 15.d0
  nmx = nint(nmx_d)
  if(nmx > 150000) then
     ierr = 2
     qext = 0.d0
     qsca = 0.d0
     qback = 0.d0
     gsca = 0.d0
     s1 = (0.d0, 0.d0)
     s2 = (0.d0, 0.d0)
     return
  end if

  xr = real(x)
  refrelr = cmplx(refrel)
  call bhmie(xr, refrelr, nang, s1r, s2r, qextr, qscar, qbackr, gscar)

  qext = dble(qextr)
  qsca = dble(qscar)
  qback = dble(qbackr)
  gsca = dble(gscar)
  do j = 1, 2*nang-1
     s1(j) = dcmplx(s1r(j))
     s2(j) = dcmplx(s2r(j))
  end do

end subroutine bhmie_f77_ref

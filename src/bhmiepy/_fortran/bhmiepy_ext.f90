! ===========================================================================
! bhmiepy Fortran sources exposed to Python through f2py (module
! _bhmiepy_ext). This file is new code written for bhmiepy (MIT license,
! see the LICENSE file at the root of the repository) -- it is not derived
! from upstream bhmie code.
! ===========================================================================

subroutine backend_info(prec_digits, range_digits)

  ! Report the precision of the real kind used by the extension, so the
  ! Python side can verify the toolchain built double-precision code.

  implicit none

  integer, intent(out) :: prec_digits, range_digits

  prec_digits = precision(1.0d0)
  range_digits = range(1.0d0)

end subroutine backend_info

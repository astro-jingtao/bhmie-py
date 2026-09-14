! ===========================================================================
! This file is part of bhmiepy. It is a copy of a source file from the bhmie
! project (https://github.com/hyperion-rt/bhmie), snapshot taken at upstream
! commit 37c072909a67d1027120762680a1c0b350875398 (2013-10-08), and may be
! modified for bhmiepy. A pristine copy of the upstream tree is kept in the
! top-level bhmie/ directory of the bhmiepy repository for reference and
! testing.
!
! Modifications for bhmiepy (relative to the upstream file): none so far.
!
! ---------------------------------------------------------------------------
! Original copyright notice and license (BSD 2-Clause):
!
! Copyright (c) 2012-13, Thomas P. Robitaille
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! The underlying BHMIE subroutine originates from Appendix A of
! Bohren & Huffman, "Absorption and Scattering of Light by Small Particles"
! (Wiley, 1983), with extensive modifications by B. T. Draine (see the
! history log in upstream/original/bhmie.f in the upstream submodule).
! ===========================================================================

module types

  implicit none
  save

  integer,parameter :: sp = selected_real_kind(p=6,r=37)
  integer,parameter :: dp = selected_real_kind(p=15,r=307)

end module types

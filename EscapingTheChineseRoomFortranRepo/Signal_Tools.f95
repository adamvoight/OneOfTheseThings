!======================================================================
MODULE Signal_Tools

use Declare_Constants
use number_generators

implicit none

private

public   :: ConvolveDFT, CorrelateDFT, DeScramble

!======================================================================
                   CONTAINS
!======================================================================


!****************************************************************
! This function returns the circular convolution of the argument
! vectors V1 and V2 with the same dimensionality, but uses DFT
!----------------------------------------------------------------
function ConvolveDFT(V1, V2)
  integer :: i, j, n
  real*8, dimension(0:D-1) :: ConvolveDFT, sum, V1, V2
  complex*16 :: FV1(D/2), FV2(D/2)
  

  real*8 W(20000)
  n = D
  call DFFTI(n,W)
  call DFFTF(n,V1,W)
  call DFFTF(n,V2,W)
  FV1 = transfer(V1,FV1)
  FV2 = transfer(V2,FV2)
  FV1 = FV1 * FV2
  V1 = transfer(FV1,V1)
  call DFFTB(n,V1,W)
     
  ConvolveDFT = V1
 
end function ConvolveDFT
!****************************************************************


!****************************************************************
! This function returns the circular correlation of the argument
! vectors V1 and V2 with the same dimensionality using DFT
!----------------------------------------------------------------
function CorrelateDFT(V1, V2)
  integer :: i, j, n
  real*8, dimension(0:D-1) :: CorrelateDFT, sum, V1, V2
  complex*16 :: FV1(D/2), FV2(D/2)

  real*8 W(20000)
  n = D
  call DFFTI(n,W)
  call DFFTF(n,V1,W)
  call DFFTF(n,V2,W)
  FV1 = transfer(V1,FV1)
  FV2 = transfer(V2,FV2)
  FV1 = CONJG(FV1) * FV2
  V1 = transfer(FV1,V1)
  call DFFTB(n,V1,W)
     
  CorrelateDFT = V1
 
end function CorrelateDFT
!****************************************************************

!****************************************************************
function DeScramble(Vector, Scram_Indices)
  real*8, dimension(D) :: DeScramble, Vector, Temp
  integer :: Scram_Indices(D), i
 do i = 1, D
   Temp(Scram_Indices(i)) = Vector(i)
 enddo
 DeScramble = Temp

end function DeScramble
!****************************************************************



END MODULE Signal_Tools


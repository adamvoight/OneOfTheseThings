
module number_generators
 implicit none
 integer, private             :: idum, inext, inextp, inext1, inextp1
 real*8, private, dimension(55) :: ma, ma2
 logical                      :: switch

 CONTAINS

!***********************************************************************
!  ran3:  generates random values 0..1]
!  seed (idum) is set to a negative integer on entry
!  From Press et al.,(1992) Numerical recipes in FORTRAN (2nd ed.) CUP
!
!  This routine runs through the ma vector pre-established when the
!       generator is initialized
!***********************************************************************
function ran3(idum)
 implicit none
 real*8                  :: ran3, mj, mk
 real*8, parameter       :: MBIG =4000000.0, MSEED =1618033.0, MZ =0.0, FAC = 1.0/MBIG
 integer               :: i, ii, k, idum

 inext  = inext + 1
 if(inext == 56) inext = 1
 inextp = inextp + 1
 if(inextp == 56)inextp= 1

 mj = ma(inext) - ma(inextp)

 if(mj < MZ) mj = mj + MBIG
 ma(inext) = mj
 ran3 = mj*FAC
 return
end function ran3



!*********************************************************************
! rnorm:  Unit Normal distribution
!*********************************************************************
function rnorm(idum)
 implicit none
 real*8                    :: fac, r, v1, v2, rnorm
 real*8,    save           :: rnorm2
 integer                 :: idum

 if (switch) then
10  v1 = 2.0 * ran3(idum) - 1.0
    v2 = 2.0 * ran3(idum) - 1.0
    r = v1**2 + v2**2
    if ((r .ge. 1.0).or.(r .eq.0)) goto 10

    fac = sqrt(-2.0 * log(r)/r)
    rnorm2 = v1 * fac
    switch = .false.
    rnorm  = v2 * fac
 else
    switch = .true.
    rnorm  = rnorm2
 endif
 return
end function rnorm



!***********************************************************************
! gaussian:  returns a Normal deviate from N(mu, sd)
!***********************************************************************
function gaussian(mu, sd)
 real*8, intent(in)        :: mu, sd
 real*8                    :: gaussian

 gaussian = (rnorm(idum) * sd) + mu   ! Calculate Normal(mu, sigma)
 return
end function gaussian



!**********************************************************************
!  FixSeed:  Sets a fixed seed (-99) for the random-number routine
! and sets up the ma vector for the random number generator
!**********************************************************************
subroutine FixSeed
 real*8                  :: mj, mk
 real*8, parameter       :: MBIG =4000000.0, MSEED =1618033.0, MZ =0.0
 integer               :: i, ii, k
 idum = -99
 write(*,'(A, I3)') '  Seed for random-number generator = ', idum

      mj= MSEED - iabs(idum)
      mj = mod(mj, MBIG)
      ma(55) = mj
      mk=1
      do i=1,54
         ii = mod(21*i, 55)
         ma(ii) = mk
         mk = mj - mk
         if(mk < MZ) mk = mk + MBIG
         mj = ma(ii)
      enddo
      do k = 1, 4
         do i = 1, 55
            ma(i) = ma(i) - ma(1 + mod(i+30, 55) )
            if(ma(i) .lt. MZ) ma(i) = ma(i) + MBIG
         enddo
      enddo
      inext = 0
      inextp = 31
 return
end subroutine FixSeed



!**********************************************************************
!  RandSeed:  Sets a random seed for the random-number routines
!  and sets up the ma vector for the random number generator
!**********************************************************************
subroutine RandSeed
 implicit none
 real*8                  :: mj, mk
 real*8, parameter       :: MBIG =4000000.0, MSEED =1618033.0, MZ =0.0
 integer               :: i, ii, k
 call system_clock(idum)

10 if (idum >  1000000) idum = idum / 10
   if (idum > 1000000) go to 10

   if (idum > 0) idum = idum *(-1)

  ! write(*,'(A, I8)') '  Seed for random-number generator = ', idum

        mj= MSEED - iabs(idum)
        mj = mod(mj, MBIG)
        ma(55) = mj
        mk=1
        do i=1,54
          ii = mod(21*i, 55)
          ma(ii) = mk
          mk = mj - mk
          if(mk < MZ) mk = mk + MBIG
          mj = ma(ii)
        enddo
        do k = 1, 4
          do i = 1, 55
             ma(i) = ma(i) - ma(1 + mod(i+30, 55) )
             if(ma(i) .lt. MZ) ma(i) = ma(i) + MBIG
          enddo
        enddo
        inext = 0
        inextp = 31

   return
end subroutine RandSeed


!**********************************************************************
! get_ran_seed:  gets current state of the random-number generator
!**********************************************************************
subroutine get_ran_seed(dummy, MAA)
  integer, intent(out), dimension(2) :: dummy
  real*8, dimension(55)                :: MAA
    dummy(1) = inext
    dummy(2) = inextp
    MAA = ma
  return
end subroutine get_ran_seed


!**********************************************************************
! assign_seed:  restores the state of the random-number generator
!**********************************************************************
subroutine assign_seed(iseed, MAA)
  integer, dimension(2), intent(in)  :: iseed
  integer                            :: get_ran_seed
  real*8, dimension(55)                :: MAA
    inext  = iseed(1)
    inextp = iseed(2)
    ma     = MAA
    switch = .true.
  return
end subroutine assign_seed

END MODULE Number_Generators

! This is a quick program to load the learned memory matrix
! and find the N Nearest neighbors to a target word in the space
! Compile: f95 Neighbors.f95 -o Neighbors
! Syntax: Neighbors <WORD> <N>
! E.g.:  "Neighbors baseball 20" finds 20 most similar memory 
! vectors to the baseball vector
!
! For example, in context space, "Neighbors baseball 20" from the
! small Wikipedia corpus returns: 
!
! ================================================
! 20 Nearest Neighbors to [batter]:
! ================================================
! batter               1.0000001
! ball                 0.7821321
! hit                  0.76954096
! fair                 0.7437778
! pitch                0.72545194
! ground               0.70555484
! attempt              0.70481705
! base                 0.69794034
! put                  0.69120353
! runners              0.6387966
! foul                 0.6332392
! plate                0.62975984
! touch                0.6270346
! fly                  0.6193897
! strike               0.6137481
! strikes              0.60961693
! contact              0.6070126
! pitcher              0.598485
! swing                0.592292
! caught               0.58490586


!======================================================================
MODULE Declare_Constants

implicit none

 integer,   public, parameter :: D              = 2048,   &
				 MAX_WORDS      = 40000,  &
				 WORD_CHARS     = 20
				
END MODULE Declare_Constants
!======================================================================

program N_Nearest_Neighbors

USE Declare_Constants

implicit none
  
  character(len=WORD_CHARS) :: Word(MAX_WORDS), cos_file, W1, W2
  integer                   :: Words_Learned, i, j, P1(D), P2(D), Pos1, Pos2, N, n_args, iargc
  real*8                      :: Memory(D,MAX_WORDS)
  integer, allocatable      :: High_Resp(:)
  real*8  ,allocatable      :: Resp_Activation(:)			     
  logical                   :: FLAG1, FLAG2, SWITCH

!=======================================
 
 n_args = iargc()
  
 call getarg(1, W1)
 call getarg(2, W2)
 read(W2,*) N 
 write(*,*) '================================================' 
 write(*,*) N, 'Nearest Neighbors to [', trim(W1), ']:'
 write(*,*) '================================================'

 call Read_Matrix(Memory, Words_Learned, Word)   

 Pos1 = 0
 do i = 1, MAX_WORDS
   if (W1 == Word(i)) then 
      Pos1 = i
      exit
   endif
 enddo
  
 allocate(High_Resp(N), Resp_Activation(N))
 call Highest_N_Activated(Memory(:,Pos1), Memory, High_Resp, Resp_Activation, Words_Learned, N)
 
 do i = 1, N
   write(*,*) Word(High_Resp(i)), Resp_Activation(i)
 enddo

 deallocate(High_Resp, Resp_Activation)

!======================================= 
 CONTAINS
!=======================================  


!*************************************************************************************
subroutine Read_Matrix(Memory, Words_Learned, Word)
  character(len=WORD_CHARS) :: Word(MAX_WORDS)
  real*8                  :: Memory(D,MAX_WORDS)
  integer                   :: Words_Learned, i

 
 open(unit=10, file='word_labels.txt', status='old') 

 open(unit=11, file='matrix.mat', status='old', form='unformatted')

 read(10,*) Words_Learned 
 
 do i = 1, Words_Learned
   read(10,'(a20)') Word(i)
   read(11) Memory(:,i)
 enddo

 close(10)
 close(11)

 end subroutine Read_Matrix 
!*************************************************************************************


!****************************************************************
subroutine Highest_N_Activated(Context, Memory, High_Resp, Resp_Activation, Words_Learned, num)
  integer  :: i, j, num, High_Resp(num), Words_Learned, loc
  real*8 :: Memory(D,MAX_WORDS), Context(D), Resp_Activation(num), &
              Cosine_Array(Words_Learned), Sim, Max

do i = 1, Words_Learned
  Cosine_Array(i) = Vector_Cosine(Context, Memory(:,i),D)
enddo

do i = 1, num
  Max = -1.0
  do j = 1, Words_Learned
    if(Cosine_Array(j) > Max) then
      Max = Cosine_Array(j)
      loc = j
    endif
  enddo
  High_Resp(i) = loc
  Resp_Activation(i) = Max
  Cosine_Array(loc) = 0.0
enddo

end subroutine Highest_N_Activated
!****************************************************************


!================================================================
!          PROXIMITY ROUTINES:
!================================================================
!****************************************************************
function Vector_Length (Vector, n)
   integer :: i, n 
   real*8    :: Vector_Length, Vector(n), SS

 SS = 0.0
  SS = dot_product(Vector, Vector)
 Vector_Length = sqrt(SS)
end function Vector_Length 
!****************************************************************



!****************************************************************
subroutine Normalize (Vector_In, n)
   integer             :: i, n
   real*8              :: Vector_In(n), Vect_Length, Test_Val
  
  Test_Val = dot_product(Vector_In, Vector_In)
  if (Test_Val == 0.0) then 
    return
  endif
  
  
  Vect_Length = 0.0 
  Vect_Length = Vector_Length(Vector_In, n)
  Vector_In = Vector_In / Vect_Length
end subroutine Normalize
!****************************************************************


!****************************************************************
function Vector_Cosine (Vector1, Vector2, n)
   
   real*8    :: Vector1(n), Vector2(n)
   real*8    :: Vector_Cosine
   integer :: n, i
 
 call Normalize(Vector1, n) 
 call Normalize(Vector2, n) 
 Vector_Cosine = dot_product(Vector1, Vector2)
 
end function Vector_Cosine
!****************************************************************



end program N_Nearest_Neighbors

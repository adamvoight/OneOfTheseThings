!===================================================================== 
!
! Syntax: Neighbors <WORD> <N>
!
!======================================================================
MODULE Declare_Constants

implicit none

 integer,   public, parameter :: D = 2048,   &
					MAX_WORDS      = 41000,  &
					WORD_CHARS     = 20
				
END MODULE Declare_Constants
!======================================================================

program N_Nearest_Neighbors

USE Declare_Constants

implicit none
  
  character(len=WORD_CHARS) :: Word(MAX_WORDS), Vocab(MAX_WORDS), W1, W2
  integer                   :: i, j, Pos1, N, Words_Learned, KnownWords
  real*8                    :: Memory(D,MAX_WORDS)
  integer, allocatable      :: High_Resp(:)
  real*8  ,allocatable      :: Resp_Activation(:)

!=======================================
  
 call getarg(1, W1)
 call getarg(2, W2)
 read(W2,*) N
 write(*,*) '=======================================================' 
 write(*,*) '  Nearest Well-Known Neighbors to [', trim(W1), ']:'
 write(*,*) '======================================================='

 call Read_Matrix(Memory, Words_Learned, Word, Vocab, KnownWords)   

 Pos1 = 0
 do i = 1, MAX_WORDS
   if (W1 == Word(i)) then 
      Pos1 = i
      exit
   endif
 enddo
  
 allocate(High_Resp(KnownWords), Resp_Activation(KnownWords))
 call Highest_N_Activated(Memory(:,Pos1), Memory, High_Resp, Resp_Activation, Words_Learned, KnownWords)
 
 j=1
 do i = 1, KnownWords
   if(isWellKnown(Word(i), Vocab, KnownWords) .AND. Resp_Activation(i) .GT. 0.1) then
   	write(*,*) Word(High_Resp(i)), Resp_Activation(i)
   	j=j+1
   endif
   if(j==N) then
   	exit
   endif
 enddo

 deallocate(High_Resp, Resp_Activation)

!======================================= 
 CONTAINS
!=======================================  

function isWellKnown(wordIn, Vocab, KnownWords)
 character(len=WORD_CHARS) :: Vocab(MAX_WORDS), wordIn
 integer                   :: i, j, KnownWords
 logical                   :: isWellKnown
 
 isWellKnown = .False.
 do i=0, KnownWords
 	if ( Vocab(i) == wordIn ) then
 		isWellKnown = .True.
 		exit
 	endif
 enddo
 
 end function isWellKnown

!*************************************************************************************
subroutine Read_Matrix(Memory, Words_Learned, Word, Vocab, KnownWords)
  character(len=WORD_CHARS) :: Word(MAX_WORDS), Vocab(MAX_WORDS)
  real*8                    :: Memory(D,MAX_WORDS)
  integer                   :: Words_Learned, KnownWords, i

 
 open(unit=10, file='word_labels.txt', status='old') 
 open(unit=11, file='well_known.txt', status='old')
 open(unit=12, file='matrix.mat', status='old', form='unformatted')

 read(10,*) Words_Learned 
 read(11,*) KnownWords
 
  do i=1, KnownWords
 	read(11,'(a20)') Vocab(i)
 enddo
 
 close(11)
 
 do i = 1, Words_Learned
   read(10,'(a20)') Word(i)
   read(12) Memory(:,i)
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

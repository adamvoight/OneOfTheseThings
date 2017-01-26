!====================================================================================
!
! This program returns the word with least semantic similarity
! Syntax: Neighbors <WORD> <WORD> <WORD> <WORD>
!
!====================================================================================

MODULE Declare_Constants

implicit none

 integer,   public, parameter :: D              = 2048,   &
				 MAX_WORDS      = 40000,  &
				 WORD_CHARS     = 20
				
END MODULE Declare_Constants


!************************************************************************************

program OneOfTheseThings

USE Declare_Constants

implicit none
  
  character(len=WORD_CHARS) :: Word(MAX_WORDS), args(5)
  integer                   :: Words_Learned, probnum, ind(4)
  real*8                    :: Memory(D,MAX_WORDS)	     

!take in arguments to args
 call getarg(1, args(1))
 call getarg(2, args(2))
 call getarg(3, args(3))
 call getarg(4, args(4))
 call getarg(5, args(5))

 read(args(5), *) probnum


!read matrix into memory
 call Read_Matrix(Memory, Words_Learned, Word)   

!get location of arguments in matrix
 call Get_Indicies(args, ind, Word)

!write least similar word to screen and append to file
 write(*,*) "Least similar word: ", args(Least_Similar(ind, Memory))
 open(unit=40,file='results.txt',Access = 'append',Status='old')
 write(40,*) probnum, Least_Similar(ind, Memory)

 contains

!************************************************************************************

subroutine Get_Indicies(args, ind, Word)
  character(len=WORD_CHARS) :: args(4), Word(MAX_WORDS)
  integer                   :: ind(4), i

  do i = 1, MAX_WORDS
     if(args(1)==Word(i))then
	ind(1)=i
     endif
     if(args(2)==Word(i))then
	ind(2)=i
     endif
     if(args(3)==Word(i))then
	ind(3)=i
     endif
     if(args(4)==Word(i))then
	ind(4)=i
     endif
  enddo

end subroutine Get_Indicies

!*************************************************************************************

subroutine Read_Matrix(Memory, Words_Learned, Word)
  character(len=WORD_CHARS) :: Word(MAX_WORDS)
  real*8                    :: Memory(D,MAX_WORDS)
  integer                   :: Words_Learned, i

 
 open(unit=10, file='word_labels.txt', status='old') 

 open(unit=11, file='matrix.mat', status='old', form='unformatted')

!read first line of word_labels (number of unique words) into Words_Learned
 read(10,*) Words_Learned 
 
!read matrix and word labels into memory
 do i = 1, Words_Learned
!read next 20 chars of word_labels into word array at current index
   read(10,'(a20)') Word(i)
!read next entry of matrix into memory array
   read(11) Memory(:,i)
 enddo

!close files
 close(10)
 close(11)

 end subroutine Read_Matrix 

!*************************************************************************************

function Least_Similar(ind, Memory)
  real*8    :: omission(4), Memory(D,MAX_WORDS), maxo
  integer   :: ind(4), i, Least_Similar

  omission(1) = Similarity(Memory(:,ind(2)), Memory(:,ind(3)), Memory(:,ind(4)))
  omission(2) = Similarity(Memory(:,ind(1)), Memory(:,ind(3)), Memory(:,ind(4)))
  omission(3) = Similarity(Memory(:,ind(1)), Memory(:,ind(2)), Memory(:,ind(4)))
  omission(4) = Similarity(Memory(:,ind(1)), Memory(:,ind(2)), Memory(:,ind(3)))

  maxo = -1.0
  Least_Similar = 0
  do i = 1, 4
     if(omission(i)>maxo)then
	maxo=omission(i)
	Least_Similar=i
     endif
  enddo

end function Least_Similar

!************************************************************************************

function Similarity(Vector1, Vector2, Vector3)
   real*8    :: Similarity, Vector1(D), Vector2(D), Vector3(D)

   Similarity = 1.0

   Similarity = Similarity * Vector_Cosine(Vector1, Vector2, D)
   Similarity = Similarity * Vector_Cosine(Vector2, Vector3, D)
   Similarity = Similarity * Vector_Cosine(Vector3, Vector1, D)

end function Similarity

!====================================================================================
!          PROXIMITY ROUTINES:
!====================================================================================

function Vector_Length (Vector, n)
   integer :: i, n 
   real*8  :: Vector_Length, Vector(n), SS

 SS = 0.0
  SS = dot_product(Vector, Vector)
 Vector_Length = sqrt(SS)
end function Vector_Length 

!************************************************************************************

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

!************************************************************************************

function Vector_Cosine (Vector1, Vector2, n)
   
   real*8  :: Vector1(n), Vector2(n)
   real*8  :: Vector_Cosine
   integer :: n, i
 
 call Normalize(Vector1, n) 
 call Normalize(Vector2, n) 
 Vector_Cosine = dot_product(Vector1, Vector2)
 
end function Vector_Cosine

!************************************************************************************

end program OneOfTheseThings

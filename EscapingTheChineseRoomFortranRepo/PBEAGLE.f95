! This code implements the encoding of the BEAGLE model from:
! Jones, M. N., & Mewhort, D. J. K. (2007). Representing word meaning
! and order information in a composite holographic lexicon. Psychological
! Review, 114, 1-37. 
!
! As written, it uses the NAG Fast-Fourier Transform to perform convolution
! in parallel. If you do not have access to NAG, you can substitute the 
! inline routine "circ_conv" on line 274 to perform circular convolution,
! however, this method is not computationally efficient for large text corpora
!
! Output files are: 
!  word_labels.txt: the text label for each vector in the order they are encoded
!  matrix.mat: the binary vectors 1..D representing the context and order information
!              for each word in the order they are encoded
!  visual.mat: the "environment" vectors. These are the vectors representing the 
!              enviornmental characteristics of the word used in encoding. Each 
!              vector is randomly initialized each time the program is run
!ownership.mat: The placeholder phi vector used in order encoding; this vector
!               is randomly initialized each time the program is run, and is 
!               a "key" required to decode positional information at a later time
!scrambles.mat: This file contains the scrambled order of P1 and P2 (see Jones & 
!               Mewhort appendix) corresponding to a particular run of the model. 
!               P1 and P2 are randomly initialized at each program run
!
!
!  This demo uses a small text corpus, "Wiki.txt," borrowed from Wikipedia. To use a
!  custom text corpus, replace "wiki.txt" with your filename in line 170. The file must
!  already be parsed with one sentence per line, and all in lowercase with no formatting
!
!  The Wikipedia corpus only contains 3944 unique terms, so change the MAX_WORDS parameter
!  in the Declare_Constants module to this value. The program will then allocate a 3944-by-
!  [vector dimensionality] matrix to represent the words, rather than the 90K it is presently
!  set to. For larger corpora, set MAX_WORDS to as large a value as your system will permit 
!  (or equal to the number of unique terms in the corpus, if you know that a priori). If there
!  are more unique terms in the corpus than MAX_WORDS, the first MAX_WORDS of them will 
!  gain lexical entries, and any others will not.

!  For more information, please contact Mike Jones: jonesmn@indiana.edu


program PBEAGLE

USE number_generators
USE Declare_Constants
USE Reading_Tools
USE Signal_Tools


implicit none
  character(len=WORD_CHARS) :: Word(MAX_WORDS), Stop_List(STOP_WORDS)
  character(len=WORD_CHARS), &
                allocatable :: Sentence_Word(:) 
  character(len=SENT_CHARS) :: Sentence  
  real*8, allocatable     :: Memory(:,:), Visual(:,:), W(:,:), Context(:,:), Ownership(:), &
  			       Order(:,:), Probe_Vector(:)
  integer                   :: N, endpoint, i, j, Words_Learned, High_Resp(N_NEIGHBOR),&
                               blank, idum, line, blank_pos, P1(D), P2(D), IFAIL
  integer, allocatable      :: loc(:)
  real*8                  :: Resp_Activation(N_NEIGHBOR), Sigma
  logical                   :: CLEAN_WORD, BAD_SENT
  logical, allocatable      :: Ignore(:)
  real                      :: T1, T2, learn_time, Frequency(MAX_WORDS)
  real*8, parameter       :: Mu = 0.0D0

  
 EXTERNAL C06EKF
   
!**********************BEGIN MAIN*********************************************
 call cpu_time(T1)
 allocate(Memory(D,MAX_WORDS), Visual(D,MAX_WORDS), Ownership(D))
 call Init(Memory, Visual, Word, Frequency, Stop_List)
 call Make_Scramble(P1, P2)
 Sigma = sqrt(1.0/D)
 Ownership = Random_Vector(D, Mu, Sigma)
 write(*,*) Ownership(1000)
 Words_Learned = 0
 line = 0
 call cpu_time(T2)
 learn_time = (T2-T1) / 60 
 write(*,*) 'Memory allocated... ', learn_time, ' minutes to allocate memory.'
 call cpu_time(T1)
 write(*,*) 'Penfield BEAGLE - Learning textbase...' 
 
 do !Sentence
   9 BAD_SENT = .false.
   call Read_Sentence(Sentence, 1, BAD_SENT)
   if (BAD_SENT) goto 9
   line = line + 1
   if (E_FILE .or. (line > TIME_OUT)) goto 10
   call Count_Words(Sentence, BAD_SENT, N)
   if (BAD_SENT) goto 9
   allocate (Sentence_Word(N), W(D,N), Order(D,N), Context(D,N), loc(N), Ignore(N))
   W = 0.0
   Order = 0.0
   loc = 0
   call Parse_Sentence(Sentence, Sentence_Word, Stop_List, Ignore, N)
   do i = 1, N
      call Does_Word_Exist(Sentence_Word(i), Word, loc(i), endpoint)
      if (loc(i) /= 0) then
        W(:,i) = Visual(:,loc(i))
      else
        W(:,i) = Random_Vector(D, Mu, Sigma)
	loc(i) = endpoint
        if (Words_Learned == MAX_WORDS) goto 8
	call Create_New_Word(Sentence_Word(i), Word, W(:,i), Visual, loc(i))
	Words_Learned = Words_Learned + 1
    8 endif
    Frequency(loc(i)) = Frequency(loc(i)) + 1.0
   enddo
   
  !Coding Item info: 
  call Compute_Item_Info(Context, W, N, Frequency, loc, Ignore)
  
  !Coding Order info: 
  call Compute_Convolutions(Ownership, Order, W, N)   
  
  !Normalizing Item and Order info to same scale: 
  do i = 1, N
    call Normalize(Context(:,i),D)
    call Normalize(Order(:,i),D)
  enddo
  
  ! Update memory:
  do i = 1, N
    Memory(:,loc(i)) = Memory(:,loc(i)) + Context(:,i)  + Order(:,i)
  ! Try using just the context information just for fun
  !  Memory(:,loc(i)) = Memory(:,loc(i)) + Context(:,i)
  enddo
  

   deallocate(Sentence_Word, W, Order, Context, loc, Ignore)
   write(*,*) 'Learned Line ', line
 enddo !Sentence
 10 write(*,*) 'Processing Complete'
   write(*,*) 'Words Learned: ', Words_Learned
   
   !Normalizing Memory: 
 do i = 1, Words_Learned
   call Normalize(Memory(:,i), D)
 enddo

 close(2)
!-----------------------------------------------------------

 call cpu_time(T2)
 learn_time = (T2-T1) / 60
 write(*,*) 'Time taken to learn: ', learn_time, ' minutes with lexicon of ', Words_Learned, &
            ' words and ', line, ' lines of text.'

 write(*,*) 'Storing output matrix...'
 call Store_Matrix(Memory, Visual, Word, Words_Learned, Ownership, P1, P2)
 
 call Wrap_Up

!**********************END   MAIN*********************************************

!======================================================================
                 CONTAINS
!======================================================================


!****************************************************************
subroutine Init(Memory, Visual, Word, Frequency, Stop_List)
  real(wp)                  :: Memory(D,MAX_WORDS), Visual(D,MAX_WORDS)
  character(len=WORD_CHARS) :: Word(MAX_WORDS), Stop_List(STOP_WORDS)
  real                   :: Frequency(MAX_WORDS) 
  integer                :: i
  
 Memory = 0.0
 Visual = 0.0
 Frequency = 1.0
 Word(:) = ' ' 
 
  !Loading Stop List words: 
 open(unit=10, file='stoplist', status='old')
 do i = 1, STOP_WORDS
   read(10,*) Stop_List(i)
 enddo
 close(10)
  
 call RandSeed
 open(unit = 1, file = 'wiki.txt', status = 'old', action = 'read')
 open(unit=2, file='output', status='replace')
end subroutine Init
!****************************************************************


!**************************************************************** 
subroutine Compute_Item_Info(Context, W, N, Frequency, loc, Ignore)
  integer   :: N, incl_left, incl_right, i, pos, loc(N)
  real(wp)  :: Context(D,N), W(D,N), log_n
  real      :: Frequency(MAX_WORDS)
  logical   :: Ignore(N)
  
 Context = 0.0 

 do pos = 1, N
  do i = 1, N
    ! Adding entrophy weighted word vector:
    if ((.not.(pos==i)) .and. (.not.(Ignore(i))) .and. (.not.(Ignore(pos)))) then
      !Context(:,pos) = Context(:,pos) + (W(:,i) / (abs(log(Frequency(loc(i))) - log(Frequency(loc(pos))))+1)) 
      Context(:,pos) = Context(:,pos) + W(:,i)
    endif
  enddo
 enddo


end subroutine Compute_Item_Info
!**************************************************************** 



!**************************************************************** 
subroutine Compute_Convolutions(Ownership, Order, W, N)
  integer  :: i, j, k, N, Temp_Conv_Order(N), L, R, Big, Small, pos, N_Conv(N), Max_Conv, kount, &
              reps, CHUNKS
  real(wp) :: Ownership(D), Order(D,N), W(D,N)
  real(wp), allocatable :: In_Vect(:), Temp(:), Temp2(:)
!  real*8 :: In_Vect(D), Temp(D), Temp2(D)
  integer, allocatable :: Order_Vector(:,:,:), n_reps(:,:) 
  logical :: FLAG

 allocate(In_Vect(D), Temp(D), Temp2(D))
 
 CHUNKS = MAX_NGRAM
 if (N < MAX_NGRAM) CHUNKS = N
 
 Max_Conv = 0
 do i = 2, CHUNKS
   Max_Conv = Max_Conv + i
 enddo


 allocate(Order_Vector(Max_Conv,N,N), n_reps(Max_Conv,N))

Order_Vector = 0
N_Conv = 0
n_reps = 0

do i = 1, N
  kount = 0
  do j = 2, CHUNKS
      pos = 1
      do 
        if ((pos+j-2) == N) exit
	FLAG = .false.
        Temp_Conv_Order = 0
	reps = 0	
	do k = pos, (pos+j-1)
	  reps = reps + 1 
	  if (k==i) then
	    if (.not. FLAG) kount = kount + 1
	    FLAG = .true.
	    Temp_Conv_Order(reps) = -1
	  else
	    Temp_Conv_Order(reps) = k
	  endif
	enddo
        if (FLAG) then
	  Order_Vector(kount,:,i) = Temp_Conv_Order
	  n_reps(kount, i) = reps
	endif
	pos = pos + 1
      enddo !  position
  enddo ! j : conv_vector (Max_Conv)
  N_Conv(i) = kount
enddo ! i : word(N)


!Calculating convolutions for this item:
do i = 1, N
  do j = 1, N_Conv(i)
    Temp = 0
    if (Order_Vector(j,1,i) == -1) then 
       Temp2 = Ownership
    else
       Temp2 = W(:,Order_Vector(j,1,i))
    endif
    do k = 2, n_reps(j,i)
       if (Order_Vector(j,k,i) == -1) then 
         In_Vect = Ownership
       else
         In_Vect = W(:,Order_Vector(j,k,i))
       endif
       Temp2 = Scramble(Temp2, P1)
       In_Vect = Scramble(In_Vect, P2)
       !call C06EKF(1,Temp2, In_Vect, D, IFAIL)
       Temp2 = ConvolveDFT(Temp2, In_Vect)
       Temp = Temp2
    enddo
    Order(:,i) = Order(:,i) + Temp
  enddo
enddo

deallocate(Order_Vector, n_reps, Temp, Temp2, In_Vect)

end subroutine Compute_Convolutions
!**************************************************************** 

!****************************************************************
function Scramble(Vector, Scram_Indices)
  real(wp), dimension(D) :: Vector, Temp, Scramble
  integer :: Scram_Indices(D), i 
 do i = 1, D
   Temp(i) = Vector(Scram_Indices(i))
 enddo
 Scramble = Temp

end function Scramble
!****************************************************************



!****************************************************************
subroutine Make_Scramble(P1, P2)
  integer, dimension(D) :: P1, P2
  integer :: i, j, x, id, raw_order(D)

!********P1:
 do i = 1, D
  raw_order(i) = i
 enddo 
 

 do i = 1, D
  do 
    x = int(ran3(id) * D) + 1
    if (raw_order(x) /= -1) goto 10 
  enddo
  10 raw_order(x) = -1
  P1(i) = x
 enddo

!********P2:
 do i = 1, D
  raw_order(i) = i
 enddo 
 
 do i = 1, D
  do 
    x = int(ran3(id) * D) + 1
    if (raw_order(x) /= -1) goto 20 
  enddo
  20 raw_order(x) = -1
  P2(i) = x
 enddo


end subroutine Make_Scramble
!****************************************************************



!****************************************************************
subroutine Store_Matrix(Memory, Visual, Word, Words_Learned, Ownership, P1, P2)
  real(wp)  :: Memory(D,MAX_WORDS), Visual(D,MAX_WORDS), Ownership(D), V1(D), M1(D, MAX_WORDS)
  character(len=WORD_CHARS) :: Word(MAX_WORDS)
  integer :: i, j, Words_Learned, P1(D), P2(D)
  
 open(unit=10, file='word_labels.txt', status='replace')
 open(unit=11, file='matrix.mat', status='replace', form='unformatted')
 open(unit=12, file='ownership.mat', status='replace', form='unformatted')
 open(unit=13, file='visual.mat', status='replace', form='unformatted')
 open(unit=14, file='scrambles.mat', status='replace', form='unformatted')
 
 write(12) Ownership
 write(10, *) Words_Learned
 write(14) P1, P2
 
 do i = 1, Words_Learned
   write(10,'(a20)') Word(i)
   write(11) Memory(:,i)     ! unformatted write; read exactly same
   write(13) Visual(:,i)
 enddo
 
 close(10)
 close(11)
 close(12)
 close(13)
 close(14) 

 
end subroutine Store_Matrix
!****************************************************************


!****************************************************************
subroutine Filter_Item_Info(Memory, Visual, Words_Learned)
  real(wp)  :: Memory(D,MAX_WORDS), Visual(D,MAX_WORDS), V_Temp(D), M_Temp(D)
  integer   :: i, j, Words_Learned
  
 do i = 1, Words_Learned
   V_Temp = Visual(:,i)
   M_Temp = Memory(:,i)
   call Normalize(V_Temp, D)
   call Normalize(M_Temp, D)
   Memory(:,i) = M_Temp - V_Temp 
 enddo

end subroutine 
!****************************************************************

!****************************************************************
! This function returns the circular convolution of the argument
! vectors V1 and V2 with the same dimensionality
!----------------------------------------------------------------
function Convolve(V1, V2)
  integer :: i, j
  real*8, dimension(0:D-1) :: Convolve, sum, V1, V2
 
 sum = 0.0
     
 do i = 0, (D-1)
  do j = 0, (D-1)
    sum(i) = sum(i) + V1(modulo(j,D))*V2(modulo(i-j,D))
  enddo
 enddo
 
 Convolve = sum
 
end function Convolve
!****************************************************************

!****************************************************************
subroutine Wrap_Up
 close(1)
 deallocate(Memory, Visual, Ownership) 
end subroutine Wrap_Up
!****************************************************************


end program PBEAGLE

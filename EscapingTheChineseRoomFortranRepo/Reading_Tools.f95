!======================================================================
MODULE Declare_Constants

implicit none

private

intrinsic KIND 
 integer,   public, parameter           :: wp = KIND(1.0D0)

 integer,   public, parameter :: D              = 2048,      &
				 MAX_WORDS      = 40000,     &
				 SENT_CHARS     = 3000,      &
				 WORD_CHARS     = 20,	     &
				 N_NEIGHBOR     = 6 ,	     &
				 ORD_CANDIDATES = 10,        &
				 TIME_OUT       = 1000000,   &
				 MAX_NGRAM      = 7,         &
				 STOP_WORDS     = 439
 logical,   public            :: E_FILE         = .false.

END MODULE Declare_Constants
!======================================================================



!======================================================================
MODULE Reading_Tools

use Declare_Constants
use number_generators

implicit none

private

public   :: Random_Vector, Vector_Cosine, Normalize, Read_Sentence, Count_Words, &
	    Parse_Sentence, Does_Word_Exist, Create_New_Word, Highest_N_Activated, &
	    Add_Noise

!======================================================================
                   CONTAINS
!======================================================================



!****************************************************************
subroutine Read_Sentence(Sentence, unit_num, BAD_SENT)
  character(len=SENT_CHARS):: Sentence
  character (len=1) :: a, dummy
  integer   :: pos, ios, k, unit_num
  logical   :: BAD_SENT 

Sentence(:) = ' '  
k = 0
 do 
   k = k + 1
   if (k > SENT_CHARS) then
      BAD_SENT = .true.
      read(unit_num, *)
      goto 97
   endif
   read(unit_num,'(a1)', iostat=ios, err=98, end=99, advance='no') a
   if (ios == -2) goto 97
   Sentence(k:k) = a
 enddo

97  return

98  write(*,*) 'Encountered Read Error ', ios
    return

99 E_FILE = .true.

end subroutine Read_Sentence
!****************************************************************


!****************************************************************
subroutine Count_Words(Sentence, BAD_SENT, N)
   integer :: N, pos, w_kount
   character(len=SENT_CHARS) :: Sentence
   character :: a
   logical :: ON_WORD, N_Empty = .false., BAD_SENT 

N = 0
pos = 0
! Counting words: 
do 
 pos = pos + 1
 if (pos > SENT_CHARS) exit
 a = Sentence(pos:pos)
 if (a /= ' ') then 
   ON_WORD = .true.
   N_Empty = .true.
 endif
 if ((ON_WORD) .and. (a == ' ')) then 
   N = N + 1
   ON_WORD = .false.
 endif  
enddo
 if (.not. N_Empty) BAD_SENT = .true.
end subroutine Count_Words
!****************************************************************


!****************************************************************
subroutine Parse_Sentence(Sentence, Sentence_Word, Stop_List, Ignore, N)
   integer :: N, pos, w_kount, k, i, j
   character(len=SENT_CHARS) :: Sentence
   character(len=WORD_CHARS) :: Sentence_Word(N), Stop_List(STOP_WORDS)
   character :: a
   logical :: ON_WORD, Ignore(N)

Sentence_Word(:) = ' '
Ignore = .false.

pos = 0
w_kount = 0
k = 1
do 
 pos = pos + 1
 if (pos > SENT_CHARS) exit
 a = Sentence(pos:pos)
 if (a /= ' ') then
   w_kount = w_kount + 1
   if (w_kount > WORD_CHARS) exit 
   ON_WORD = .true.
   Sentence_Word(k)(w_kount:w_kount) = a
 else 
   if ((ON_WORD) .and. (a == ' ')) then 
      w_kount = 0
      k = k + 1
      ON_WORD = .false.
   endif 
 endif 
enddo


!Filtering out stop words and adjusting N
 do i = 1, STOP_WORDS
  do j = 1, N
    if (Sentence_Word(j) == Stop_List(i)) then 
    Ignore(j) = .true.
    endif
  enddo
 enddo

end subroutine Parse_Sentence
!****************************************************************


!****************************************************************
subroutine Does_Word_Exist(Target_Word, Word, loc, endpoint)
   integer :: loc, endpoint, i,j 
   character(len=WORD_CHARS) :: Target_Word, Word(MAX_WORDS)
   
i = 0
loc = 0   
do
  if (i >= MAX_WORDS) return 
  i = i + 1
  if (Target_Word == Word(i)) then 
    loc = i
    exit
  endif
  if (Word(i) == ' ') then 
    endpoint = i 
    exit
  endif
enddo 
end subroutine Does_Word_Exist   
!****************************************************************

!****************************************************************
subroutine Highest_N_Activated(Context, Memory, High_Resp, Resp_Activation, Words_Learned, num)
  integer  :: i, j, num, High_Resp(num), Words_Learned, loc
  real(wp) :: Memory(D,MAX_WORDS), Context(D), Resp_Activation(num), &
              Cosine_Array(Words_Learned), Sim, Max

do i = 1, Words_Learned
  Cosine_Array(i) = dot_product(Context, Memory(:,i))
enddo

do i = 1, num
  Max = -99.0
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

 do i = 1, num
   Resp_Activation(i) = Vector_Cosine(Context, Memory(:,High_Resp(i)),D)
 enddo

end subroutine Highest_N_Activated
!****************************************************************

!****************************************************************
subroutine Create_New_Word(Target_Word, Word, Vector, Visual, endpoint)
   integer :: loc, endpoint, i,j 
   character(len=WORD_CHARS) :: Target_Word, Word(MAX_WORDS)
   real(wp) :: Visual(D,MAX_WORDS), Vector(D) 

 Word(endpoint) = Target_Word
 Visual(:,endpoint) = Vector

end subroutine Create_New_Word
!****************************************************************


!****************************************************************
function Random_Vector(D, Mu, Sigma)
 integer  :: D
 real(wp) :: Random_Vector(D), Mu, Sigma
 integer  :: i
 
Random_Vector(1) = gaussian(Mu, Sigma) !emptying first value of 0/0
do i = 1, D
  Random_Vector(i) = gaussian(Mu, Sigma)
enddo
end function Random_Vector
!****************************************************************

!****************************************************************
function Add_Noise(Vector, N_Flip, Mu, Sigma)
  integer :: N_Flip, i, j, i_dum
  real(wp):: Vector(D), Temp(D), Add_Noise(D), Mu, Sigma

 Temp = Vector
 do i = 1, N_Flip
  j = int(ran3(i_dum)*D)+1
  Temp(j) = gaussian(Mu, Sigma)  
 enddo
 Add_Noise = Temp

end function Add_Noise
!****************************************************************


!================================================================
!          PROXIMITY ROUTINES:
!================================================================
!****************************************************************
function Vector_Length (Vector, n)
   integer :: i, n 
   real(wp)    :: Vector_Length, Vector(n), SS

 SS = 0.0
  SS = dot_product(Vector, Vector)
 Vector_Length = sqrt(SS)
end function Vector_Length 
!****************************************************************



!****************************************************************
subroutine Normalize (Vector_In, n)
   integer             :: i, n
   real(wp)              :: Vector_In(n), Vect_Length, Test_Val
  
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


!****************************************************************
subroutine Max_Value (Vector, n, location, Max_V)
   real(wp)    :: Max_V, Vector(n), Max
   integer :: count, n, location
  
  location = 0 
  Max = -1.0
  Max_V = 0.0
  do count = 1, n 
    if ((Vector(count)) > Max) then 
      Max = (Vector(count))
      location = count
    end if
  end do 
  
  Max_V = Max
end subroutine Max_Value   
!****************************************************************



END MODULE Reading_Tools




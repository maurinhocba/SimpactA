 CHARACTER(len=*) FUNCTION inttoch(num,digit)
 !=======================================================================
 ! Puts an integer into a string
 !=======================================================================
 IMPLICIT NONE

   INTEGER(kind=4),INTENT(IN):: num    !Number to translate to a string
   INTEGER(kind=4),INTENT(IN):: digit  !Number of digit writed

   !Local variables
   INTEGER (kind=4) :: i, nn, dg

   inttoch = ''                !initializes
   IF (num == 0) inttoch='0'   !case num is 0
   nn = IABS(num)              !valor absoluto

   ! FIRST task, convert number NUM to string of characters
   DO
     IF (nn == 0) EXIT         !
     inttoch = CHAR(MOD(nn,10)+48)//TRIM(inttoch) !add the digit
     nn = nn/10                                   !go to next digid
   END DO
   IF (num < 0) inttoch='-'//TRIM(inttoch)        !add sign

   ! SECOND task, write according to required format

   IF (digit > 0) THEN         !number of digits to use is prescribed
     dg = MIN0(digit,LEN(inttoch))    !minimum number of digits
     nn = LEN_TRIM(inttoch)           !present length of string
     IF (dg > nn) THEN                !if 0´s must be added
       DO i=nn,1,-1
         IF (inttoch(i:i) == '-') THEN  !TRUE only possible when i==1
           nn = nn - 1                  !
           EXIT
         END IF
         inttoch(dg-nn+i:dg-nn+i) = inttoch(i:i) !shift chars
       END DO
       ! add 0s when necessary
       DO i=1,dg-nn
         IF (inttoch(i:i) == '-') CYCLE  !TRUE only possible when i==1
         inttoch(i:i) = '0'
       END DO
     END IF
   END IF

 RETURN
 END FUNCTION inttoch

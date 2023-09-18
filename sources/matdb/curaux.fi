   SUBROUTINE find_cr(chead,c,name,nc,k)
   ! find the associated curve
   !
   IMPLICIT NONE
   TYPE (curve), POINTER :: chead  !curve head positions
   CHARACTER (len=*), INTENT(IN) :: name   !curve name to search
   INTEGER(kind=4), INTENT(IN) :: nc        !number of curves
   INTEGER(kind=4), INTENT(OUT) :: k        !curve position
   TYPE (curve), POINTER :: c               !pointer to curve

   c => chead
   k = 1
   DO
     IF( TRIM(c%name) == TRIM(name) ) EXIT
     k = k+1
     IF(k > nc)THEN
        WRITE(55,"('curve_name ',a,' not found')",ERR=9999) trim(name)
        CALL runen3('FINDCR: curve name not found')
     END IF
     c => c%next
   END DO

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE find_cr

   FUNCTION inte_cr (val,np,t,i)
   ! interpolates in a multi-linear curve
   ! VAL(3,NP) are the tabulated values
   ! T is the abcisa data input
   ! INTE_CR is the output ordinate
   ! I is the interval to begin with
   IMPLICIT NONE
   REAL (kind=8) :: inte_cr
   INTEGER (kind=4), INTENT(IN) :: np
   REAL (kind=8), INTENT(IN) :: val(3,np),t
   INTEGER (kind=4), INTENT(IN OUT) :: i

   i = MIN(i,np-1)
   DO
     IF( val(1,i) <= t )THEN
       IF( val(1,i+1) >= t )THEN
         inte_cr = val(2,i) + (val(2,i+1) - val(2,i))*(t - val(1,i))/ &
                              (val(1,i+1) - val(1,i))
         EXIT
       ELSE
         i = i+1
         IF( i == np ) THEN
           i = i-1
           inte_cr = val(2,np)
           EXIT
         END IF
       END IF
     ELSE
       i= i - 1
       IF( i == 0 ) THEN
         i = 1
         inte_cr = val(2,1)
         EXIT
       END IF
     END IF
   END DO
   RETURN
   END FUNCTION  inte_cr

   SUBROUTINE cur_point( npos,head,point )
   ! point to a curve in a list
   IMPLICIT NONE
   ! dummy arguments
   INTEGER (kind=4), INTENT (IN OUT) :: npos !position in list
   TYPE (curve), POINTER :: head             !(IN)first curve
   TYPE (curve), POINTER :: point            !(OUT)curve
   ! local variables
   INTEGER (Kind=4) :: i

   IF( npos == 0 )RETURN
   point => head
   DO i=1,npos-1
     point => point%next
   END DO
   npos = point%np

   RETURN
   END SUBROUTINE cur_point

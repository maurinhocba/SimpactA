 SUBROUTINE secd07 ( heade, nelem, secs )
 !***********************************************************************
 !
 !*** extract SECTIONs used in the set
 !
 !***********************************************************************
 IMPLICIT NONE
 !     Arguments
 INTEGER (kind=4), INTENT(IN) :: nelem
 INTEGER (kind=4), INTENT(OUT) :: secs(:)

 TYPE (ele07), POINTER  :: heade      !pointer to first element

 INTEGER (kind=4) :: i,j,n,ns,olds
 TYPE (ele07), POINTER  :: e      !pointer to an element

 !initializes
 ns = 1                !one section at least
 secs(1) = heade%matno !first element section
 olds = secs(1)
 e => heade%next       !point to second element

 DO i=2,nelem
   n = e%matno
   IF( n /= olds )THEN
     j = 1
     DO
       IF( n == secs(j) )EXIT
       j = j+1
       IF( j > ns )THEN
         ns = ns+1
         secs(ns) = n
         olds = n
         EXIT
       END IF
     END DO
   END IF
   e => e%next
 END DO

 END SUBROUTINE secd07

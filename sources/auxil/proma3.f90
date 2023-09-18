      SUBROUTINE proma3(a,b,c,n1,n3,zeros)
!********************************************************************
!
!***this routine evaluates a matrix-vector product adding over a
!
!           a(i) = b(i,k) * c(k)    a = b * c
!
!********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: n1,n3
      REAL (kind=8),INTENT(IN) :: b(n1,n3), c(n3)
      REAL (kind=8),INTENT(OUT) :: a(n1)
      LOGICAL(kind=4),INTENT(IN) :: zeros

      IF (zeros) THEN
        a = MATMUL(b,c)
      ELSE
        a = a + MATMUL(b,c)
      END IF
      RETURN

      END SUBROUTINE proma3

      SUBROUTINE proma1(a,b,c,n1,n2,n3)
!********************************************************************
!
!***THIS ROUTINE EVALUATES A MATRIX PRODUCT
!
!           A(I,J) = B(I,K) * C(K,J)    A = B * C
!
!********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: n1,n2,n3
      REAL (kind=8),INTENT(IN) :: b(n1,n3), c(n3,n2)
      REAL (kind=8),INTENT(OUT):: a(n1,n2)

      a = MATMUL(b,c)
      RETURN

      END SUBROUTINE proma1

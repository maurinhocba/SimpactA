      SUBROUTINE proma2(a,b,c,n1,n2,n3)
!********************************************************************
!
!***this routine evaluates a matrix product
!                                                t
!           a(i,j) = b(i,k) * c(j,k)    a = b * c
!
!********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: n1,n2,n3
      REAL (kind=8),INTENT(IN) :: b(1:n1,1:n3), c(1:n2,1:n3)
      REAL (kind=8),INTENT(OUT) :: a(1:n1,1:n2)

      a(1:n1,1:n2) = MATMUL(b(1:n1,1:n3),TRANSPOSE(c(1:n2,1:n3)))
      RETURN

      END SUBROUTINE proma2

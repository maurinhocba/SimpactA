      SUBROUTINE proma4(a,b,c,n1,n2,n3)
!********************************************************************
!
!***this routine evaluates a matrix product
!                                            t
!           a(i,j) = b(k,i) * c(k,j)    a = b * c
!
!********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: n1,n2,n3
      REAL (kind=8),INTENT(IN) :: b(1:n3,1:n1), c(1:n3,1:n3)
      REAL (kind=8),INTENT(OUT) :: a(1:n1,1:n3)

      a(1:n1,1:n2) = MATMUL( TRANSPOSE(b(1:n3,1:n1)), c(1:n3,1:n2) )
      RETURN

      END SUBROUTINE proma4

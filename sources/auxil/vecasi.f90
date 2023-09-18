      SUBROUTINE vecasi(n,v1,v2)
!*****************************************************************************
!
!***  vector assign:    v1(n) -> v2(n)
!
!*****************************************************************************
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: n
      REAL (kind=8),INTENT(IN) :: v1(n)
      REAL (kind=8),INTENT(OUT):: v2(n)

      v2 = v1
      RETURN

      END SUBROUTINE vecasi

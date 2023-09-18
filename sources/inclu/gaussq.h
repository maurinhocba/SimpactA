      SUBROUTINE gaussq    (ngaus ,posgp ,weigp )
!********************************************************************
!
!*** Gauss - Legendre numerical integration constants
!
!********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4), INTENT(IN) :: ngaus
      REAL    (kind=8), INTENT(OUT) :: posgp(ngaus),weigp(ngaus)

      END SUBROUTINE gaussq

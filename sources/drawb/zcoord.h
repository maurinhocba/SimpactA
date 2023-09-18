      SUBROUTINE zcoord (z2, z1, y2, y1, x2, x1, zint, yint, xint)
!-----------------------------------------------
!   m o d u l e s
!-----------------------------------------------
      USE kind_param, ONLY: double
      IMPLICIT NONE
!-----------------------------------------------
!   d u m m y   a r g u m e n t s
!-----------------------------------------------
      REAL (double), INTENT (IN) :: z2
      REAL (double), INTENT (IN) :: z1
      REAL (double), INTENT (IN) :: y2
      REAL (double), INTENT (IN) :: y1
      REAL (double), INTENT (IN) :: x2
      REAL (double), INTENT (IN) :: x1
      REAL (double), INTENT (OUT) :: zint
      REAL (double), INTENT (IN) :: yint
      REAL (double), INTENT (IN) :: xint
!-----------------------------------------------
      END SUBROUTINE zcoord

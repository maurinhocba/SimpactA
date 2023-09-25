      SUBROUTINE linint (x1, y1, x2, y2, x3, y3, x4, y4, xc, yc,            &
     &                   intersect, coincide)

      !   m o d u l e s

      USE kind_param, ONLY: double
      IMPLICIT NONE

      !   d u m m y   a r g u m e n t s

      REAL (double) :: x1
      REAL (double) :: y1
      REAL (double) :: x2
      REAL (double) :: y2
      REAL (double) :: x3
      REAL (double) :: y3
      REAL (double) :: x4
      REAL (double) :: y4
      REAL (double) :: xc
      REAL (double) :: yc
      LOGICAL, INTENT (OUT) :: intersect
      LOGICAL, INTENT (OUT) :: coincide

      END SUBROUTINE linint

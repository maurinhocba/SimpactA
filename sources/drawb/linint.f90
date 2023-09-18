SUBROUTINE linint (x1, y1, x2, y2, x3, y3, x4, y4, xc, yc, &
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

  !   i n t e r f a c e   b l o c k s

  INTERFACE
  INCLUDE 'linequ.h'
  INCLUDE 'insidf.h'
  END INTERFACE

  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s

  REAL (double) :: w, a1, b1, c1, a2, b2, c2, wx, wy
  LOGICAL :: between1, between2
  !-----------------------------------------------

  intersect = .FALSE.
  coincide = .FALSE.

  ! equation of the 1st line: a1*x + b1*y = c1
  !                 2nd line: a2*x + b2*y = c2

  IF (x1 == x2 .AND. y1 == y2) RETURN ! line is reduced to a point

  CALL linequ (x1, y1, x2, y2, a1, b1, c1)
  CALL linequ (x3, y3, x4, y4, a2, b2, c2)

  !     solution of the set of the two equations
  !       a1*x + b1*y = c1
  !       a2*x + b2*y = c2
  !
  !     w, wx, wy - determinants
 
  w  = a1 * b2 - a2 * b1
  wx = c1 * b2 - c2 * b1
  wy = a1 * c2 - a2 * c1

  IF (w /= 0.) THEN         ! not parallel

     wy = a1 * c2 - a2 * c1
     xc = wx / w
     yc = wy / w
     CALL insidf (x1, y1, x2, y2, xc, yc, between1)
     CALL insidf (x3, y3, x4, y4, xc, yc, between2)
     IF (between1 .AND. between2) intersect = .TRUE.

  ELSE IF (w == 0 .AND. wx == 0 .AND. wy == 0) THEN !coinciding

     coincide = .TRUE.

  END IF

  RETURN

END SUBROUTINE linint

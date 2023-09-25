SUBROUTINE insidf (x1, y1, x2, y2, xc, yc, between)
  !---------------------------------------------------------
  !      checks, if (xc,yc) lies between (x1,y1) and (x2,y2)
  !---------------------------------------------------------
  !   m o d u l e s

  USE kind_param, ONLY: double

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  REAL (double), INTENT (IN) :: x1
  REAL (double), INTENT (IN) :: y1
  REAL (double), INTENT (IN) :: x2
  REAL (double), INTENT (IN) :: y2
  REAL (double), INTENT (IN) :: xc
  REAL (double), INTENT (IN) :: yc
  LOGICAL, INTENT (OUT) :: between
  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s

  REAL (double) :: dot
  !-----------------------------------------------

  dot = (x1-xc) * (x2-xc) + (y1-yc) * (y2-yc)
  IF (dot <= 0.) THEN
     between = .TRUE.
  ELSE
     between = .FALSE.
  END IF
  RETURN
END SUBROUTINE insidf

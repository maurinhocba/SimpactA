SUBROUTINE compar (xe1, ye1, xe2, ye2,  xd1, yd1, xd2, yd2, xint, yint)
  !-----------------------------------------------------------------
  !  subroutine to find a common part of the two coinciding segments
  !  (xe1,ye1)-(xe2,ye2)  and (xd1,yd1)-(xd2,yd2)
  !-----------------------------------------------------------------

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  REAL (kind=8), INTENT (IN) :: xe1
  REAL (kind=8), INTENT (IN) :: ye1
  REAL (kind=8), INTENT (IN) :: xe2
  REAL (kind=8), INTENT (IN) :: ye2
  REAL (kind=8), INTENT (IN) :: xd1
  REAL (kind=8), INTENT (IN) :: yd1
  REAL (kind=8), INTENT (IN) :: xd2
  REAL (kind=8), INTENT (IN) :: yd2
  REAL (kind=8), INTENT (OUT) :: xint (2)
  REAL (kind=8), INTENT (OUT) :: yint (2)
  !-----------------------------------------------

  !   l o c a l   v a r i a b l e s

  INTEGER :: lefte, leftd
  REAL (kind=8) :: xle, yle, xre, yre, xld, yld, xrd, yrd
  !------------------------------------------------------

  !  lefte - indicator which point is the left or lower one
  !          for the element edge
  !  leftd - indicator which point is the left or lower one
  !          for the drawbead segment
  !  leftd = 1, 2 ;   leftd = 1, 2
  !          1 or 2 point indicated by lefte or leftd

  lefte = 2
  leftd = 2
  IF (xe1 == xe2) THEN !edge parallel to y-axis
     IF (ye1 < ye2) lefte = 1
     IF (yd1 < yd2) leftd = 1
  ELSE
     IF (xe1 < xe2 .OR. xd1 < xd2) lefte = 1
  END IF

  !  xle - left x extremum of the element edge
  !  yle - y-coord of the left extremum
  !  xre - right x extremum of the element edge
  !  yre - y-coord of the right extremum
  !  xld - left x extremum of the drawbead segment
  !  yld - y-coord of the left extremum
  !  xrd - right x extremum of the drawbead segment
  !  yrd - y-coord of the right extremum

  IF (lefte == 1) THEN
     xle = xe1
     yle = ye1
     xre = xe2
     yre = ye2
  ELSE
     xle = xe2
     yle = ye2
     xre = xe1
     yre = ye1
  END IF
  IF (leftd == 1) THEN
     xld = xd1
     yld = yd1
     xrd = xd2
     yrd = yd2
  ELSE
     xld = xd2
     yld = yd2
     xrd = xd1
     yrd = yd1
  END IF

  !   lefte - indicates the element edge (lefte=2)
  !           or drawbead segment (lefte=1) (the left extremum
  !           of the indicated line will be the left extremum
  !           of their common part
  !   leftd - indicates the element edge (lefte=2)
  !           or drawbead segment (lefte=1) (the right extremum
  !           of the indicated line will be the right extremum
  !           of their common part

  lefte = 2 !2-indicates edge
  leftd = 2
  IF (xe1 == xe2) THEN !edge parallel to y-axis
     IF (yle < yld) lefte = 1
     IF (yre > yrd) leftd = 1
  ELSE
     IF (xle < xld) lefte = 1
     IF (xre > xrd) leftd = 1
  END IF
  IF (lefte == 1) THEN
     xint (1) = xld
     yint (1) = yld
  ELSE
     xint (1) = xle
     yint (1) = yle
  END IF
  IF (leftd == 1) THEN
     xint (2) = xrd
     yint (2) = yrd
  ELSE
     xint (2) = xre
     yint (2) = yre
  END IF
  RETURN
END SUBROUTINE compar

SUBROUTINE qshfun (s, t, h1, h2, h3, h4)

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  REAL (kind=8), INTENT (IN) :: s
  REAL (kind=8), INTENT (IN) :: t
  REAL (kind=8), INTENT (OUT) :: h1
  REAL (kind=8), INTENT (OUT) :: h2
  REAL (kind=8), INTENT (OUT) :: h3
  REAL (kind=8), INTENT (OUT) :: h4
  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s

  REAL (kind=8) :: sp, sm, tp, tm
  !-----------------------------------------------

  tp = .25 * (1.0+t)
  tm = .25 * (1.0-t)
  sp = 1.0 + s
  sm = 1.0 - s
  h1 = tm * sm
  h2 = tm * sp
  h3 = tp * sp
  h4 = tp * sm

  RETURN
END SUBROUTINE qshfun

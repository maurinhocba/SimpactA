      SUBROUTINE proje3 (xs, ys, x,   y,  z, n,  sn, tn)

      !   m o d u l e s

      USE kind_param, ONLY: double

      !   i n t e r f a c e   b l o c k s

      IMPLICIT NONE

      !   d u m m y   a r g u m e n t s

      REAL (double), INTENT (IN) :: xs
      REAL (double), INTENT (IN) :: ys
      REAL (double), INTENT (INOUT) :: sn
      REAL (double), INTENT (INOUT) :: tn
      REAL (double), INTENT (IN) :: x (4)
      REAL (double), INTENT (IN) :: y (4)
      REAL (double), INTENT (IN) :: z (4)
      REAL (double) :: n (3)
      !-----------------------------------------------
      !   l o c a l   v a r i a b l e s

      REAL (double), DIMENSION (3) :: a, b
      REAL (double) :: area, a1, a2, b1, b2
      !-----------------------------------------------

      a1 = x (2) - x (1)
      a2 = y (2) - y (1)
      b1 = x (3) - x (1)
      b2 = y (3) - y (1)
      area = a1 * b2 - a2 * b1

      a1 = x (3) - x (2)
      a2 = y (3) - y (2)
      b1 = xs - x (2)
      b2 = ys - y (2)
      sn = a1 * b2 - a2 * b1

      a1 = x (1) - x (3)
      a2 = y (1) - y (3)
      b1 = xs - x (3)
      b2 = ys - y (3)
      tn = a1 * b2 - a2 * b1

      sn = sn / area
      tn = tn / area

      a (1) = x (2) - x (1)
      a (2) = y (2) - y (1)
      a (3) = z (2) - z (1)
      b (1) = x (3) - x (1)
      b (2) = y (3) - y (1)
      b (3) = z (3) - z (1)
      CALL vecpro (a, b, n)
      CALL vecuni (3, n, area)

      RETURN
      END SUBROUTINE proje3

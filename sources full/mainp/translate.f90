      SUBROUTINE translate (xx1, xx2, xx3, xs, ys, zs)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE kind_param, ONLY: double
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL (double), INTENT (INOUT) :: xs
      REAL (double), INTENT (INOUT) :: ys
      REAL (double), INTENT (INOUT) :: zs
      REAL (double), INTENT (INOUT) :: xx1 (6)
      REAL (double), INTENT (INOUT) :: xx2 (6)
      REAL (double), INTENT (INOUT) :: xx3 (6)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i
      REAL (double) :: toler, x0, y0, z0, scale, dd, dist2
      REAL (double), DIMENSION (3) :: a, b
      REAL (double) :: dot, toler2
      REAL (double), DIMENSION (3) :: c
      REAL (double) :: xl, pxs, pys, pzs, toler3
!-----------------------------------------------
      DATA toler / 1.E-10 /
      DATA toler2 / 1.E-10 /
      DATA toler3 / 1.E-15 /
      xx1 (5) = xx1 (1)
      xx2 (5) = xx2 (1)
      xx3 (5) = xx3 (1)
      xx1 (6) = xx1 (2)
      xx2 (6) = xx2 (2)
      xx3 (6) = xx3 (2)
      DO i = 1, 4
        dd = dist2 (xx1(i), xx1(i+1), xx2(i), xx2(i+1), xx3(i),xx3(i+1))
        IF (dd == 0.) THEN
           dd = dist2 (xx1(i), xs, xx2(i), ys, xx3(i), zs)
           IF (Abs(dd) < toler3) THEN
              x0 = (xx1(1)+xx1(2)+xx1(3)+xx1(4)) / 4.
              y0 = (xx2(1)+xx2(2)+xx2(3)+xx2(4)) / 4.
              z0 = (xx3(1)+xx3(2)+xx3(3)+xx3(4)) / 4.
              scale = 1.D0 - toler
              xs = x0 + scale * (xs-x0)
              ys = y0 + scale * (ys-y0)
              zs = z0 + scale * (zs-z0)
           ELSE
              a (1) = xx1 (i+2) - xx1 (i)
              a (2) = xx2 (i+2) - xx2 (i)
              a (3) = xx3 (i+2) - xx3 (i)
              b (1) = xs - xx1 (i)
              b (2) = ys - xx2 (i)
              b (3) = zs - xx3 (i)
              CALL vecasi (3, b, c)
              CALL vecuni (3, a, xl)
              CALL vecuni (3, c, xl)
              dot = DOT_PRODUCT (a, c)
              IF (Abs(dot) < toler2) THEN
                 x0 = (xx1(1)+xx1(2)+xx1(3)+xx1(4)) / 4.
                 y0 = (xx2(1)+xx2(2)+xx2(3)+xx2(4)) / 4.
                 z0 = (xx3(1)+xx3(2)+xx3(3)+xx3(4)) / 4.
                 pxs = toler * (x0-xx1(i))
                 pys = toler * (y0-xx2(i))
                 pzs = toler * (z0-xx3(i))
                 xs = xx1 (i) + b (1) + pxs
                 ys = xx2 (i) + b (2) + pys
                 zs = xx3 (i) + b (3) + pzs
              END IF
           END IF
           RETURN
        END IF
      END DO
      RETURN
      END SUBROUTINE translate
!
      REAL (KIND(0.0d0)) FUNCTION dist2 (x1, x2, y1, y2, z1, z2)
      DOUBLE PRECISION :: x1
      DOUBLE PRECISION :: x2
      DOUBLE PRECISION :: y1
      DOUBLE PRECISION :: y2
      DOUBLE PRECISION :: z1
      DOUBLE PRECISION :: z2
      dist2 = (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) + (z1-z2) * (z1-z2)
      RETURN
      END FUNCTION dist2

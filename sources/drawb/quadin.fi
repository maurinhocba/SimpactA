SUBROUTINE quadin (xx1, xx2,  xx3,  xl1,  yl1,   xl2, &
                     yl2, xint, yint, zint, lcros, alpha)

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT (INOUT) :: lcros
  REAL (kind=8) :: xl1
  REAL (kind=8) :: yl1
  REAL (kind=8) :: xl2
  REAL (kind=8) :: yl2
  REAL (kind=8), INTENT (OUT) :: alpha
  REAL (kind=8) :: xx1 (5)
  REAL (kind=8) :: xx2 (5)
  REAL (kind=8) :: xx3 (5)
  REAL (kind=8) :: xint (2)
  REAL (kind=8) :: yint (2)
  REAL (kind=8) :: zint (2)

  !   i n t e r f a c e   b l o c k s
  INTERFACE
    INCLUDE 'linint.h'
    INCLUDE 'zcoord.h'
  END INTERFACE
  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s

  INTEGER :: i, k
  REAL (kind=8) :: dist2
  LOGICAL :: intersect, coincide
  !-----------------------------------------------

  alpha = 1.
  DO i = 1, 4
     CALL linint (xx1(i), xx2(i), xx1(i+1), xx2(i+1), xl1,  &
                  yl1,  xl2, yl2, xint(lcros+1), yint(lcros+1),  &
                  intersect, coincide)
     IF (coincide) THEN
        alpha = 0.5
        lcros = 2
        CALL compar (xx1(i), xx2(i), xx1(i+1), xx2(i+1), xl1,  &
                     yl1, xl2, yl2, xint, yint)
        DO k = 1, 2
           CALL zcoord (xx3(i+1), xx3(i), xx2(i+1), xx2(i),xx1(i+1), &
                        xx1(i),   zint(k),yint(k), xint(k))
        END DO
        RETURN
     END IF
     IF (intersect) THEN
        IF (lcros == 1) THEN

        ! check if the 2nd intersection point coincides with the 1st one

           dist2 = (xint(lcros)-xint(lcros+1)) ** 2 +               &
                   (yint(lcros)-yint(lcros+1)) ** 2

           IF (dist2 > 0.) THEN
              lcros = lcros + 1
              CALL zcoord  &
                   (xx3(i+1), xx3(i), xx2(i+1), xx2(i), xx1(i+1),  &
                   xx1(i), zint(lcros), yint(lcros), xint(lcros))
           END IF
        ELSE
           lcros = lcros + 1
           CALL zcoord  &
                (xx3(i+1), xx3(i), xx2(i+1), xx2(i), xx1(i+1), &
                 xx1(i), zint(lcros), yint(lcros), xint(lcros))
        END IF
     END IF
     IF (lcros == 2) RETURN
  END DO
  RETURN
END SUBROUTINE quadin

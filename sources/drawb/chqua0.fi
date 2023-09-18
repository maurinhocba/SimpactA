SUBROUTINE chqua0 (x, xl, nstatold, ldbs, irect, ndbpoin, fdat, nnode)
  !-------------------------------------------------
  !   m o d u l e s

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT (OUT) :: nstatold
  INTEGER, INTENT (IN) :: ndbpoin
  INTEGER, INTENT (IN) :: nnode
  INTEGER, INTENT (OUT) :: ldbs (2)
  INTEGER, INTENT (IN) :: irect (4)
  REAL (kind=8), INTENT (IN) :: x (:,:)
  REAL (kind=8), INTENT (IN) :: xl (:, :)
  REAL (kind=8), INTENT (OUT) :: fdat (3)
  !-----------------------------------------------

  !   l o c a l   v a r i a b l e s

  INTEGER :: i, nn, llow, lupp, lcros, ierr, lcrosold, jj
  REAL (kind=8), DIMENSION (5) :: xx1, xx2, xx3
  REAL (kind=8) :: alpha, xqmax, yqmax, xqmin, yqmin, xdmax, ydmax, &
                   xdmin, ydmin, xl1, xl2, yl1, yl2
  REAL (kind=8), DIMENSION (2) :: xint, yint, zint
  REAL (kind=8) :: xdmid, ydmid, zdmid, ss, tt
  REAL (kind=8), DIMENSION (3) :: n
  LOGICAL :: outside
  !-----------------------------------------------

  DO i = 1, 5
     jj = i / 5 + Mod (i, 5)
     nn = irect (jj)
     xx1 (i) = x (1, nn)
     xx2 (i) = x (2, nn)
     xx3 (i) = x (3, nn)
  END DO
  xqmax = dmax1 (xx1(1), xx1(2), xx1(3), xx1(4))
  yqmax = dmax1 (xx2(1), xx2(2), xx2(3), xx2(4))
  xqmin = dmin1 (xx1(1), xx1(2), xx1(3), xx1(4))
  yqmin = dmin1 (xx2(1), xx2(2), xx2(3), xx2(4))

  llow = 1
  lupp = ndbpoin - 1
  lcros = 0

  DO i = llow, lupp
     IF (lcros /= 2) THEN
        xl1 = xl (1, i)
        yl1 = xl (2, i)
        xl2 = xl (1, i+1)
        yl2 = xl (2, i+1)
        xdmax = dmax1 (xl1, xl2)
        ydmax = dmax1 (yl1, yl2)
        xdmin = dmin1 (xl1, xl2)
        ydmin = dmin1 (yl1, yl2)

        ! preliminary checking of a possibility of intersection

        outside = xdmin > xqmax .OR. ydmin > yqmax .OR. xdmax <     &
                  xqmin .OR. xdmax < xqmin
        IF ( .NOT. outside .AND. lcros /= 2) THEN
           lcrosold = lcros

           ! exact checking of intersection

           CALL quadin (xx1, xx2,  xx3,  xl1,  yl1,   xl2, &
                        yl2, xint, yint, zint, lcros, alpha)
           IF (lcros == lcrosold+1) ldbs (lcros) = i
           IF (lcros == lcrosold+2) THEN
              ldbs (1) = i
              ldbs (2) = i
           END IF
        END IF !if(.not.outside)then
     END IF !if(lcros/=2)then
  END DO
  IF (lcros == 2) THEN

     !  intersection found out

     nstatold = 1
     xdmid = (xint(1)+xint(2)) / 2.
     ydmid = (yint(1)+yint(2)) / 2.
     zdmid = (zint(1)+zint(2)) / 2.

     !  projection of the midpoint of the intersecting line

     IF (nnode == 3) THEN
        CALL proje3 (xdmid, ydmid, xx1,  xx2, xx3, &
                     n,     ss,    tt)
     ELSE IF (nnode == 4) THEN
        CALL projec (xdmid, ydmid, zdmid, xx1,  xx2, xx3, &
                     n,     ss,    tt,    ierr)
     END IF
     fdat (1) = ss
     fdat (2) = tt
     fdat (3) = 0.
  ELSE

     ! intersection not found

     nstatold = 0
     fdat (1) = 0.
     fdat (2) = 0.
     fdat (3) = 0.
  END IF
  RETURN
END SUBROUTINE chqua0

SUBROUTINE chquad (x,     xl,      nstatnew, nstatold, ldbs,  &
                   irect, ndbpoin, fdat, daz, elmodul,  fmax, &
                   fcont,   tdbf,     ihflag,   foldh,  nnode)
  !------------------------------------------------------------

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT (OUT) :: nstatnew
  INTEGER, INTENT (IN) :: nstatold
  INTEGER, INTENT (IN) :: ndbpoin
  INTEGER, INTENT (IN) :: ihflag
  INTEGER, INTENT (IN) :: nnode
  REAL (kind=8), INTENT (IN) :: elmodul
  REAL (kind=8), INTENT (IN) :: fmax
  REAL (kind=8), INTENT (IN) :: foldh
  INTEGER, INTENT (INOUT) :: ldbs (2)
  INTEGER, INTENT (IN) :: irect (4)
  REAL (kind=8), INTENT (IN) :: x (:,:)
  REAL (kind=8) :: xl (:, :)
  REAL (kind=8) :: fdat (3)
  REAL (kind=8), INTENT (INOUT) :: fcont (:, :)
  REAL (kind=8), INTENT (INOUT) :: tdbf (4)
  LOGICAL,INTENT(INOUT):: daz


  !-----------------------------------------------

  !   l o c a l   v a r i a b l e s

  LOGICAL :: outside
  INTEGER :: i, nn, llow, lupp, lcros, ierr, nn1, nn2, nn3, nn4,    &
             jj, lcrosold, ld
  REAL (kind=8), DIMENSION (2) :: xint, yint, zint
  REAL (kind=8), DIMENSION (3) :: n, s, t
  REAL (kind=8), DIMENSION (4) :: h
  REAL (kind=8), DIMENSION (5) :: xx1, xx2, xx3
  REAL (kind=8) :: alpha, dist2, dist1, dx, dy, dz, ds, frold,  &
                   force, forcex, forcey,  forcez, h1, h2, h3, h4, &
                   ss, tt, xqmax, yqmax, xqmin, yqmin, xdmax, ydmax, &
                   xdmin, ydmin, xdmid, ydmid, zdmid, xl1, xl2, yl1, &
                   yl2, xlc, xold, yold, zold
  !-----------------------------------------------------------------

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
  IF (nstatold == 1) THEN
     llow = max0 (1, ldbs(1)-2)
     lupp = min0 (ndbpoin-1, ldbs(2)+2)
  ELSE
     llow = 1
     lupp = ndbpoin - 1
  END IF
  lcros = 0

  ! checking of the start point of the drawbead line

  IF (llow == 1) THEN
     IF (nnode == 3) THEN
        CALL proje3 (xl(1, 1), xl(2, 1), xx1,  xx2, xx3, &
                     n,        ss,       tt)
        outside = ss < 0.d0 .OR. ss > 1.d0 .OR. tt < 0.d0 .OR. tt > &
        1.d0 .OR. 1.d0 - ss - tt < 0.d0 .OR. 1.d0 - ss - tt > 1.d0
     ELSE IF (nnode == 4) THEN
        CALL projec (xl(1, 1), xl(2, 1), 0.d0, xx1, xx2, xx3, &
                     n,  ss, tt, ierr)
        outside = Abs (ss) > 1.d0 .OR. Abs (tt) > 1.d0
     END IF
     IF ( .NOT. outside) THEN
        lcros = 1
        ldbs (1) = 1
        alpha = 1.
        xint (1) = xl (1, 1)
        yint (1) = xl (2, 1)
        IF (nnode == 4) THEN
           CALL qshfun (ss, tt, h1, h2, h3, h4)
           zint (1) = h1*xx3(1) + h2*xx3(2) + h3*xx3(3) + h4*xx3(4)
        ELSE
           zint (1) = ss*xx3(1) + tt*xx3(2) + (1.-ss-tt)*xx3(3)
        END IF
     END IF
  END IF
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

        !  preliminary checking of a possibility of intersection

        outside = xdmin > xqmax .OR. ydmin > yqmax .OR. xdmax <     &
        xqmin .OR. xdmax < xqmin
        IF ( .NOT. outside .AND. lcros /= 2) THEN
           lcrosold = lcros

           ! exact checking of intersection

           CALL quadin (xx1, xx2,  xx3,  xl1,  yl1,   xl2,  &
                        yl2, xint, yint, zint, lcros, alpha)
           IF (lcros == lcrosold+1) ldbs (lcros) = i
           IF (lcros == lcrosold+2) THEN
              ldbs (1) = i
              ldbs (2) = i
           END IF
        END IF !if(.not.outside)then
     END IF !if(lcros/=2)then
  END DO

  ! checking of the end point of the drawbead line

  IF (lcros == 1 .AND. lupp == ndbpoin-1) THEN
     IF (nnode == 3) THEN
        CALL proje3 (xl(1, ndbpoin), xl(2, ndbpoin), xx1, xx2, xx3, &
              n, ss, tt)
        outside = ss < 0.d0 .OR. ss > 1.d0 .OR. tt < 0.d0 .OR. tt > &
        1.d0 .OR. 1.d0 - ss - tt < 0.d0 .OR. 1.d0 - ss - tt > 1.d0
     ELSE IF (nnode == 4) THEN
        CALL projec (xl(1, ndbpoin), xl(2, ndbpoin), 0d0, xx1, xx2, xx3, &
              n, ss, tt, ierr)
        outside = Abs (ss) > 1.d0 .OR. Abs (tt) > 1.d0
     END IF
     IF ( .NOT. outside) THEN

        ! check if the 2nd intersection point coincides with the first one

        dist2 = (xint(lcros)-xl(1, ndbpoin)) ** 2 +                 &
                (yint(lcros)-xl(2, ndbpoin)) ** 2
        IF (dist2 > 0.) THEN
           lcros = 2
           ldbs (2) = ndbpoin - 1
           alpha = 1.
           xint (2) = xl (1, ndbpoin)
           yint (2) = xl (2, ndbpoin)
           IF (nnode == 4) THEN
              CALL qshfun (ss, tt, h1, h2, h3, h4)
              zint (2) = h1 * xx3 (1) + h2 * xx3 (2) + h3 * xx3 (3) &
                       + h4 * xx3 (4)
           ELSE
              zint (2) = ss * xx3 (1) + tt * xx3 (2) + (1.-ss-tt) * &
                         xx3 (3)
           END IF
        END IF
     END IF
  END IF
  IF (lcros == 2) THEN

     ! intersection found out

     nstatnew = 1
     xdmid = (xint(1)+xint(2)) / 2.
     ydmid = (yint(1)+yint(2)) / 2.
     zdmid = (zint(1)+zint(2)) / 2.

     ! projection of the midpoint of the intersecting line

     IF (nnode == 3) THEN
        CALL proje3 (xdmid, ydmid, xx1, xx2, xx3, n, ss, tt)
     ELSE IF (nnode == 4) THEN
        CALL projec (xdmid, ydmid, zdmid, xx1, xx2, xx3, n, ss, tt, ierr)
     END IF

     ! local system: (n,t,s)  s = n x t
     ! t is directed in accordance with the numbering order of the
     ! drawbead line

     IF (ldbs(1) <= ldbs(2)) THEN
        ld = ldbs (1)
     ELSE
        ld = ldbs (2)
     END IF
     dist1 = (xint(1)-xl(1, ld)) ** 2 + (yint(1)-xl(2, ld)) ** 2
     dist2 = (xint(2)-xl(1, ld)) ** 2 + (yint(2)-xl(2, ld)) ** 2
     IF (dist1 < dist2) THEN
        t (1) = xint (2) - xint (1)
        t (2) = yint (2) - yint (1)
        t (3) = zint (2) - zint (1)
     ELSE
        t (1) = xint (1) - xint (2)
        t (2) = yint (1) - yint (2)
        t (3) = zint (1) - zint (2)
     END IF
     CALL vecuni (3, t, xlc)
     CALL vecpro (n, t, s)
     IF (nnode == 4) THEN
        CALL qshfun (ss, tt, h(1), h(2), h(3), h(4))
     ELSE IF (nnode == 3) THEN
        h (1) = ss
        h (2) = tt
        h (3) = 1. - h (1) - h (2)
     END IF

     IF (nstatold == 0) THEN !for the first time intersected
        frold = foldh
     ELSE
        IF (nnode == 4) CALL qshfun (fdat(1), fdat(2), h1, h2, h3, h4)
        nn1 = irect (1)
        nn2 = irect (2)
        nn3 = irect (3)
        nn4 = irect (4)
        IF (nnode == 4) THEN
           xold = h1 * x (1, nn1) + h2 * x (1, nn2) + h3 * x (1,nn3) + h4 * x (1, nn4)
           yold = h1 * x (2, nn1) + h2 * x (2, nn2) + h3 * x (2,nn3) + h4 * x (2, nn4)
           zold = h1 * x (3, nn1) + h2 * x (3, nn2) + h3 * x (3,nn3) + h4 * x (3, nn4)
        ELSE IF (nnode == 3) THEN
           h1 = fdat (1)
           h2 = fdat (2)
           h3 = 1. - h1 - h2
           xold = h1 * x (1, nn1) + h2 * x (1, nn2) + h3 * x (1,nn3)
           yold = h1 * x (2, nn1) + h2 * x (2, nn2) + h3 * x (2,nn3)
           zold = h1 * x (3, nn1) + h2 * x (3, nn2) + h3 * x (3,nn3)
        END IF
        dx = xdmid - xold
        dy = ydmid - yold
        dz = zdmid - zold
        ds = s (1) * dx + s (2) * dy + dz * s (3)

        ! restraining force on the unit length

        frold = fdat (3)
        force = ds * elmodul
        frold = frold + force
        IF (Abs(frold) > Abs(fmax)) frold = dsign (1.0d0, frold) * Abs (fmax)
     END IF

     ! restraining force for the element

     force = frold * alpha * xlc
     forcex = force * s (1)
     forcey = force * s (2)
     forcez = force * s (3)
     DO jj = 1, nnode
        nn = irect (jj)
        fcont (1, nn) = fcont (1, nn) - forcex * h (jj)
        fcont (2, nn) = fcont (2, nn) - forcey * h (jj)
        fcont (3, nn) = fcont (3, nn) - forcez * h (jj)
     END DO
     daz = .TRUE.   !Drawbead affected zone
     fdat (1) = ss
     fdat (2) = tt
     fdat (3) = frold

     ! total restraining force for the drawbead line

     IF (ihflag /= 0) THEN
        tdbf (1) = tdbf (1) + forcex
        tdbf (2) = tdbf (2) + forcey
        tdbf (3) = tdbf (3) + forcez
        tdbf (4) = tdbf (4) + force
     END IF
  ELSE

     ! intersection not found

     nstatnew = 0
     fdat (1) = 0.
     fdat (2) = 0.
     fdat (3) = 0.
  END IF
  RETURN

END SUBROUTINE chquad

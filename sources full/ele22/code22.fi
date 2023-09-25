SUBROUTINE code22(nnode,nelem,ndime,dtime,lnods,matno,coord)

  !  determines the maximum time increment for 3/4 2D thermal element

  IMPLICIT NONE

  INTEGER (kind=4) :: nnode,nelem,ndime,lnods(:,:),matno(:)
  REAL (kind=8) :: dtime,coord(:,:)

  INTEGER*4 j,ielem,isec,osec,ln(nnode),nvert

  REAL (kind=8) :: cp,cond,dens,deltm,xlong2,xl(12),x(3,4),hmin
  TYPE (section), POINTER :: sec
  TYPE (mater  ), POINTER :: mat

  IF( nnode == 8) nvert = 12
  IF( nnode == 4) nvert = 6
  osec = -1
  DO ielem=1,nelem

    isec = matno(ielem)
    IF( isec /= osec )THEN
       sec => psecs(isec)%p
       mat => sec%mtbas
       dens   = mat%prope(5)          !mass density (not scaled)
       cp     = mat%prope(7)          !heat capacity
       cond   = mat%prope(8)          !heat conductivity
      osec = isec
    END IF

    xl(1:nvert) = 0d0
    IF(nnode == 8) THEN
      ln(1:nnode) = lnods(1:nnode,ielem)
      DO j = 1,ndime
        xl( 1) = xl( 1) + (coord(j,ln(1))-coord(j,ln(2)))**2
        xl( 2) = xl( 2) + (coord(j,ln(2))-coord(j,ln(3)))**2
        xl( 3) = xl( 3) + (coord(j,ln(3))-coord(j,ln(4)))**2
        xl( 4) = xl( 4) + (coord(j,ln(1))-coord(j,ln(4)))**2
        xl( 5) = xl( 5) + (coord(j,ln(5))-coord(j,ln(6)))**2
        xl( 6) = xl( 6) + (coord(j,ln(6))-coord(j,ln(7)))**2
        xl( 7) = xl( 7) + (coord(j,ln(7))-coord(j,ln(8)))**2
        xl( 8) = xl( 8) + (coord(j,ln(5))-coord(j,ln(8)))**2
        xl( 9) = xl( 9) + (coord(j,ln(1))-coord(j,ln(5)))**2
        xl(10) = xl(10) + (coord(j,ln(2))-coord(j,ln(6)))**2
        xl(11) = xl(11) + (coord(j,ln(3))-coord(j,ln(7)))**2
        xl(12) = xl(12) + (coord(j,ln(4))-coord(j,ln(8)))**2
      END DO
      xlong2 = MINVAL(xl(1:nvert), mask= xl(1:nvert) > 0d0)
    ELSE IF(nnode == 4) THEN
      ln(1:nnode) = lnods(1:nnode,ielem)
      x = coord(1:3,ln(1:4))
      xlong2 = hmin (x)
      xlong2 = xlong2*xlong2
    END IF

    IF (cond /= 0d0) THEN
      deltm=0.5*xlong2*(dens*cp)/cond
      IF(deltm < dtime) dtime=deltm
    ENDIF
  END DO

END SUBROUTINE code22

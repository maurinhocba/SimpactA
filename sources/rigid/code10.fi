SUBROUTINE code10(nelem,nnode,ndime,lnods,matno,dtime,coord)

  !  determines the maximum time increment for thermal element

  USE ctrl_db, ONLY : tmscal
  IMPLICIT NONE

  INTEGER (kind=4) :: nnode,nelem,ndime,lnods(:,:),matno(:)
  REAL (kind=8) :: dtime,coord(:,:)

  INTEGER(kind=4):: j,ielem,isec,osec,ln(nnode),nvert

  REAL (kind=8) :: deltm,xlong2,xl(12),x(3,nnode),hmin, &
                   densi,capac,condu,fac

  TYPE (section), POINTER :: sec

  IF( ndime == 2 )THEN
    IF( nnode == 3) nvert = 3
    IF( nnode == 4) nvert = 6
  ELSE
    IF( nnode == 4) nvert = 4
    IF( nnode == 8) nvert = 12
  END IF

  osec = -1
  DO ielem=1,nelem

    isec = matno(ielem)                !associated material
    IF( isec /= osec )THEN             !present and previous material diff
      sec => psecs(isec)%p
      ! properties may be temperature dependant
      densi  = sec%mtbas%prope(5)*tmscal   !mass density (scaled)
      capac  = sec%mtbas%prope(7)      !heat capacity
      condu  = sec%mtbas%prope(8)      !heat conductivity
      fac    = capac*densi/condu       !factor
      osec  = isec                     !update previous material
    END IF

    ln(1:nnode) = lnods(1:nnode,ielem)
    x(1:ndime,1:nnode) = coord(1:ndime,lnods(1:nnode,ielem))
    xl(1:nvert) = 0d0
    IF( ndime == 2 )THEN
      IF(nnode == 4) THEN
        DO j = 1,ndime
          xl( 1) = xl( 1) + (x(j,1)-x(j,2))**2
          xl( 2) = xl( 2) + (x(j,2)-x(j,3))**2
          xl( 3) = xl( 3) + (x(j,3)-x(j,4))**2
          xl( 4) = xl( 4) + (x(j,1)-x(j,4))**2
          xl( 5) = xl( 5) + (x(j,1)-x(j,3))**2
          xl( 6) = xl( 6) + (x(j,2)-x(j,4))**2
        END DO
        xlong2 = MINVAL(xl(1:nvert), mask= xl(1:nvert) > 0d0)
      ELSE IF(nnode == 3) THEN
        DO j = 1,ndime
          xl( 1) = xl( 1) + (x(j,1)-x(j,2))**2
          xl( 2) = xl( 2) + (x(j,2)-x(j,3))**2
          xl( 3) = xl( 3) + (x(j,3)-x(j,1))**2
        END DO
        xlong2 = MINVAL(xl(1:nvert), mask= xl(1:nvert) > 0d0)
      END IF
    ELSE
      IF(nnode == 8) THEN
        DO j = 1,ndime
          xl( 1) = xl( 1) + (x(j,1)-x(j,2))**2
          xl( 2) = xl( 2) + (x(j,2)-x(j,3))**2
          xl( 3) = xl( 3) + (x(j,3)-x(j,4))**2
          xl( 4) = xl( 4) + (x(j,1)-x(j,4))**2
          xl( 5) = xl( 5) + (x(j,5)-x(j,6))**2
          xl( 6) = xl( 6) + (x(j,6)-x(j,7))**2
          xl( 7) = xl( 7) + (x(j,7)-x(j,8))**2
          xl( 8) = xl( 8) + (x(j,5)-x(j,8))**2
          xl( 9) = xl( 9) + (x(j,1)-x(j,5))**2
          xl(10) = xl(10) + (x(j,2)-x(j,6))**2
          xl(11) = xl(11) + (x(j,3)-x(j,7))**2
          xl(12) = xl(12) + (x(j,4)-x(j,8))**2
        END DO
        xlong2 = MINVAL(xl(1:nvert), mask= xl(1:nvert) > 0d0)
      ELSE IF(nnode == 4) THEN
        xlong2 = hmin (x)
        xlong2 = xlong2*xlong2
      END IF
    END IF

    IF (condu /= 0d0) THEN
      deltm=0.5*xlong2*(densi*capac)/condu
      IF(deltm < dtime) dtime=deltm
    ENDIF
  END DO

END SUBROUTINE code10

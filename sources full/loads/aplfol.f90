SUBROUTINE aplfol (ndime,coora,resid,lcurv,dtime,ttime,ntype,fltype,flparm,headf,factor)
!     applies follower load
USE loa_db
IMPLICIT NONE

  !--- Dummy varialbes
  CHARACTER(len=*):: fltype
  INTEGER(kind=4):: lcurv,ndime,ntype
  REAL(kind=8):: dtime,ttime,resid(:,:),coora(:,:),factor
  TYPE(flpar):: flparm
  TYPE(foll_seg),POINTER:: headf
  !--- Function
  REAL(kind=8):: functs
  !--- Local variables
  REAL (kind=8), PARAMETER  :: pi=3.141592653589793
  LOGICAL:: dummy
  INTEGER(kind=4):: n,inode,nnode,nodcn
  REAL(kind=8):: area2,dp,dvol,facts,facts0,pforc,press,q,temp,temp0,   &
                 timin,vol,px(3),x(ndime,4),xcen(3),trans(3,3),vn(3)
  REAL (kind=8) l1(3),l2(3),t1(3),t3(3)
  TYPE (foll_seg),POINTER :: seg
   REAL(kind=8):: x0,y0,gamma,xc,h,xl,xr,w,dd,dy,rr
   LOGICAL :: water,full,axisy
  !---------------------------------------------------------------

  water = .FALSE.
  facts = functs(lcurv,ttime)*factor
  timin = 0d0 !default

  IF (TRIM(fltype) == 'INFLAT' .OR. TRIM(fltype) == 'TUBHYD') THEN
    timin = flparm%timini
    nodcn = flparm%nodcen
    xcen(1:ndime) = coora(1:ndime,nodcn)
    vol   = 0d0
    seg => headf
    DO
      IF (.NOT.ASSOCIATED(seg)) EXIT
      nnode = seg%nnode
      x(1:ndime,1:nnode) = coora(1:ndime,seg%lnofl(1:nnode))
      IF ( nnode == 2) THEN
        CALL vol2d(dvol,x(1:2,1:2),xcen,ntype)
      ELSE IF ( nnode == 3) THEN
        CALL tetr(dvol,x,xcen)
      ELSE IF (nnode == 4) THEN
        CALL pyrmd4 (dvol, x, xcen)
      END IF
      vol = vol + dvol
      seg => seg%next
    END DO

    IF (ttime >= timin ) THEN
      IF (TRIM(fltype) == 'INFLAT') THEN
        temp0 = functs(flparm%lc,0d0)
        temp  = functs(flparm%lc,ttime)
        IF (flparm%volin == 0d0 ) flparm%volin = vol
        facts0 = functs(lcurv,timin)*factor
      ELSE IF (TRIM(fltype) == 'TUBHYD') THEN
        q = functs(flparm%lc,ttime)*factor
        q = q * flparm%qref
        IF (flparm%vol == 0d0 .AND. ttime >= timin) flparm%vol = vol
        dvol = vol - flparm%vol
        dp = flparm%kvol * (q*dtime - dvol) / vol
      END IF
    END IF

  ELSE IF( TRIM(fltype) == 'WATERB' )THEN
    water = .TRUE.
    gamma = flparm%qref   !reference specific weight
    y0 = flparm%press     !reference heigth for water pressure
    full = .TRUE.         !bag is full of water
    rr = 1d0              !Initializes radius
    axisy = MOD(flparm%nodcen,10) == 3
    IF(flparm%nodcen < 4 )THEN    !
      w = flparm%pext       !width at free surface
      !  compute new water level
      dy = 0d0; dd = 0d0    !initializes
      IF(     w > 0d0 ) dy = (flparm%volin-flparm%vol)/w                      !in terms of width
      IF( facts > 0d0 ) dd = (1d0-flparm%vol/flparm%volin)*flparm%kvol*dtime/gamma  !in terms of compressibility
      IF( ABS(dd) > ABS(dy) .AND. dy /= 0d0 )THEN
        y0 = y0 + dy
      ELSE
        y0 = y0 + dd
      END IF
      flparm%press = y0     !update water level
      xl = 0d0  ; xr = 0d0  !initializes coordinates for free surfaces
      full = .FALSE.        !bag may contain air
    ELSE IF(flparm%nodcen < 14 )THEN  !water volume is constant
      IF( facts > 0d0 ) dd = (1d0-flparm%vol/flparm%volin)*flparm%kvol*dtime/gamma  !in terms of compressibility
      y0 = y0 + dd
      flparm%press = y0     !update water level
    ELSE ! IF(flparm%nodcen < 24 )THEN !water volume changes with curve
      flparm%volin = flparm%volin + facts*dtime
      !dd = (1d0-flparm%vol/flparm%volin)*flparm%kvol*dtime/gamma  !in terms of compressibility
      dd = (1d0-flparm%vol/flparm%volin)            !volumetric strain (-)
      IF( facts > 0.d0 )THEN
        dd = dd*(1d0+1d2*ABS(dd))*flparm%kvol*dtime/gamma  !filling process
      ELSE
        dd = dd*flparm%kvol*dtime/gamma  !constant volume or decreasing
      END IF
      y0 = y0 + dd
      flparm%press = y0     !update water level
    END IF
    x0    = flparm%timini !x-coordinate for volume computation
    vol = 0d0             !Initializes volume
  END IF

  seg => headf
  DO
    IF (.NOT.ASSOCIATED (seg) ) EXIT

    ! computes new water pressure and volume
    IF( water )THEN                   !for water bags
      IF (ndime == 3) THEN !only for triangular element
        x(1,1:3) = coora(1,seg%lnofl(1:3)) - x0 !distance to reference line
        x(2,1:3) = y0 - coora(2,seg%lnofl(1:3)) !depth at nodes
        !*** evaluate the first two side vectors
        l1 = coora(:,seg%lnofl(1)) - coora(:,seg%lnofl(2))       !side 1
        l2 = coora(:,seg%lnofl(3)) - coora(:,seg%lnofl(1))       !side 2
        !*** evaluate the cross product => plane normal
        CALL vecpro(l1,l2,t3)                             !normal * area2  !plano yz
        CALL vecuni(3,t3,area2)                           !computes twice area
        t1 = (/-1d0,0d0,0d0/)
        area2 = area2*DOT_PRODUCT(t3,t1)
      ELSE
        x(1,1:2) = coora(1,seg%lnofl(1:2)) - x0 !distance to reference line
        x(2,1:2) = y0 - coora(2,seg%lnofl(1:2)) !depth at nodes
        IF( axisy ) rr = pi*(x(1,1)+x(1,2))     !2 Pi * r
      END IF

      IF( full )THEN
        IF(ndime == 3) THEN
          vol = vol+((x(1,1)+x(1,2)+x(1,3))/3d0)*area2/2d0
          seg%fload = (x(2,1)+x(2,2)+x(2,3))/3d0*gamma    !depth at centroid
        ELSE
          vol = vol+((x(1,1)+x(1,2))/2d0)*(x(2,2)-x(2,1))*rr
          seg%fload = (x(2,1)+x(2,2))/2d0*gamma   !depth at centroid
        END IF
      ELSE
        IF(ndime == 3) THEN   !only for triangular element
          !IF( x(2,1) >= 0d0 .AND. x(2,2)>= 0d0 .AND. x(2,3) >= 0d0 )THEN all sumerged
            vol = vol+((x(1,1)+x(1,2)+x(1,3))/3d0)*area2/2d0
            seg%fload = (x(2,1)+x(2,2)+x(2,3))/3d0*gamma   !depth at centroid
        ELSE
          IF( x(2,1) >= 0d0 .AND. x(2,2) >= 0d0 )THEN    !both sumerged
            vol = vol+((x(1,1)+x(1,2))/2d0)*(x(2,2)-x(2,1))*rr
            seg%fload = (x(2,1)+x(2,2))/2d0*gamma   !depth at centroid

          ELSE IF( x(2,1) > 0d0 .OR. x(2,2) > 0d0 )THEN !if one node is below line
            h  = x(2,1)/(x(2,1) - x(2,2))           !normalized position from node 1 (+)
            xc = (1d0-h)*x(1,1) + h*x(1,2)          !x-intersection with free surface
            IF( x(2,2) > 0 )THEN                    !node 1 is above line
              IF( axisy ) rr = pi*(xc+x(1,2))       !2 Pi * r
              vol = vol+(xc    +x(1,2))/2d0*x(2,2)*rr
              seg%fload = x(2,2)*(1d0-h)/2d0*gamma  !equivalent pressure at centroid
              IF( xc > xr ) xr = xc
            ELSE                                    !node 2 is above line
              IF( axisy ) rr = pi*(x(1,1)+xc )       !2 Pi * r
              vol = vol-(x(1,1)+xc    )/2d0*x(2,1)*rr
              seg%fload = x(2,1)*h/2d0              !equivalent pressure at centroid(-)
              IF( xc < xl ) xl = xc
            END IF
          ELSE
            seg%fload = 0d0                         !pressure at centroid
          END IF
        END IF
      END IF

    END IF      ! end water pressure computations

    press = seg%fload

    IF (ttime >= timin ) THEN
      IF (TRIM(fltype) == 'INFLAT') THEN
        press = flparm%volin*temp*(press*facts0+flparm%pext)/(vol*temp0) - flparm%pext
      ELSE IF (TRIM(fltype) == 'TUBHYD') THEN
        press = flparm%press + dp
      ELSE IF (water) THEN
        IF (flparm%nodcen < 21 ) press = press * facts
      ELSE
        press = press * facts
      END IF
    ELSE
      press = press * facts
    END IF

    nnode = seg%nnode
    x(1:ndime,1:nnode) = coora(1:ndime,seg%lnofl(1:nnode))
    IF ( nnode == 2) THEN
      CALL area2d (x(1:2,1:2), area2, vn, ntype)
      pforc = press*area2/nnode
      px(1:2) = pforc*vn(1:2)
    ELSE IF ( nnode == 3) THEN
      CALL axepld(x,trans,dummy,dummy,area2,.FALSE.)
      pforc = press*area2/2./nnode
      px = pforc*trans(1:3,3)
    ELSE IF (nnode == 4) THEN
      CALL areaq (x, area2, vn)
      pforc = press*area2/nnode
      px = pforc*vn
    END IF

    DO inode=1,nnode
      n = seg%lnofl(inode)
      resid(1:ndime,n) = resid(1:ndime,n) - px(1:ndime)
    END DO
    seg => seg%next
  END DO

  IF ( TRIM(fltype) == 'INFLAT' .OR. TRIM(fltype) == 'TUBHYD' ) THEN
    flparm%vol = vol
    flparm%press = press
  END IF
  IF ( TRIM(fltype) == 'WATERB'  ) THEN
    flparm%vol = vol      !keep present volume
    !keep width at free surface
    IF( flparm%nodcen == 1 ) THEN
      flparm%pext  = xr          !symmetric problems
    ELSE IF ( flparm%nodcen == 2 ) THEN
      flparm%pext  = xr - xl     !closed bags
    END IF
  END IF

RETURN
END SUBROUTINE aplfol

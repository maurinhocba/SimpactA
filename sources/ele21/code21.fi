      SUBROUTINE code21(nnode,nelem,ndime,dtime,lnods,matno,coord)

      !  determines the maximum time increment for 3/4 2D thermal element

      IMPLICIT NONE

      INTEGER (kind=4) :: nnode,nelem,ndime,lnods(:,:),matno(:)
      REAL (kind=8) :: dtime,coord(:,:)

      INTEGER*4 j,n,ielem,isec,osec,ln(nnode)

      REAL (kind=8) :: cp,cond,dens,deltm,xlong2,xl(nnode)
      TYPE (section), POINTER :: sec
      TYPE (mater  ), POINTER :: mat

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

        DO n=1,nnode
          ln(n)=lnods(n,ielem)
        END DO
        xl=0d0
        DO j = 1,ndime
          DO n=1,nnode-1
            xl(n)=xl(n)+(coord(j,ln(n+1))-coord(j,ln(n)))**2
          END DO
          xl(nnode)=xl(nnode)+(coord(j,ln(1))-coord(j,ln(nnode)))**2
        END DO
        xlong2=xl(1)
        DO n=2,nnode
          IF(xl(n) > 0d0 .AND. xl(n) < xlong2) xlong2=xl(n)
        END DO
        IF (cond /= 0d0) THEN
          deltm=0.5*xlong2*(dens*cp)/cond
          IF(deltm < dtime) dtime=deltm
        ENDIF
      END DO

      END SUBROUTINE code21

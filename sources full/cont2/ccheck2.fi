      SUBROUTINE ccheck2(nsegm,lcseg,nhseg,cnode,x,xc)
!
!     check surface normals  (conforming surfaces only)
!
!.... input
!....   nsegm = number of segments of the master surface considered
!....   nnseg = number of nodes per segment of the surface considered
!....   lcseg(inseg,icseg) = global node number for the local element node
!....                        [inseg] of the segment [icseg]
!....   nhseg(nnseg,nsegm) = neighbours segments to each side
!....   cnode  : external node
!....   x     = nodal coordinates
!....   xc    = coordinates of element centers

      IMPLICIT NONE
!     Dummy arguments
      INTEGER (kind=4), INTENT(IN) :: nsegm,cnode
      INTEGER (kind=4), INTENT(IN OUT) :: nhseg(:,:),lcseg(:,:)
      REAL (kind=8), INTENT(IN) :: x(:,:)
      REAL (kind=8), INTENT(OUT) :: xc(:,:)
!     local variables
      INTEGER (kind=4), PARAMETER :: jo(2) = (/ 2,1 /)
      INTEGER (kind=4) icseg,near,j,k
      REAL(kind=8) :: d(2),g,eta,t(2),side(2),lg,mind
      LOGICAL :: inrange,ok
      LOGICAL :: check(nsegm)

!     search the nearest segment

      mind = 1e20
      DO icseg = 1,nsegm                  !for each segment
        xc(:,icseg) = (x(:,lcseg(1,icseg)) + x(:,lcseg(2,icseg)) )/2d0  !center coordinates
        d = x(:,cnode) - xc(:,icseg)                                    !vector to center
        lg = DOT_PRODUCT(d,d)                                           !distance to center
        IF( lg < mind )THEN                     !check minimun distance
          mind = lg
          near = icseg
        END IF
      END DO

      !check nearest segment
      side(:) = (x(:,lcseg(2,near)) - x(:,lcseg(1,near)) )  !side vector
      CALL vecuni(2,side,lg)                                !unit vector along side
      t = (/ -side(2),side(1) /)                            !unit norma
      d = x(:,cnode) - xc(:,near)                           !vector to center
      g = DOT_PRODUCT(d,t)                                  !gap
      eta = 2d0*DOT_PRODUCT(d,side)/lg                      !local coordinates [-1:1]
      inrange = ABS(eta) < 2d0                              !extended in both directions

      check = .FALSE.
      IF( g >= 0d0 )THEN !node is outside the surface
        ok = .TRUE.
      ELSE
        IF(inrange) THEN
          ok = .FALSE.
          lcseg(:,near) = (/ lcseg(2,near),lcseg(1,near) /)
          nhseg(:,near) = (/ nhseg(2,near),nhseg(1,near) /)
        ELSE
          WRITE(*,*) 'check can not be performed for this surface'
          RETURN
        END IF
      END IF
      check(near) = .TRUE.
      ! this have sense for a simple surface (closed or open)
      icseg = near            !seed segment
      j = 1                   !begin with left neighbour
      DO
        k = nhseg(j,icseg)
        IF( k /= 0 )THEN
          IF( check(k) )EXIT               !it means that all the segments have been checked
          IF( lcseg(j,icseg) /= lcseg(jo(j),k) )THEN      !check continuity
            IF( lcseg(j,icseg) /= lcseg(j,k) )THEN          !for conforming surfaces
              lcseg(:,k) = (/ lcseg(2,k),lcseg(1,k) /)        !inverse connectivities
              nhseg(:,k) = (/ nhseg(2,k),nhseg(1,k) /)        !inverse neighbour
            ELSE  !non-conforming surface 
              ! compare normal with previos element
            END IF
          END IF
          check(k) = .TRUE.
          icseg = k
        ELSE   !go in the other sense
          j = jo(j)          !change to right neighbour
          icseg = near       !back to seed element
        END IF
      END DO

      RETURN
      END SUBROUTINE ccheck2

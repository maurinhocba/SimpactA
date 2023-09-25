      SUBROUTINE mastd2 (ns,ni,x,nhseg)
!
!     Compares sides to detect internal interfaces. If matched sides
!     found assign neighbour segments for contact search
!
!     INPUT
!       ns    : number of boundary sides
!       ni(4,ns): boundary sides data
!         (1) : associated segment
!         (2) : associated side
!         (3) : node of the side
!         (4) : opposite node of the side
!       x(2,:): nodal coordinates
!
!    OUTPUT
!       nhseg(2,:) neighbours segments (connected segments to each side)

      IMPLICIT NONE

!     Dummy arguments
      INTEGER (kind=4), INTENT(IN) :: ns,ni(:,:)
      INTEGER (kind=4), INTENT(OUT) :: nhseg(:,:)
      REAL (kind=8), INTENT (IN) :: x(:,:)

!     local variables
      INTEGER (kind=4) i,j,ii,jj
      REAL (kind=8) :: v1(2),di,d(2)
      REAL (kind=8), ALLOCATABLE :: l(:),y(:,:)

!     Generate tangential vectors and lengths of each side (first-second)

      ALLOCATE ( l(ns),   &      !length of segments
                 y(2,ns) )       !coordinates of nodes

      DO i=1,ns                              !for each side
        ii = ni(3,i)                         !first node of the side
        jj = ni(4,i)                         !second node of the side
        y(1:2,i) = x(1:2,ii)                 !coordinates of the node
        d = x(1:2,jj) - x(1:2,ii)            !oriented vector
        CALL vecuni(3,d,l(i))                !unit vector and length ==> l(i)
      END DO

      !     Search for a neighbour almost coincident node

      DO i=1,ns                              !for each boundary node
        DO j=i+1,ns                            !check for every other side
          ! .. first see if they are neighbours
          v1 = y(1:2,j) - y(1:2,i)           !vector I->J between centers
          di = SQRT(DOT_PRODUCT(v1,v1))      !distance
          IF( di > l(i)/10d0 .OR. di > l(j)/10d0 ) CYCLE  !they are far away CYCLE
          nhseg(ni(2,i),ni(1,i)) = ni(1,j)     ! neighbour found
          nhseg(ni(2,j),ni(1,j)) = ni(1,i)     ! neighbour found
        END DO
      END DO
      DEALLOCATE ( l,y )

      RETURN
      END SUBROUTINE mastd2

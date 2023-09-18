 SUBROUTINE mastda (ns,ni,x,nhseg)
!
!     Compares sides to detect internal interfaces. If matched sides
!     found assign neighbour segments for contact search
!
!     INPUT
!       ns    : number of boundary sides
!       ni(4,ns): boundary sides data
!         (1) : associated segment
!         (2) : associated side
!         (3) : first node of the side
!         (4) : second node of the side
!       x(3,:): nodal coordinates
!
!    OUTPUT
!       nhseg(3,:) neighbours segments (connected segments to each side)

 IMPLICIT NONE

!     Dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ns,ni(:,:)
 INTEGER (kind=4), INTENT(OUT) :: nhseg(:,:)
 REAL (kind=8), POINTER  :: x(:,:)

!     local variables
 REAL (kind=8), PARAMETER :: tolc=0.94d0, told=0.1d0
 INTEGER (kind=4) i,ii,j,jj,k,kk
 REAL (kind=8) pro,proi,cosa,di,v1(3),d(3)
 REAL (kind=8), ALLOCATABLE :: t(:,:),l(:),y(:,:),r(:)

!     Generate tangential vectors and lengths of each side (first-second)

 ALLOCATE ( t(3,ns), &      !unit vector of sides
            l(ns),   &      !length of sides
            y(3,ns), &      !coordinates of center of side
            r(ns) )         !best projection coordinate

 DO i=1,ns                              !for each side
   ii = ni(3,i)                         !first node of the side
   jj = ni(4,i)                         !second node of the side
   y(1:3,i) = (x(1:3,ii)+x(1:3,jj))/2d0 !coordinates of the center
   d = x(1:3,jj) - x(1:3,ii)            !oriented vector
   CALL vecuni(3,d,l(i))                !unit vector and length ==> l(i)
   t(1:3,i) = d                         !assign unit vector ==> t(i)
   r(i) = 1d1                           !anything large
 END DO

 !     Search for a neighbour "parallel" side

 DO i=1,ns                              !for each boundary side
   ii = ni(1,i)                         !associated segment to side I
   k  = ni(2,i)                         !associated side to side I
   DO j=i+1,ns                            !check for every other side
     ! .. first see if they are neighbours
     v1 = y(1:3,j) - y(1:3,i)           !vector I->J between centers
     di = SQRT(DOT_PRODUCT(v1,v1))      !distance
     IF( di > (l(i)+l(j))/2d0 ) CYCLE   !they are far away CYCLE
     ! .. second check if segments are parallel
     cosa = -DOT_PRODUCT(t(1:3,i),t(1:3,j))   !angle cosine (must be neg)
     IF( cosa < tolc ) CYCLE            !if angle too large CYCLE
     !check if center of I side projects over side J
     jj = ni(1,j)                       !associated segment to side J
     kk = ni(2,j)                       !associated side to side J
     !next 2 lines, t(i) is choosen but t(i) and t(j) are almost parallel
     pro = DOT_PRODUCT(t(1:3,i),v1)    !proyection vector I-J over side
     d  = v1 - pro*t(1:3,i)            !Normal vector from side I to side J
     CALL vecuni(3,d,di)               !distance between sides
     IF( di < told*MAX( l(i),l(j)) )THEN !if in tolerance
       ! check side I
       proi = ABS(pro/l(i))     !normalized local coordinate
       IF( proi < r(i))THEN
         nhseg(k,ii) = jj    !'better' neighbour found
         r(i) = proi          !keep local coordinate
       END IF
       ! check side J
       proi = ABS(pro/l(j))     !normalized local coordinate
       IF( proi < r(j))THEN
         nhseg(kk,jj) = ii    !'better' neighbour found
         r(j) = proi          !keep local coordinate
       END IF
     END IF
   END DO
 END DO
!      WRITE(55,"(10f8.5)")r
 DEALLOCATE ( t,l,y,r )

 RETURN
 END SUBROUTINE mastda

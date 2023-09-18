 SUBROUTINE gaus15(heade,coord,iffix,istop,gauss, &
                   angdf,nrf,rhead,nelem,         &
                   shear,shears,factors,ninv)
 USE ctrl_db, ONLY : npoin
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for shell element 15
 !--------------------------------------------------------------------
 ! for the new version:
 ! 1-compute local systems
 ! 2-compute bending stiffness normal to each side
 ! 3-compute F matrix for branching sides
 ! 4-compute initial angles
 ! 5-compute side factors using side elements data
 !

 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: iffix(:),nrf,nelem
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf !,      & !default angle X_1 - Ort_1
                                    !maxcn         !maximum normal curvature
 TYPE (ele15), POINTER :: heade   !pointer to first element in the set
 TYPE (reinf), POINTER :: rhead   !pointer to first branching side
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists
 INTEGER(kind=4), INTENT (IN) :: shear      !/= 0 if shear are to be computed
 REAL(kind=8), POINTER :: shears(:,:),factors(:)
 INTEGER(kind=4), POINTER :: ninv(:)

 INTEGER (kind=4), PARAMETER :: nextn(3) = (/2,3,1/) !cycling list
 INTEGER (kind=4) i,j,k,l,m,jj,kk,osec,isec,ni,sides(3),ielem,mtype
 REAL (kind=8) x(3,6),a2(4),stra0(6),t(3,15),c(4),y(3)
 REAL (kind=8) angle,factor,ri,fi,lh,ls,a,b,n1,n2,t3,nu
 TYPE (ele15), POINTER :: e,ej
 TYPE (reinf), POINTER :: si !pointer to a branching side
 LOGICAL :: is(3),newmt

 TYPE (section), POINTER  :: sec                  !section pointers
 TYPE (mater), POINTER  :: mat                    !material pointers

 !     ---------------------------------------------------------
 factor = ASIN(1d0)/90d0         !pi/180

 IF( .NOT.gauss )THEN
   !  associate other element to clamped stiffeners
   si => rhead            !point to first side
   ! search each side
   DO i=1,nrf             !for each side
     l = si%lnods(3)    ! first side element
     m = si%lnods(4)    ! second side element
     IF( l == 0 .OR. m == 0 )THEN
       j = si%lnods(1)      !internal node numeration (first node of the side)
       k = si%lnods(2)      !internal node numeration (second node)
       y = (coord(:,j)+coord(:,k))/2d0 !reinforcement center
       j = 3                !left element missing
       IF( m == 0 ) j = 4   !rigth element missing
       l = MAX(l,m)     !Element
       m = MOD(l,10)    !associated side (J-K)
       l = l/10         !first side element number
       e => heade       !point to first
       DO jj=1,l-1         !advance until element found
         e => e%next
       END DO
       lh = HUGE(1d0)        !a large number
       DO kk=1,2               !for each side
         m = nextn(m)        !move ahead on the side
         IF( e%lside(m) > 0 )THEN  !if side element exist
           x(:,1:3) = coord(:,e%lnods(hh(:,m)))   !element coordinates
           DO jj=1,3               !compute distance to center
             x(jj,4) = y(jj) - SUM(x(jj,1:3))/3d0
           END DO
           ls = DOT_PRODUCT(x(:,4),x(:,4))
           IF( ls < lh )THEN
             l = m           !side candidate
             lh = ls
           END IF
         END IF
       END DO
       l = e%lside(l)     !nearest side element
       ej => heade        !point to first
       DO jj=1,l-1
         ej => ej%next          !advance until element found
       END DO
       DO kk=1,3               !find first free position
         IF( .NOT.ASSOCIATED(ej%si(kk)%p) )THEN
           ej%si(kk)%p => si   !point to reinforcement
           si%lnods(j) = -l*10-kk    !keep element and side (- is a flag)
           EXIT
         END IF
       END DO
     END IF
     si => si%next
   END DO
   !     compute local systems and bending stiffness at each side
   e => heade       !point to first element
   DO ielem=1,nelem           !loop over the elements in the set
     isec = e%matno                   !associated material
     newmt = osec /= isec             !same material than previous ?
     IF(newmt)THEN                    !if new material, get properties
       sec => psecs(isec)%p
       mat => sec%mtbas
       osec = isec                       ! keep present material
       mtype = mat%mtype
       t3    = sec%rprop(1)**3
       IF( mtype == 1)THEN
         ri = mat%prope(7)*t3              !C1 * thickness^3
         nu = mat%prope(2) ! mat%prope(8)/mat%prope(7)    !Poisson ratio
       ELSE IF( mtype == 5)THEN            !elastic orthotropic
         ! plane stress elasticity matrix
         c(1:4) = mat%prope(16:19)*t3
       !ELSE IF( mtype == 10 )THEN               ! Elastic Rubber
       !  nto = mat%matdef(6)                    !number of terms in Ogden model
       !  eprop = RESHAPE(mat%prope(7:12),(/2,3/))  ! elastic properties
       END IF
     END IF

     x(1:3,1:3) = coord(1:3,e%lnods(1:3))    !get original coordinates of main triangle
     DO m =1,3                       !for the side nodes in the patch
       i = m+3                       !associated node
       ni = e%lnods(i)               !node number
       SELECT CASE (ni)
       CASE (1:)                   !one adjacent element
         x(1:3,i) = coord(1:3,ni)    !get original coordinates of adjacent node
         sides(m) = 1                !side element exists
       CASE (0)                    !boundary node
         sides(m) = 0                !side element does not exists
         j = e%lnods(hh(3,m))        !first node of the side
         k = e%lnods(hh(2,m))        !second node of the side
         is(m) = ( iffix(j) == 1 .AND. iffix(k) == 1 )  !clamped side ?
         a2(m+1)  = 1d0              !initializes side area (unnecessary perhaps)
         x(1:3,i) = 0d0              !null values (unnecessary perhaps)
       END SELECT
     END DO

     ! compute element local systems, side proyections and element areas
     angle = ( angdf + e%angle )*factor
     CALL axep15(x,e%a,e%b,e%c,e%cd,a2,sides,angle,e%a0)
     ! Check area
     e%area1 = a2(1)/2.d0              !store area of main triangle
     IF( e%area1 <= 0d0 ) THEN         !check that area is possitive
       WRITE(*,1000)e%numel,e%area1
       istop = 1
     END IF
     !  compute side bending stiffness
     IF( mtype == 1 .OR. mtype == 10 ) fi = ri/a2(1)   !element stiffness ri/(ls*hi)
     DO i=1,3                          !for each side
       IF(mtype == 5 ) THEN !non-isotropic material
         n1 =  e%b(i)                  !side normal
         n2 = -e%a(i)
         ls = n1**2+n2**2        !side length
         n1 = n1*n1/ls           !square cosine
         n2 = n2*n2/ls           !square sine
         ! bending stiffness in the normal direction
         fi = (c(1)*n1*n1 + 2d0*(c(2)+2d0*c(4))*n1*n2 + c(3)*n2*n2 )/a2(1) !
         IF( sides(i) == 0 .AND. .NOT.is(i) )THEN    !simple supported side
           ! bending stiffness in the normal direction
           nu = ((c(1)+c(3))*n1*n2 + c(2)*(n1*n1+n2*n2) - 4d0*n1*n2*c(4))/a2(1) !
           nu = nu/fi
         END IF
       END IF
       IF( sides(i) == 0 )THEN         !if a non existing side
         IF( is(i) )THEN          !for clamped sides
           e%ci(i) = 1.0d0        !used by clamped (and branching sides)
         ELSE
           e%ci(i) = nu           !nu for free sides
         END IF
       ELSE
         e%ci(i) = fi           !side bending stiffness ri/(ls*hi)
       END IF
       IF( ASSOCIATED(e%si(i)%p))THEN
         a = e%a(i)          !t1 proyection /2A
         b = e%b(i)          !t1 proyection /2A
         ls = SQRT(a**2+b**2)      ! factor =  (1/2A)^2
         IF( e%si(i)%p%lnods(1) == e%lnods(hh(3,i)) .AND. e%si(i)%p%lnods(3) > 0 )THEN
           e%si(i)%p%fc(1) = 1d0/e%area1  !left triangle
           e%si(i)%p%dc = (/ a/ls,b/ls /)
         ELSE IF(e%si(i)%p%lnods(1) == e%lnods(hh(2,i)) .AND. e%si(i)%p%lnods(4) > 0  )THEN
           e%si(i)%p%fc(2) = 1d0/e%area1  !rigth triangle
           e%si(i)%p%dc = (/ a/ls,b/ls /)
         ELSE !not an adjacent triangle
           IF( e%si(i)%p%lnods(3) < 0 )THEN  !if left triangle is missing
             e%si(i)%p%fc(1) = 1d0/e%area1 /1.2d0  !left triangle
           ELSE
             e%si(i)%p%fc(2) = 1d0/e%area1 /1.2d0  !rigth triangle
           END IF
         END IF
       END IF
     END DO

     e => e%next          !point to next element
   END DO

   !  compute length of reinforced sides
   si => rhead   !point to first side
   DO i=1,nrf    !for each branching side
     x(:,1:2) = coord(:,si%lnods(1:2))  !side coordinates
     t(:,9)   = x(:,2) - x(:,1)         !side vector j ==> k  (s)
     CALL vecuni(3,t(1,9),ls)           !unit side vector
     si%l0 = ls                         !initial length (inverse)
     a = SUM(si%fc)                     !sum of inverse areas
     si%fc = si%fc /a                   !contribution factors
     si => si%next      !point to next side
   END DO

   e => heade       !point to first element
   DO ielem=1,nelem           !loop over the elements in the set
     x(1:3,1:3) = coord(1:3,e%lnods(1:3))    !get original coordinates of main triangle
     DO m =1,3                       !for the side nodes in the patch
       i = m+3                       !associated node
       ni = e%lnods(i)               !node number
       SELECT CASE (ni)
       CASE (1:)                   !one adjacent element
         x(1:3,i) = coord(1:3,ni)    !get original coordinates of adjacent node
         sides(m) = 1                !side element exists
         ! locate element data of the neighbour element
         j = e%lside(m)              !element order
         IF(j < ielem ) CYCLE        !element already processed
         ej => e                     !point to present element
         DO k=ielem+1,j              !loop to reach j element
           ej => ej%next
         END DO
         k = 1                       !find the common side
         DO
           IF(ej%lside(k) == ielem )EXIT  !side found
           k = k+1
         END DO
         !  compute side factors
         e%ci(m)  = ej%ci(k)/(e%ci(m)+ej%ci(k))     !side factor    rj/(ri+rj)
         ej%ci(k) = 1d0 - e%ci(m)                   !side factor    ri/(ri+rj)
         e%mems(m)   = e%matno == ej%matno   !use quadratic membrane approach for this side
         ej%mems(k) = e%mems(m)
       CASE (0)                    !boundary node
         sides(m) = 0                !side element does not exists
         j = e%lnods(hh(3,m))        !first node of the side
         k = e%lnods(hh(2,m))        !second node of the side
         is(m) = ( iffix(j) == 1 .AND. iffix(k) == 1 )  !clamped side ?
         x(1:3,i) = 0d0              !null values (unnecessary perhaps)
       END SELECT
     END DO
     ! compute second fundamental form of the undeformed structure
     e%gamma = 0d0       !initializes distorsion angles
     CALL stra15(e%a,e%b,e%c,e%cd,e%ci,e%a0,e%area1,t, &
                 e%lb,x,sides,is,e%mems,.TRUE.,stra0,e%gamma)
     e%stra1 = stra0     !store actual metric tensor

     ! compute initial curvatures of reinforcement (if exists)
     DO i=1,3            !for each element side
       IF( ASSOCIATED(e%si(i)%p))THEN    !if side is reinforced
         l = 1                                    !left triangle
         IF(ABS(e%si(i)%p%lnods(4))/10 == ielem) l = 2
         !compute curvature and twist (initial, must be substracted from total curvatures)
         a = e%si(i)%p%dc(1) !t1 proyection
         b = e%si(i)%p%dc(2) !t2 proyection
         n1 = stra0(4)*a + stra0(6)*b   !twice curvature vector /2A
         n2 = stra0(6)*a + stra0(5)*b   !
         e%si(i)%p%curv(1,l) = (a*n1+b*n2)/2d0      !bending curvature
         e%si(i)%p%curv(2,l) = (b*n1-a*n2)/2d0      !torsion curvature
       END IF
     END DO

     e => e%next          !point to next element
   END DO
   gauss = .TRUE.
   !     ---------------------------------------------------------
   IF( shear < 0 )THEN  !initializes values for shear computation
     ALLOCATE( ninv(npoin) )
     ninv = 0
     e => heade       !point to first element
     DO ielem=1,nelem           !loop over the elements in the set
       ninv(e%lnods(1:3)) = 1
       e => e%next
     END DO
     k = 0
     DO i=1,npoin
       IF( ninv(i) == 0 )CYCLE
       k = k + 1
       ninv(i) = k
     END DO
     ALLOCATE( shears(2,k),factors(k) )
     shears = 0d0
     factors = 0d0
     e => heade       !point to first element
     DO ielem=1,nelem           !loop over the elements in the set
       DO i=1,3
         j = ninv(e%lnods(i))
         factors(j) = factors(j) + 1d0/e%area1
       END DO
       e => e%next
     END DO
     !  factors
     DO i=1,k
       factors(i) = 1d0/factors(i)
     END DO
   END IF
   !     ---------------------------------------------------------
 END IF
 1000 FORMAT(' NEGATIVE OR NUL JACOBIAN AT SHELL ELEMENT (15) ',I5,/ &
      &      ' JACOBIAN = ', E12.4 )

 RETURN
 END SUBROUTINE gaus15

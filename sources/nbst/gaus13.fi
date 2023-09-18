 SUBROUTINE gaus13(heade,coord,iffix,istop,gauss, &
                   angdf,nbs,bhead,nelem,shearc,        &
                   moments,factors,ninv,locax)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for shell element 13
 !--------------------------------------------------------------------
 ! for the new version:
 ! 1-compute local systems
 ! 2-compute bending stiffness normal to each side
 ! 3-compute F matrix for branching sides
 ! 4-compute initial angles
 ! 5-compute side factors using side elements data
 !
 ! It is not clear when actual configuration is used!
 USE npo_db, ONLY : label

 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: iffix(:),nbs,nelem
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf !,      & !default angle X_1 - Ort_1
                                    !maxcn         !maximum normal curvature
 TYPE (ele13), POINTER :: heade   !pointer to first element in the set
 TYPE (sideb), POINTER :: bhead   !pointer to first branching side
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists
 INTEGER(kind=4), INTENT (IN) :: shearc     !/= 0 if shear are to be computed
 INTEGER(kind=4), INTENT (IN) :: locax      ! local axis
  REAL(kind=8), POINTER :: moments(:,:),factors(:)
 INTEGER(kind=4), POINTER :: ninv(:)

 INTEGER (kind=4) i,j,k,m,nn,mm,osec,isec,ni,sides(3),ielem,secty
 REAL (kind=8) x(3,6),a2(4),stra0(6),t(3,15),n(3,8),c(6)
 REAL (kind=8) angle,factor,fi,ls,lh,a,b,cosa,sina,n1,n2,n12,t3,nu,thick
 TYPE (ele13), POINTER :: e,ej
 TYPE (sideb), POINTER :: si !pointer to a branching side
 LOGICAL :: is(3),newmt

 TYPE (section), POINTER  :: sec                  !section pointers
 TYPE (mater), POINTER  :: mat                    !material pointers

 !     ---------------------------------------------------------
 factor = ASIN(1d0)/90d0         !pi/180

 IF( .NOT.gauss )THEN
   !     compute local systems and bending stiffness at each side
   osec = -1        !initializes
   e => heade       !point to first element
   DO ielem=1,nelem           !loop over the elements in the set
     isec = e%matno                   !associated material
     newmt = osec /= isec             !same material than previous ?
     IF(newmt)THEN                    !if new material, get properties
       sec => psecs(isec)%p
       mat => sec%mtbas
       osec = isec                       ! keep present material
       secty = sec%secty
       thick = sec%rprop(1)
       t3    = thick**3
       IF( sec%secty == 12 )THEN
         c(1:4) = sec%rprop(9:12)             !D11,D12,D22,DS
       ELSE !secty = 13
         c(1:6) = sec%rprop(12:17)            !D11,D12,D13,22,D23,DS
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
         j = e%lnods(hh(3,m))        !first node of the side
         k = e%lnods(hh(2,m))        !second node of the side
         IF( iffix(j) < 0 .AND. iffix(k) < 0 )THEN !node on a side
           mm = -iffix(j)  !compute position of the i+1 node
           t(:,1) =  coord(:,nndpd(2,mm))  !node 1
           t(:,2) =  coord(:,nndpd(3,mm))  !node 2
           nn = -iffix(k)  !compute position of the i+2 node
           t(:,3) =  coord(:,nndpd(3,nn))  !node 3
           t(:,4) =  coord(:,nndpd(2,nn))  !node 4 of the quadrilateral
           t(:,5) = -t(:,1)+t(:,2)+t(:,3)-t(:,4) !phi,xita at the center
           t(:,10) = x(:,hh(2,m))-x(:,hh(3,m))    !side vector
           CALL vecpro(t(1,5),t(1,10),x(1,i)) !outward (to the solid) normal vector
           CALL vecuni(3,x(1,i),ls)          !unit normal to the quad
           CALL vecpro(t(1,10),x(1,i),t(1,9)) !normal vector to the element at the side
           CALL vecuni(3,t(1,9),ls)          !unit normal to the element
           t(:,9) = t(:,9)*thick/2d0         !t * h/2
           t(:,5) = x(:,hh(3,m)) - t(:,9)    !
           t(:,6) = x(:,hh(3,m)) + t(:,9)
           t(:,7) = x(:,hh(2,m)) + t(:,9)
           t(:,8) = x(:,hh(2,m)) - t(:,9)
           e%cd(:,:,m) = RESHAPE( (/ -1, 1, 1, -1, -1, -1, 1, 1 /), (/ 4,2 /) ) !initial guess
           DO i=1,4
             CALL projq4(t(1,i+4),t(1,1),e%cd(i,1,m),e%cd(i,2,m),a2(1),.FALSE.)
           END DO
           sides(m) = 2              !side element is a solid
           ! check orientation
           t(:,11) = x(:,m)-x(:,hh(3,m))    !side vector j-->i
           CALL vecpro(t(1,10),t(1,11),t(1,12)) !outward (to the solid) normal vector
           CALL vecuni(3,t(1,12),ls)           !unit normal to the triangle
           ls = DOT_PRODUCT(t(:,9),t(:,12))
           IF( ls < 0d0 )THEN
             WRITE(55,9001)e%numel,label(e%lnods(1:3)),  &
                           label(j),label(nndpd(2,mm)) ,label(nndpd(2,mm)), &
                           label(k),label(nndpd(2,nn)) ,label(nndpd(2,nn))
             WRITE(* ,9001)e%numel,label(e%lnods(1:3)),  &
                           label(j),label(nndpd(2,mm)) ,label(nndpd(2,mm)), &
                           label(k),label(nndpd(2,nn)) ,label(nndpd(2,nn))
             WRITE(lures,9001)e%numel,label(e%lnods(1:3)),  &
                              label(j),label(nndpd(2,mm)) ,label(nndpd(2,mm)), &
                              label(k),label(nndpd(2,nn)) ,label(nndpd(2,nn))

           END IF
         ELSE
           x(1:3,i) = 0d0            !null values (unnecessary perhaps)
           e%c(:,m,2)  = 0d0
           e%a0(m)     = 0d0
           e%cd(:,:,m) = 0d0
           sides(m) = 0              !side element is a solid
         END IF
         is(m) = ( iffix(j) /= 0 .AND. iffix(k) /= 0 )  !clamped side ?
         !a2(m+1)  = 1d0              !initializes side area (unnecessary perhaps)
       CASE (:-1)                  !branching side
         sides(m) = -1               !side element is a branching side
       END SELECT
     END DO

     ! compute element local systems, side proyections and element areas
     angle = ( angdf + e%angle )*factor
     CALL axep13(x,e%a,e%b,e%c,e%cd,a2,sides,angle,e%a0,locax)
     ! Check area
     e%area1 = a2(1)/2.d0              !store area of main triangle
     IF( e%area1 <= 0d0 ) THEN         !check that area is possitive
       WRITE(*,1000)e%numel,e%area1
       istop = 1
     END IF
     !  compute side bending stiffness
     !IF( mtype == 1 .OR. mtype == 8 .OR. mtype == 10 ) fi = ri/a2(1)   !element stiffness ri/(ls*hi)
     DO i=1,3                          !for each side
       n1 =  e%b(i)                  !side normal
       n2 = -e%a(i)
       ls = n1**2+n2**2        !side length
       n12= n1*n2/ls           !cosine sine
       n1 = n1*n1/ls           !square cosine
       n2 = n2*n2/ls           !square sine
       ! bending stiffness in the normal direction
       IF( secty == 12 ) THEN
         fi = (c(1)*n1*n1 + 2d0*(c(2)+2d0*c(4))*n1*n2 + c(3)*n2*n2 )/a2(1) !
       ELSE  !secty == 13
         fi = (c(1)*n1*n1 + 2d0*(c(2)+2d0*c(6))*n1*n2 + c(4)*n2*n2 &
            - 4d0*(c(3)*n1+c(5)*n2)*n12)/a2(1) !
       END IF
       IF( sides(i) == 0 .AND. .NOT.is(i) )THEN    !simple supported side
         ! bending stiffness in the transversal direction => equivalent poisson ratio
         IF( secty == 12 ) THEN
           nu = ((c(1)+c(3))*n1*n2 + c(2)*(n1*n1+n2*n2) - 4d0*n1*n2*c(4))/a2(1) !
         ELSE  !secty == 13
           nu = ((c(1)+c(4))*n1*n2 + c(2)*(n1*n1+n2*n2) - 4d0*n1*n2*c(6)  &
              + 4d0*(c(3)-c(5))*n12*(n1-n2))/a2(1) !
         END IF
         nu = nu/fi
       END IF

       IF( sides(i) == 0 .OR. sides(i) == 2)THEN         !if a non existing side
         IF( is(i) )THEN          !for clamped sides
           e%ci(i) = 1.0d0        !used by clamped (and branching sides)
         ELSE
           e%ci(i) = nu           !nu for free sides
         END IF
       ELSE
         IF( sides(i) < 0 )THEN
           si => e%si(i)%p          !point SI according to side
           si%fc(-e%lnods(i+3),si%nn) = fi
           e%ci(i) = 1d0
         ELSE
           e%ci(i) = fi           !side bending stiffness ri/(ls*hi)
         END IF
       END IF
     END DO

     e => e%next          !point to next element
   END DO

   !  compute constants for branching nodes (original angles, normal derivatives et al)
   si => bhead   !point to first side
   DO i=1,nbs    !for each branching side
     nn = si%nn                         !number of connected elements
     x(:,1:2) = coord(:,si%lnods(-1:0)) !branching side coordinates (ordered)
     t(:,9)   = x(:,2) - x(:,1)         !side vector j ==> k  (s)
     CALL vecuni(3,t(1,9),ls)           !unit side vector
     si%alph0(nn) = 1d0/ls              !initial length (inverse)
     a  = 0d0                           !initializes stiffness sum
     DO j=1,nn                          !for each connected element
       ni = ABS(si%lnods(j))            !node number
       x(:,3) = coord(:,ni)             !coordinates of the node
       n(:,j) = x(:,3) - x(:,1)         !side (-k) vector j ==> i
       b = DOT_PRODUCT(n(:,j),t(:,9))   !side j-i (-k) proyection over S (note sign)
       CALL vecpro(t(1,9),n(1,j),t(1,j))!normal vector to the plane j-k x j-i = t(i)
       CALL vecuni(3,t(1,j),lh)         !unit normal vector and triangle heigth (lh)
       si%fc(j,nn) = si%fc(j,nn)*ls     !stiffness / heigth
       a  = a + si%fc(j,nn)             !stiffness sum
       CALL vecpro(t(1,9),t(1,j),n(1,j))!normal vector to the side (in-plane) s x t(i) = n(i)
       si%c(3,j) = -1d0/lh              !node (i) -b(i)/2A = 1/h     (opposite node)
       si%c(2,j) =    b/lh/ls           !node (k) -b(k)/2A =+b/h/ls
       si%c(1,j) = -si%c(2,j)-si%c(3,j) !node (j)
     END DO
     si%fc(:,nn) = si%fc(:,nn) / a      !relative stiffness factors
     a = 0d0                            !initializes factor
     DO j=1,nn-1                        !for each consecutive pair (not ordered)
       a = a + si%fc(j,nn)              !update factor a = sum(ri(1:j))
       si%fc(1  :j ,j) = 1d0-a          !column top (until diagonal)
       si%fc(j+1:nn,j) = -a             !column bot (after diagonal)
       cosa = DOT_PRODUCT(t(:,j),t(:,j+1))  !cosine alpha (elem i & elem i+1)
       sina = DOT_PRODUCT(n(:,j),t(:,j+1))  !sine alpha (elem i & elem i+1)
       si%alph0(j) = ATAN2(sina,cosa)   !original angle (elem i & elem i+1)
     END DO
     si%gamma = 0d0     !initializes distorsions (Null)
     si%bb    = 0d0     !initializes B matrix (Null)
     si => si%next      !point to next side
   END DO
   !     ---------------------------------------------------------

   !  compute plastic strains from stresses
   !IF( plstr /= 0)THEN
   !  CALL bran13(nbs,bhead,coora)
   !  e => heade       !point to first element
   !  DO ielem=1,nelem           !loop over the elements in the set
   !    CALL stst13( e, iffix, coora, .FALSE.)
   !    e => e%next
   !  END DO
   !  plstr = 0
   !END IF
   !     compute constants for elements
   e => heade       !point to first element
   DO ielem=1,nelem              !loop over the elements in the set

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
         j = e%lnods(hh(3,m))        !first node of the side
         k = e%lnods(hh(2,m))        !second node of the side
         IF( iffix(j) < 0 .AND. iffix(k) < 0 )THEN !node on a side
           nn = -iffix(j)  !compute position of the i+1 node
           t(:,1) =  coord(:,nndpd(2,nn))  !node 1
           t(:,2) =  coord(:,nndpd(3,nn))  !node 2
           nn = -iffix(k)  !compute position of the i+2 node
           t(:,3) =  coord(:,nndpd(3,nn))  !node 3
           t(:,4) =  coord(:,nndpd(2,nn))  !node 4 of the quadrilateral
           t(:,7) = -t(:,1)+t(:,2)+t(:,3)-t(:,4) !phi,xita
           t(:,8) = x(:,hh(2,m))-x(:,hh(3,m))    !side vector
           CALL vecpro(t(1,7),t(1,8),x(1,i)) !outward normal vector
           CALL vecuni(3,x(1,i),ls)          !unit normal to the quad
           e%mems(m)  = .FALSE.      ! do not use quadratic membrane approach for this side
           sides(m) = 2              ! solid side element
         ELSE
           x(1:3,i) = 0d0              !null values (unnecessary perhaps)
           sides(m) = 0                !side element does not exists
         END IF
         is(m) = ( iffix(j) /= 0 .AND. iffix(k) /= 0 )  !clamped side ?
         a2(m+1)  = 1d0              !initializes side area (unnecessary perhaps)
       CASE (:-1)                  !branching side
         sides(m) = -1               !side element is a branching side
         k = -ni                     !position in side connectivities
         si => e%si(m)%p
         ! what for if si%gamma = 0
         IF( si%lnods(k) > 0 )THEN  !according to order
           e%gamma(m) =  si%gamma(k)  !order is correct
         ELSE
           e%gamma(m) = -si%gamma(k)  !inverse order
         END IF
       END SELECT
     END DO

     ! compute second fundamental form of the undeformed structure
     e%gamma = 0d0       !initializes distorsion angles
     CALL stra13(e%a,e%b,e%c,e%cd,e%ci,e%a0,e%area1,t, &
                 e%lb,x,sides,is,e%mems,.TRUE.,stra0,e%gamma)
     e%stra1 = stra0      !store actual fundamental forms

     e => e%next          !point to next element
   END DO
   gauss = .TRUE.
   !     ---------------------------------------------------------
   IF( shearc < 0 )THEN  !initializes values for shear computation
     IF( ASSOCIATED(ninv) )DEALLOCATE(ninv)
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
     IF( ASSOCIATED(factors) )DEALLOCATE(factors,moments)
     ALLOCATE( moments(3,k),factors(k) )
     moments = 0d0
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
     ! nullify factors for nodes on branching sides
     si => bhead   !point to first side
     DO i=1,nbs    !for each branching side
       DO j=-1,0                          !for each node on the side
         m = ninv(si%lnods(j))            !node position
         factors(m) = 0d0
       END DO
       si => si%next      !point to next side
     END DO
   END IF
   !     ---------------------------------------------------------
 END IF
 1000 FORMAT(' NEGATIVE OR NUL JACOBIAN AT SHELL ELEMENT (13) ',I5,/ &
      &      ' JACOBIAN = ', E12.4 )
 9001 FORMAT(' There is an error in the definition of NARIS nodes',/ &
      &      ' Check element = ', i8, 5x, 3i7 /, &
      &      ' nodes :', i8, 5x, 2i7,/           &
      &      ' nodes :', i8, 5x, 2i7)

 RETURN
 END SUBROUTINE gaus13

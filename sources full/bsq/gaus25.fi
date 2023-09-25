 SUBROUTINE gaus25(heade,coord,iffix,istop,gauss, &
                   angdf,nbs,bhead,nelem,locax)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for shell element 25
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
 TYPE (ele25), POINTER :: heade   !pointer to first element in the set
 TYPE (sideb), POINTER :: bhead   !pointer to first branching side
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists
 INTEGER(kind=4), INTENT (IN) :: locax      ! local axis

 INTEGER (kind=4) i,j,k,l,m,nn,osec,isec,ni,sides(4),ielem,secty
 REAL (kind=8) x(3,12),a2(0:4),t0(3),tb(3,8),tns(3,2,8), &
               t(3,2,0:8),h(3,3),hb(3,2),ta(3,9),ti(3,8),n(3,8),c(6),  &
               fff(3,4),stran(9),strag(4),dt(4)
 REAL (kind=8) angle,factor,fi,ls,a,sina,cosa,n1,n2,n12,th3
 TYPE (ele25), POINTER :: e,ej
 TYPE (sideb), POINTER :: si !pointer to a branching side
 LOGICAL :: is(4),newmt

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
       th3    = sec%rprop(1)**3
       IF( sec%secty == 12 )THEN
         c(1:4) = sec%rprop(9:12)             !D11,D12,D22,DS
       ELSE !secty = 13
         c(1:6) = sec%rprop(12:17)            !D11,D12,D13,22,D23,DS
       END IF
     END IF

     x(1:3,1:4) = coord(1:3,e%lnods(1:4))    !get original coordinates of main quadrilateral
     DO m =1,4                       !for the side nodes in the patch
       i = fn(m)                     !first node of associated side element (local)
       ni = e%lnods(i)               !first node of associated side element (global)
       SELECT CASE (ni)
       CASE (1:)                   !one adjacent element
         x(1:3,i) = coord(1:3,ni)    !get first original coordinates of adjacent node
         i=i+1
         ni = e%lnods(i)             !second node of associated side element (global)
         x(1:3,i) = coord(1:3,ni)    !get second original coordinates of adjacent node
         sides(m) = 1                !side element exists
       CASE (0)                    !boundary node
         sides(m) = 0                !side element does not exists
         j = e%lnods(hh(3,m))        !first node of the side
         k = e%lnods(hh(2,m))        !second node of the side
         x(1:3,m) = coord(1:3,j)        !get first node original coordinates
         x(1:3,nextn(m)) = coord(1:3,k) !get second node original coordinates
         IF( iffix(j) < 0 .AND. iffix(k) < 0 )THEN !node on a side
           l = -iffix(j)  !compute position of the i+1 node
           t(:,1,1) =  coord(:,nndpd(2,l))  !node 1
           t(:,2,1) =  coord(:,nndpd(3,l))  !node 2
           l = -iffix(k)  !compute position of the i+2 node
           t(:,1,2) =  coord(:,nndpd(3,l))  !node 3
           t(:,2,2) =  coord(:,nndpd(2,l))  !node 4 of the quadrilateral
           t(:,1,3) = -t(:,1,1)+t(:,2,1)+t(:,1,2)-t(:,2,2) !phi,xita at the center
           t(:,2,3) = x(:,nextn(m))-x(:,m)  !side vector
           CALL vecpro(t(1,1,3),t(1,2,3),x(1,i)) !outward (to the solid) normal vector
           CALL vecuni(3,x(1,i),ls)          !unit normal to the quad
           e%normd(1:3,2,m+4) = -x(:,i)         !store unit normal to the solid side
           e%normd(4,2,m+4) = 1d0               !flag
         ELSE
           x(1:3,i:i+1) = 0d0           !null values (unnecessary)
           e%normd(4,2,m+4) = 0d0       !null values (4,2,ii) is a flag
         END IF
         is(m) = ( iffix(j) /= 0 .AND. iffix(k) /= 0 )  !clamped side ?
         !a2(m+1)  = 1d0              !initializes side area (unnecessary perhaps)
         !x(1:3,i) = 0d0              !null values (unnecessary perhaps)
       CASE (:-1)                  !branching side
         sides(m) = -1               !side element is a branching side
       END SELECT
     END DO

     ! compute element local systems, side proyections and element areas
     angle = ( angdf + e%angle )*factor
     CALL axep25(x,t,t0,tb,e%cartd,e%normd,e%ns,e%area,sides,is,angle,e%a0,locax)
     ! Check area
     a2(0) = e%area(0)                   !area of main quadrilateral
     IF( a2(0) <= 0d0 ) THEN         !check that area is possitive
       WRITE(*,1000)e%numel,e%area
       istop = 1
     END IF
     !  compute side bending stiffness
    !IF( mtype == 1 .OR. mtype == 10 ) fi = ri/a2(0)   !element stiffness ri/(ls*hi)

    DO i=1,4                          !for each side
       n1 =  e%ns(1,i)                  !side normal
       n2 =  e%ns(2,i)
       ls = n1**2+n2**2        !side length
       n12= n1*n2/ls           !cosine sine
       n1 = n1*n1/ls           !square cosine
       n2 = n2*n2/ls           !square sine
       ! bending stiffness in the normal direction
       IF( secty == 12 ) THEN
         fi = (c(1)*n1*n1 + 2d0*(c(2)+2d0*c(4))*n1*n2 + c(3)*n2*n2 )/a2(0) !
       ELSE  !secty == 13
         fi = (c(1)*n1*n1 + 2d0*(c(2)+2d0*c(6))*n1*n2 + c(4)*n2*n2 &
            - 4d0*(c(3)*n1+c(5)*n2)*n12)/a2(0) !
       END IF
       IF( sides(i) < 0 )THEN
          si => e%si(i)%p          !point SI according to side
          si%fc(-e%lnods(fn(i)),si%nn) = fi
          e%ci(i) = 1d0
          si%t0(:,-e%lnods(fn(i))) = t0
          si%tb(:,-e%lnods(fn(i))) = tb(:,i)
          si%c(1,-e%lnods(fn(i)),si%nn) = e%normd(ln(2,i),1,i)
          si%c(2,-e%lnods(fn(i)),si%nn) = e%normd(ln(1,i),1,i)
          si%c(3,-e%lnods(fn(i)),si%nn) = e%normd(nextn(i+1),1,i)
          si%c(4,-e%lnods(fn(i)),si%nn) = e%normd(nextn(i+2),1,i)
        ELSE
          e%ci(i) = fi/e%area(i)           !side bending stiffness ri/(ls*hi)
        END IF
    END DO

     e => e%next          !point to next element
   END DO

   !  compute constants for branching nodes (original angles, normal derivatives et al)
   si => bhead   !point to first side
   DO i=1,nbs    !for each branching side
     nn = si%nn                         !number of connected elements
     x(:,1:2) = coord(:,si%lnods(-1:0)) !branching side coordinates (ordered)
     ta(:,9)   = x(:,2) - x(:,1)         !side vector j ==> k  (s)
     CALL vecuni(3,ta(1,9),ls)           !unit side vector
     si%alph0(nn) = 1d0/ls              !initial length (inverse)
     a  = 0d0                           !initializes stiffness sum
     k = 0
     DO j=1,nn                          !for each connected element
       ni = si%lnods(j+k)            !extra node number
       si%fc(j,nn) = si%fc(j,nn)*ls     !stiffness / heigth
       a  = a + si%fc(j,nn)             !stiffness sum
       IF (ni > 0) THEN
         ti(:,j) =+si%tb(:,j)             !normal at the side of element
       ELSE
         ti(:,j) =-si%tb(:,j)             !normal at the side of element
       END IF
       CALL vecpro(ta(1,9),ti(1,j),n(1,j)) !unit normal side vector
       k = k +1
     END DO
     si%fc(:,nn) = si%fc(:,nn) / a      !relative stiffness factors
     a = 0d0                            !initializes factor
     DO j=1,nn-1                        !for each consecutive pair (not ordered)
       a = a + si%fc(j,nn)              !update factor a = sum(ri(1:j))
       si%fc(1  :j ,j) = 1d0-a          !column top (until diagonal)
       si%fc(j+1:nn,j) = -a             !column bot (after diagonal)
       !cosa = -DOT_PRODUCT(n(:,j),n(:,j+1))  !sine alpha (elem i & elem i+1)
       !sina = +DOT_PRODUCT(n(:,j+1),ti(:,j))  !cosine alpha (elem i & elem i+1)
       !si%alph0(j) = ATAN2(sina,cosa)   !original angle (elem i & elem i+1)

       cosa = +DOT_PRODUCT(n(:,j),n(:,j+1)) !n_i . n_i+1
       sina = +DOT_PRODUCT(n(:,j+1),ti(:,j))!n_i+1 . t_i
       si%alph0(j) = ATAN2(sina,cosa)   !original angle (elem i & elem i+1)
     END DO
     si%gamma = 0d0     !initializes distorsions (Null)
     si%bbr    = 0d0     !initializes B matrix (Null)
     si => si%next      !point to next side
   END DO
   !     ---------------------------------------------------------

   !  compute actual distorsions for branching nodes
   !     compute constants for elements
   e => heade       !point to first element
   DO ielem=1,nelem              !loop over the elements in the set

     x(1:3,1:4) = coord(1:3,e%lnods(1:4))    !get original coordinates of main quadrilateral
     DO m =1,4                       !for the side nodes in the patch
       i = fn(m)                     !first node of associated side element (local)
       ni = e%lnods(i)               !first node of associated side element (global)
       SELECT CASE (ni)
       CASE (1:)                   !one adjacent element
         x(1:3,i) = coord(1:3,ni)    !get first original coordinates of adjacent node
         i=i+1
         ni = e%lnods(i)             !second node of associated side element (global)
         x(1:3,i) = coord(1:3,ni)    !get second original coordinates of adjacent node
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
       CASE (0)                    !boundary node
         sides(m) = 0                !side element does not exists
         j = e%lnods(hh(3,m))        !first node of the side
         k = e%lnods(hh(2,m))        !second node of the side
         x(1:3,m) = coord(1:3,j)        !get first node original coordinates
         x(1:3,nextn(m)) = coord(1:3,k) !get second node original coordinates
         IF( iffix(j) < 0 .AND. iffix(k) < 0 )THEN !node on a side
           l = -iffix(j)  !compute position of the i+1 node
           t(:,1,1) =  coord(:,nndpd(2,l))  !node 1
           t(:,2,1) =  coord(:,nndpd(3,l))  !node 2
           l = -iffix(k)  !compute position of the i+2 node
           t(:,1,2) =  coord(:,nndpd(3,l))  !node 3
           t(:,2,2) =  coord(:,nndpd(2,l))  !node 4 of the quadrilateral
           t(:,1,3) = -t(:,1,1)+t(:,2,1)+t(:,1,2)-t(:,2,2) !phi,xita at the center
           t(:,2,3) = x(:,nextn(m))-x(:,m)  !side vector
           CALL vecpro(t(1,1,3),t(1,2,3),x(1,i)) !outward (to the solid) normal vector
           CALL vecuni(3,x(1,i),ls)          !unit normal to the quad
           e%normd(1:3,2,m+4) = -x(:,i)         !store unit normal to the solid side
           e%normd(4,2,m+4) = 1d0               !flag
         ELSE
           x(1:3,i:i+1) = 0d0           !null values (unnecessary)
           e%normd(4,2,m+4) = 0d0       !null values (4,2,ii) is a flag
         END IF
         is(m) = ( iffix(j) /= 0 .AND. iffix(k) /= 0 )  !clamped side ?
         !a2(m+1)  = 1d0              !initializes side area (unnecessary perhaps)
         !x(1:3,i) = 0d0              !null values (unnecessary perhaps)
       CASE (:-1)                  !branching side
         sides(m) = -1               !side element is a branching side
         k = -ni                     !position in side connectivities
         si => e%si(m)%p          !point SI according to side
       END SELECT
     END DO

     ! compute second fundamental form of the undeformed structure
     CALL stra25(e%cartd,e%normd,e%ns,t,t0,tb,tns,e%a0,e%ci,x,sides,is, &
                 h,hb,dt,e%area,e%gamma,e%lb,stran,strag,fff)
     e%stra0 = stran(4:6)          !keep original curvatures
     e%stra1 = stran(1:6)          !first and second fundamental forms
     e => e%next          !point to next element
   END DO
   gauss = .TRUE.
 END IF
 1000 FORMAT(' NEGATIVE OR NUL JACOBIAN AT SHELL ELEMENT (25) ',I5,/ &
      &      ' JACOBIAN = ', E12.4 )

 RETURN
 END SUBROUTINE gaus25

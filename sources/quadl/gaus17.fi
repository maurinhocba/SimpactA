 SUBROUTINE gaus17(heade,coord,istop,gauss,angdf,ntype,ngaus,as)
 !-------------------------------------------------------------------
 !--- evaluates gauss points constants for shell element 17
 !--------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: ntype,ngaus
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !initial coordinates
                              angdf,      & !default angle X_1 - Ort_1
                              as            !stabilization parameter (NGAUS=1)
 TYPE (ele17), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists

 !local variables
 INTEGER (kind=4) i,j,k,isec,osec,ngaud
 REAL (kind=8) :: x(2,nnode),cartd(nnode,2),shape(nnode,ngaus), &
                  deriv(nnode,2,ngaus),posgp(2),weigp(2),j1(4)

 REAL (kind=8) :: angle,factor,xita,eta,r,sa,ca,thick,young
 REAL (kind=8), PARAMETER :: twopi=6.283185307179586
 TYPE (ele17), POINTER :: e
 TYPE (section), POINTER :: sec

 !     ---------------------------------------------------------

 IF( .NOT.gauss )THEN
   !       gauss points in local coordinates and weigths

   !ngaud = NINT(SQRT(REAL(ngaus)))
   ngaud = 2
   IF(ngaus == 1)ngaud=1
   CALL gaussq (ngaud ,posgp ,weigp )

   !       gauss points shape and derivatives of nodal functions

   k = 0
   DO i=1,ngaud
     xita = posgp(i)
     DO j=1,ngaud
       k = k+1
       eta = posgp(j)
       CALL shape3(deriv(1,1,k),shape(1,k),xita,eta,nnode)
     END DO
   END DO
 END IF

 !     ---------------------------------------------------------

 factor = twopi/360d0         !pi/180

 osec = -1
 e => heade       !point to first element
 DO               !loop over the elements in the set
   IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit

   IF(.NOT.gauss)THEN    !Compute initial constants if requested
     x = coord(:,e%lnods)     !get original coordinates
     IF( ntype == 3 .AND. ANY(x(1,:) < 0d0))THEN
       WRITE(55,"('For axilsymmetric problems all X1 (r) coordinates' &
              &   ' must be possitive ',/,                            &
              &   ' ELEMENT ',i6,/,(2E15.6) )",ERR=9999) e%numel,x
       CALL runen3('GAUS17: negative coord. r detected ')
     END IF

     angle = ( angdf + e%angle )*factor
     sa    = SIN(angle)
     ca    = COS(angle)

     isec = e%matno
     IF( isec /= osec )THEN
       sec => psecs(isec)%p
       IF( ntype == 1 ) THEN  !plane stress
         thick = sec%rprop(1)
       ELSE IF( ntype == 2 )THEN
         thick = 1d0
       END IF
       young = sec%mtbas%prope(1)
       osec = isec
     END IF
     ! compute element local systems, nodal functions and derivatives
     k=0
     DO i=1,ngaud
       DO j=1,ngaud
         k = k+1
         CALL jaco17(cartd,deriv(:,:,k),e%dvol(k),x,e%numel,nnode,istop,j1)
         e%dvol(k)= e%dvol(k)*weigp(i)*weigp(j)           !Gauss point area
         IF( ntype == 1 ) e%dvol(k)= e%dvol(k)*thick      !thickness por plane stress
         IF( ntype == 3 )THEN                             !2Pi * r for axilsymmetric problems
           r = DOT_PRODUCT( x(1,:),shape(:,k) )
           thick = r*twopi
           e%dvol(k)= e%dvol(k)*thick
         END IF
         ! compute cartesian derivatives in the skew system
         e%cartd(:,1,k) = cartd(:,1)*ca + cartd(:,2)*sa
         e%cartd(:,2,k) =-cartd(:,1)*sa + cartd(:,2)*ca
         e%stint(:,k) = 0d0    !initializes stresses
       END DO
     END DO
   END IF
   IF( ngaus == 1 )THEN  !compute Gamma vector ==> cartd(:,1,2)
     posgp= x(:,1)-x(:,2)+x(:,3)-x(:,4)
     xita =  posgp(1)*ca + posgp(2)*sa
     eta  = -posgp(1)*sa + posgp(2)*ca
     ! gamma vector
     e%cartd(1,1,2) = (+1d0-xita*e%cartd(1,1,1)-eta*e%cartd(1,2,1))/4d0
     e%cartd(2,1,2) = (-1d0-xita*e%cartd(2,1,1)-eta*e%cartd(2,2,1))/4d0
     e%cartd(3,1,2) = (+1d0-xita*e%cartd(3,1,1)-eta*e%cartd(3,2,1))/4d0
     e%cartd(4,1,2) = (-1d0-xita*e%cartd(4,1,1)-eta*e%cartd(4,2,1))/4d0
     ! factor Cq Cos(a) Sin(a)
     e%cartd(1,2,2) = e%dvol(1)*young*as/3d0 !Cq
     e%cartd(2,2,2) = j1(1)*j1(1)+j1(2)*j1(2)
     e%cartd(3,2,2) = j1(3)*j1(3)+j1(4)*j1(4)
     e%cartd(4,2,2) = (j1(1)*j1(3)+j1(2)*j1(4))/2d0
   END IF

   e => e%next          !point to next element
 END DO
 gauss = .TRUE.

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE gaus17

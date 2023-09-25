 SUBROUTINE gaus11(nstre,nbn,ntype,ngaus,iffix,     &
                   coord,istop,shap,          &
                   heade,gauss,nhead,strai)
 !********************************************************************
 !
 !*** evaluates gauss points constants for shell/beam element
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: nstre,nbn,ntype,ngaus,iffix(:)
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: shap(2,2),  & !shape functions
                              coord(:,:)    !initial coordinates
 TYPE (ele11), POINTER :: heade   !pointer to first element in the set
 TYPE (nodei), POINTER :: nhead   !pointer to first intersection node
 LOGICAL, INTENT (IN) :: strai    !TRUE if initial curvatures are null
 LOGICAL, INTENT (IN OUT) :: gauss    !TRUE if initial constant exists

 INTEGER (kind=4) i,j,k,ns,isec,osec,nm,ln(4)
 REAL (kind=8) x(2,nnode),y0(2),yi(2),t(2),li,a,b,pff(2,2),curv(2,2),r0(2,2),ei, a0
 LOGICAL :: fix(2),siden(2)
 TYPE (ele11), POINTER :: e
 TYPE (section), POINTER :: sec
 TYPE (nodei), POINTER :: ni,nd      !pointer to an intersection node
 !REAL (kind=8),PARAMETER ::  factor = 1.74532925199433d-02  !pi/180


 IF( .NOT.gauss )THEN
   !  compute original angles for branching nodes
   ni => nhead   !point to first node
   DO i=1,nbn    !for each node
     y0 = coord(:,ni%lnods(0))          !branching node coordinates
     a  = 0d0                           !initializes
     DO j=1,ni%nn                       !for each connected node
       yi = coord(:,ni%lnods(j))        !coordinates of the node
       t  = yi - y0                     !oriented vector 0-i
       ni%alph0(j) = ATAN2(t(2),t(1))   !original angle
       ni%fc(j) = ni%fc(j)/SQRT(DOT_PRODUCT(t,t))   !inverse length
       a  = a + ni%fc(j)
     END DO
     ni%bb = 0d0        !initializes B matrix
     ni%fc = ni%fc / a  !stiffness factors
     ni => ni%next      !point to next node
   END DO

   ns = nstre/2     !number of strains or stresses
   e => heade       !point to first element

   osec = -1
   DO               !loop over the elements in the set
     IF (.NOT.ASSOCIATED (e) ) EXIT  !if all element processed, exit
     isec = e%matno
     IF( isec /= osec )THEN
       sec => psecs(isec)%p
       ei = sec%rprop(12)    ! bending stiffness
       osec = isec
     END IF

     ! gathers the initial coordinates and fixities
     fix = .FALSE.       !initializes
     siden = .FALSE.     !initializes
     ln = e%lnods(1:4)        !copy nodes
     DO j=1,nnode        !for each node and end
       k = ln(j)    !node number
       SELECT CASE ( k )
       !  first case: Branching node
       CASE (:-1)        !find pointers to intersection node
         IF( j == 1 )THEN!first node
           ni => nhead   !point to first node
           DO i=1,-k-1   !loop
             ni => ni%next
           END DO
         ELSE            !second node
           nd => nhead   !point to first node
           DO i=1,-k-1   !loop
             nd => nd%next
           END DO
         END IF
       !  second case: border node
       CASE (0)          !find fix type at end node

         IF( j == 1 )THEN        !first node
           IF(iffix(ln(2)) < 0) THEN  !node on a side
             nm = -iffix(ln(2))  !compute side normal at left node x h
             x(:,1) = (/ coord(2,nndpd(3,nm))-coord(2,nndpd(2,nm)), &
                        -coord(1,nndpd(3,nm))+coord(1,nndpd(2,nm)) /)
             siden(1) = .TRUE.
           END IF
           fix(1) = iffix(ln(2)) /= 0    !TRUE: restraint  FALSE: free
           IF( fix(1) ) e%lnods(5) = -1  !use element position as a flag
         ELSE                    !second node
           IF(iffix(ln(3)) < 0) THEN  !node on a side
             nm = -iffix(ln(3))  !compute side normal at right node x h
             x(:,4) = (/ coord(2,nndpd(3,nm))-coord(2,nndpd(2,nm)), &
                        -coord(1,nndpd(3,nm))+coord(1,nndpd(2,nm)) /)
             siden(2) = .TRUE.
           END IF
           fix(2) = iffix(ln(3)) /= 0    !TRUE: restraint  FALSE: free
           IF( fix(2) ) e%lnods(6) = -1  !use element position as a flag
         END IF
       !  third case: inner node
       CASE (1: )        !inner node
         x(:,j) = coord(:,k)    !get coordinates
       END SELECT
     END DO

     ! compute initial geometry data
     t(:) = x(:,3) - x(:,2)           !Dx
     li = SQRT(DOT_PRODUCT(t,t))      !length
     IF(li == 0d0) THEN
       WRITE(lures,"(' null Length at beam/shell element (11):'/)",ERR=9999) e%numel
       istop = 1
     END IF
     t = t/li                         !unit tangent vector
     e%l0(2) = li                     !store initial length

     !          for first segment
     IF( e%lnods(1) > 0 )THEN         !if node exists
       yi(:) = x(:,2) - x(:,1)            !Dxi left side tangent vector
       li = SQRT(DOT_PRODUCT(yi,yi))      !length
       e%l0(1) = li                       !store initial length
       IF( strai )THEN
         e%cab(1,1) = 0d0                 !a0
       ELSE
         yi = yi/li                         !unit tangent vector
         a = DOT_PRODUCT(t,yi)              !cos a0    (t*ti)
         b = -t(2)*yi(1)+t(1)*yi(2)         !sin a0    (n*ti)
         e%cab(1,1) = ATAN2(b,a)            !a0
       END IF
       e%cab(2,1) = 1d0/li                !inverse length
       e%cab(3,1) = 0d0                   !initial distorsion
       a = ei/e%l0(2)                     !R
       b = e%cab(4,1)/li                  !Ri
       e%cab(4,1) = 2d0*b/(a+b)           !Ci
     ELSE   !keep initial angle for boundary condition or branching nodes
       IF(iffix(ln(2)) < 0) THEN  !node on a side
         a = DOT_PRODUCT(t,x(:,1))          !cos a0    (t*ti)
         b = t(2)*x(1,1)-t(1)*x(2,1)        !sin a0    (n*ti)
         a0 = ATAN2(b,a)                    !a0 initial angle with  solid
         e%cab(1,1) = ATAN2(t(2),t(1))      !a0 initial angle of shell
       ELSE                           !normal boundary
         e%cab(1,1) = ATAN2(t(2),t(1))      !a0
         IF( strai .AND. ABS(t(2)) < 0.03 )e%cab(1,1)= ATAN2(0d0,t(1))
       END IF
       e%cab(2,1) = 2d0/e%l0(2)           !twice the inverse length
       e%cab(3,1) = 0d0                   !Gamma1
       e%cab(4,1) = 1d0                   !Ci
     END IF
     !             for second segment
     IF( e%lnods(4) > 0 )THEN             !if node exists
       yi(:) = x(:,4) - x(:,3)            !Dxd rigth side tangent vector
       li = SQRT(DOT_PRODUCT(yi,yi))      !length ld
       e%l0(3) = li                       !store initial length
       IF( strai )THEN
         e%cab(1,2) = 0d0                 !a0
       ELSE
         yi = yi/li                         !unit tangent vector td
         a = DOT_PRODUCT(t,yi)              !cos a0    (t*td)
         b = t(2)*yi(1)-t(1)*yi(2)          !sin a0   -(n*td)
         e%cab(1,2) = ATAN2(b,a)            !a0
       END IF
       e%cab(2,2) = 1d0/li                !length
       e%cab(3,2) = 0d0                   !initial distortion
       a = ei/e%l0(2)                     !R
       b = e%cab(4,2)/li                  !Rd
       e%cab(4,2) = 2d0*b/(a+b)           !Cd
     ELSE   !keep initial angle for boundary condition or branching nodes
       IF(iffix(ln(3)) < 0) THEN  !node on a side
         a = DOT_PRODUCT(t,x(:,4))          !cos a0    (t*ti)
         b = t(2)*x(1,4)-t(1)*x(2,4)        !sin a0    (n*ti)
         a0 = ATAN2(b,a)                    !a0 initial angle with  solid
         e%cab(1,2) = ATAN2(t(2),t(1))      !a0 initial angle of shell
       ELSE                           !normal boundary
         e%cab(1,2) = ATAN2(t(2),t(1))      !a0
         IF( strai .AND. ABS(t(2)) < 0.03 )e%cab(1,2)= ATAN2(0d0,t(1))
       END IF
       e%cab(2,2) = 2d0/e%l0(2)           !twice the inverse length
       e%cab(3,2) = 0d0                   !Gamma2
       e%cab(4,2) = 1d0                   !Cd
     END IF
     ! compute initial curvatures
     IF( ntype == 3 )THEN
       DO i=1,ngaus
         r0(1,i) = shap(1,i)*x(1,2) + shap(2,i)*x(1,3)
         r0(2,i) = 0d0
       END DO
     END IF
     ! what do we need from STRA11
     CALL stra11(ntype,ngaus,x,ni,nd,fix,ln,e%l0,e%cab,r0, &
                 e%lambd,pff,curv,t,yi,shap)
     IF(siden(1))e%cab(1,1) = a0
     IF(siden(2))e%cab(1,2) = a0
     IF( ntype == 3 )THEN
       DO i=1,ngaus      ! for each gauss point
         e%ro(1,i) =  r0(1,i)                  !original radius
         e%ro(2,i) = -r0(1,i)*curv(2,i) / 2d0  !sin a^o
       END DO
     END IF
     e => e%next
   END DO
   gauss = .TRUE.
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE gaus11

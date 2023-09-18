 SUBROUTINE axep25(x,t,t0,tb,cd,cn,ns,area,sides,is,angle,a0,locax)
 !***********************************************************************
 !
 !    this routine compute the element local axes system
 !    for the 4 node element, and for the adjacent elements
 !    (local x-axis is directed along fiber at an Angle with
 !    standard direction (intersection with X-Y plane)
 !***********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8),INTENT(IN) :: x(3,12), &  !nodal coordinates of the patch
                             angle       !angle between standard X1 and local X1
 REAL (kind=8),INTENT(OUT) :: t(3,2,0:8), & !local x-y base vectors
                              t0(3),      & !normal at the center of main element
                              tb(3,8),    & !normal at the center of adyacent element
                              cd(4,2,0:8),& !cartesyan derivatives of shape functions
                              cn(4,2,1:8),   & !normal derivatives
                              area(0:4),  & !element area
                              a0(4),      & !initial angle
                              ns(2,4)       !side normals in the local system

 LOGICAL, INTENT(IN) :: is(4)               !type of boundary
 INTEGER(kind=4), INTENT(IN) :: locax
 INTEGER(kind=4), INTENT(IN) :: sides(4)

 ! local variables
 INTEGER (kind=4) :: i,ii
 REAL (kind=8) :: l10(3),l20(3),s,nsn(3,4)
 REAL (kind=8) :: l1(3),l2(3),t1(3),t2(3),cosa,sina,lt, &
                  jac(2,2),ji(2,2),jacn(2,2),jin(2,2)

 !     Gauss points at the center and at each mid-side
 REAL (kind=8), PARAMETER :: &  !shape function derivatives
  deriv(4,2,0:4) = RESHAPE((/-0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0,  &
                             -0.50d0, 0.50d0, 0.00d0, 0.00d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0,  &
                             -0.25d0, 0.25d0, 0.25d0,-0.25d0, 0.00d0,-0.50d0, 0.50d0, 0.00d0,  &
                              0.00d0, 0.00d0, 0.50d0,-0.50d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0,  &
                             -0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.50d0, 0.00d0, 0.00d0, 0.50d0/),&
                           (/4,2,5/))

 !*** evaluate the local system at the center

 l10 = MATMUL(x(:,1:4),deriv(:,1,0))  !x_xita_0
 l20 = MATMUL(x(:,1:4),deriv(:,2,0))  !x_eta_0

 !*** evaluate the cross product => plane normal

 CALL vecpro(l10,l20,t0)     !normal direction
 CALL vecuni(3,t0,area(0))   !t3 (unit length), length = area/4

 IF (area(0) == 0.0D0)THEN   !check
   WRITE (55,*) 'ERROR: AREA == 0 '
   WRITE(55,"(3e15.5)")x(:,1:4)
   CALL runen3('BSQ: WRONG ELEMENT DEFINITION   ')
 END IF

 SELECT CASE (locax)
 CASE (1)
   lt = (t0(2)*t0(2)+t0(3)*t0(3)) !component in th Y-Z plane
   IF( lt  < 1.0d-5) THEN         !If t0 is almost orthogonal to  Y-Z plane
     t2 = (/ -t0(3), 0d0, t0(1) /) !choose t2 orthogonal to global Y direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t0,t1)
   ELSE       !         SELECT local y=t1 in the global YZ plane
     t1 = (/ 0d0, -t0(3), t0(2)  /)
     t2 = (/ lt, -t0(2)*t0(1), -t0(3)*t0(1) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (2)
   lt = (t0(3)*t0(3)+t0(1)*t0(1)) !component in th Z-X plane
   IF( lt  < 1.0d-5) THEN         !If t0 is almost orthogonal to  Z-Y plane
     t2 = (/ t0(2), -t0(1), 0d0 /) !choose t2 orthogonal to global Z direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t0,t1)
   ELSE       !         SELECT local z=t1 in the global ZX plane
     t1 = (/ t0(3), 0d0, -t0(1)  /)
     t2 = (/  -t0(1)*t0(2), lt, -t0(3)*t0(2) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (3)
   lt = (t0(1)*t0(1)+t0(2)*t0(2)) !component in th X-Y plane
   IF( lt  < 1.0d-5) THEN         !If t0 is almost orthogonal to  X-Y plane
     t2 = (/ 0d0, t0(3), -t0(2) /) !choose t2 orthogonal to global X direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t0,t1)
   ELSE       !         SELECT local x=t1 in the global xy plane
     t1 = (/ -t0(2), t0(1) , 0d0 /)
     t2 = (/ -t0(1)*t0(3), -t0(2)*t0(3), lt /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 END SELECT

 cosa = COS(angle)                            !angle to compute
 sina = SIN(angle)                            !local X1 direction

 t(1:3,1,0) = t1*cosa + t2*sina               !local X1 direction
 t(1:3,2,0) =-t1*sina + t2*cosa               !local X2 direction

 !*** find the cartesian derivative
 !          jacobian
 jac(1,1) = DOT_PRODUCT(t(:,1,0),l10)         !x1,xi
 jac(1,2) = DOT_PRODUCT(t(:,1,0),l20)         !x1,eta
 jac(2,1) = DOT_PRODUCT(t(:,2,0),l10)         !x2,xi
 jac(2,2) = DOT_PRODUCT(t(:,2,0),l20)         !x2,eta
 !   inverse jacobian
 lt = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)     !determinant = area(0)
 ji(1,1)  =  jac(2,2)/lt                      !xi,x1
 ji(1,2)  = -jac(1,2)/lt                      !xi,x2
 ji(2,1)  = -jac(2,1)/lt                      !eta,x1
 ji(2,2)  =  jac(1,1)/lt                      !eta,x2
 cd(:,:,0) = MATMUL(deriv(:,:,0),ji)          !N(j)_x1, N(J)_x2
 CALL vecuni(3,l10,lt)                        !unit vectors for future use
 CALL vecuni(3,l20,lt)
 ! compute cartesian derivatives at each mid-side point
 DO i=1,4  !for each side

   !*** evaluate the local system at the side

   l1 = MATMUL(x(:,1:4),deriv(:,1,i))  !x_xita_i
   l2 = MATMUL(x(:,1:4),deriv(:,2,i))  !x_eta_i

   !*** evaluate the cross product => plane normal

   CALL vecpro(l1,l2,tb(1,i))               !normal direction
   area(i) = SQRT(tb(1,i)*tb(1,i)+tb(2,i)*tb(2,i)+tb(3,i)*tb(3,i)) !length = area/4

   IF (area(i) == 0.0D0)THEN
     WRITE (55,*) 'ERROR: AREA == 0 '
     WRITE(55,"(3e15.5)")x(1:3,1:4)
     CALL runen3('BSQ: WRONG ELEMENT DEFINITION   ')
   ELSE
     tb(:,i) = tb(:,i)/area(i)                   !t3 (unit length)
   END IF

   ! intersection of both planes (at center & at the side)
   IF( MOD(i,2) == 0 )THEN
     t1 = l10               !sides 2 and 4
   ELSE
     t1 = l20               !sides 1 and 3
   END IF
   CALL vecpro(tb(1,i),t1,t2)               !normal to intersection on tangent plane
   cosa = DOT_PRODUCT(t1,t(:,1,0))     !angle between intersection and t1_0
   sina = DOT_PRODUCT(t1,t(:,2,0))     !angle between intersection and t2_0
   t(:,1,i) = cosa*t1 - sina*t2        !first direction on tangent plane at g-p
   t(:,2,i) = sina*t1 + cosa*t2        !second direction on tangent plane at g-p
   !     Jacobian Matrix
   jac(1,1) = DOT_PRODUCT(t(:,1,i),l1) !x1,xi
   jac(1,2) = DOT_PRODUCT(t(:,1,i),l2) !x1,eta
   jac(2,1) = DOT_PRODUCT(t(:,2,i),l1) !x2,xi
   jac(2,2) = DOT_PRODUCT(t(:,2,i),l2) !x2,eta
   !     Jacobian Inverse
   lt = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)  !determinant = area(i)
   ji(1,1) = jac(2,2)/lt         !xi,x1
   ji(1,2) =-jac(1,2)/lt         !xi,x2
   ji(2,1) =-jac(2,1)/lt         !eta,x1
   ji(2,2) = jac(1,1)/lt         !eta,x2
   !*** find the cartesian derivative
   cd(:,:,i) = MATMUL(deriv(:,:,i),ji)
   !    side normal at the local system
   t1 = x(:,ln(1,i)) - x(:,ln(2,i)) !side
   ns(1,i) = +DOT_PRODUCT(t1,t(:,2,i))  !normal x side length
   ns(2,i) = -DOT_PRODUCT(t1,t(:,1,i))
   s = SQRT(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
   t1 = t1/s
   CALL vecpro(t1(1),tb(1,i),nsn(1,i))
   !CALL vecuni()
   jacn(1,1) = DOT_PRODUCT(nsn(:,i),l1) !x1,xi
   jacn(1,2) = DOT_PRODUCT(nsn(:,i),l2) !x1,eta
   jacn(2,1) = DOT_PRODUCT(t1,l1) !x2,xi
   jacn(2,2) = DOT_PRODUCT(t1,l2) !x2,eta
   !     Jacobian Inverse
   lt = jacn(1,1)*jacn(2,2) - jacn(1,2)*jacn(2,1)  !determinant = area(i)
   jin(1,1) = jacn(2,2)/lt         !xi,x1
   jin(1,2) =-jacn(1,2)/lt         !xi,x2
   jin(2,1) =-jacn(2,1)/lt         !eta,x1
   jin(2,2) = jacn(1,1)/lt         !eta,x2
   !*** find the cartesian derivative
   cn(:,:,i) = MATMUL(deriv(:,:,i),jin)

   ii = i + 4                   !side element position in array "t"
   IF(sides(i) == 1)THEN        !adjacent element exists
     !*** evaluate the local system at the side
     l1 = MATMUL(x(:,ln(:,i)),deriv(1:4,1,1))  !x_xita  adjacent element connectivities are such that the
     l2 = MATMUL(x(:,ln(:,i)),deriv(1:4,2,1))  !x_eta   mid-side is always the first integration point in DERIV
     !*** evaluate the cross product => plane normal
     CALL vecpro(l1,l2,tb(1,ii))                               !normal direction
     lt = SQRT(tb(1,ii)*tb(1,ii)+tb(2,ii)*tb(2,ii)+tb(3,ii)*tb(3,ii)) !length = area/4
     IF (lt == 0.0D0)THEN
       WRITE (55,*) 'ERROR: AREA == 0 '
       WRITE(55,"(3e15.5)")x(1:3,ln(:,i))
       CALL runen3('BSQ: WRONG ELEMENT DEFINITION   ')
     ELSE
       tb(:,ii) = tb(:,ii)/lt                            !t3b (unit length)
     END IF
     t1 = l1                             !intersection is the side
     CALL vecuni(3,t1,lt)                !normalizes
     CALL vecpro(tb(1,ii),t1,t2)         !normal to intersection on tangent pl
     cosa = DOT_PRODUCT(t1,t(:,1,i))     !angle between intersection and t1
     sina = DOT_PRODUCT(t1,t(:,2,i))     !angle between intersection and t2
     t(:,1,ii) = cosa*t1 - sina*t2       !first direction on tangent plane at g-p
     t(:,2,ii) = sina*t1 + cosa*t2       !second direction on tangent plane at g-p
     !     Jacobian Matrix
     jac(1,1) = DOT_PRODUCT(t(:,1,ii),l1) !x1,xi
     jac(1,2) = DOT_PRODUCT(t(:,1,ii),l2) !x1,eta
     jac(2,1) = DOT_PRODUCT(t(:,2,ii),l1) !x2,xi
     jac(2,2) = DOT_PRODUCT(t(:,2,ii),l2) !x2,eta
     !     Jacobian Inverse
     lt = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)  !determinant
     ji(1,1) = jac(2,2)/lt         !xi,x1
     ji(1,2) =-jac(1,2)/lt         !xi,x2
     ji(2,1) =-jac(2,1)/lt         !eta,x1
     ji(2,2) = jac(1,1)/lt         !eta,x2
     !*** find the cartesian derivative
     cd(:,:,ii) = MATMUL(deriv(:,:,1),ji)   !Use first integration point
     !
     !CALL vecpro(t1,tb(1,ii),nsna(1,i))
     !CALL vecpro(tb(1,ii),t1,nsna(1,i))
     !jacn(1,1) = DOT_PRODUCT(nsna(:,i),l1) !x1,xi
     !jacn(1,2) = DOT_PRODUCT(nsna(:,i),l2) !x1,eta
     jacn(1,1) = DOT_PRODUCT(t2,l1) !x1,xi
     jacn(1,2) = DOT_PRODUCT(t2,l2) !x1,eta
     jacn(2,1) = DOT_PRODUCT(t1,l1) !x2,xi
     jacn(2,2) = DOT_PRODUCT(t1,l2) !x2,eta
     !     Jacobian Inverse
     lt = jacn(1,1)*jacn(2,2) - jacn(1,2)*jacn(2,1)  !determinant
     jin(1,1) = jacn(2,2)/lt         !xi,x1
     jin(1,2) =-jacn(1,2)/lt         !xi,x2
     jin(2,1) =-jacn(2,1)/lt         !eta,x1
     jin(2,2) = jacn(1,1)/lt         !eta,x2
     !*** find the cartesian derivative
     cn(:,:,ii) = MATMUL(deriv(:,:,1),jin)        !Use first integration point
     ! compute initial curvatures and angles
     cosa = +DOT_PRODUCT(nsn(:,i),t2)       !angle between normals
     sina = -DOT_PRODUCT(nsn(:,i),tb(:,ii))
     a0(i) = ATAN2(sina,cosa)              !keep initial angle 2*alpha0
   ELSE IF(sides(i) == 0) THEN       !side is at the boundary
     IF( is(i) ) THEN !symmetry or clamped
       lt = SQRT(ns(1,i)**2+ns(2,i)**2)                      !side length
       !store cartesian derivatives in the side system & initial normal
       cn(:,1,ii) = cn(:,1,i)
       IF(cn(4,2,ii)==0d0) cn(1:3,2,ii) = MATMUL(x(:,1:4),cn(:,1,ii))            !keep original side normal
       cd(4,2,ii) = lt
     ELSE   !simple supported
       cn(:,:,ii) = 0d0 !nothing
     END IF
   ELSE !IF(sides(i) < 0 ) THEN  branching side (unnecessary perhaps)
     cn(:,2,ii)  = 0d0
     a0(i)     = 0d0
     !cd(:,:,i) = 0d0
   END IF
 END DO
 lt = SUM(area(1:4))         !sum of the jacobians at sides
 area(1:4) = area(1:4) / lt  !normalizes contributions
 area(0) = lt                !element area

 RETURN
 END SUBROUTINE axep25

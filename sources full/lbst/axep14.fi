 SUBROUTINE axep14(x,a,b,area2,sides,is,angle,ci,nonrg,locax)
 !***********************************************************************
 !
 !    this routine compute the element local axes system
 !    for the 3 node element, and for the adjacent elements
 !    (local x-axis is directed along fiber at an Angle with
 !    standard direction (intersection with X-Y plane)
 !***********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8),INTENT(IN) :: x(3,6), &  !nodal coordinates
                             angle      !angle between standard X1 and local X1
 REAL (kind=8),INTENT(OUT) :: a(3,0:3),b(3,0:3),ci(3), & !shape function derivatives
                              area2                !element areas (twice)
 LOGICAL, INTENT(IN) :: sides(3), &     !side existence
                        is(3),    &     !side boundary condition
                        nonrg           !non-regular meshes
 INTEGER(kind=4), INTENT(IN) :: locax
 ! local variables
 INTEGER (kind=4) i,j,ns,ss(2)

 REAL (kind=8) l1(3),l2(3),l3(3),t1(3),t2(3),t3(3),ll1,ll2,lt, &
               cosa,sina,cosb,sinb,cosg,sing,area(3),    &
               r11,r22,r12,t11,t22,t12,no(2,2),nd(2,2)

 !*** evaluate the side vectors
 !
 l1 = x(1:3,3) - x(1:3,2)                             !side 1
 l2 = x(1:3,1) - x(1:3,3)                             !side 2
 l3 = x(1:3,2) - x(1:3,1)                             !side 3

 !*** evaluate the cross product => plane normal

 t3(1) = l1(2)*l2(3) - l1(3)*l2(2)                    !normal * area2
 t3(2) = l1(3)*l2(1) - l1(1)*l2(3)
 t3(3) = l1(1)*l2(2) - l1(2)*l2(1)

 area2 = SQRT(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3)) !computes twice area

 IF (area2 == 0.0D0)THEN
   WRITE (55,*,ERR=9999) 'ERROR: AREA2 == 0 '
   WRITE(55,"(3e15.5)",ERR=9999) x(1:3,1:3)
   CALL runen3('AXEP14: WRONG ELEMENT DEFINITION   ')
 END IF

 t3 = t3/area2                                     !t3 (unit length)
 SELECT CASE (locax)
 CASE (1)
   lt = (t3(2)*t3(2)+t3(3)*t3(3)) !component in th Y-Z plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  Y-Z plane
     t2 = (/ -t3(3), 0d0, t3(1) /) !choose t2 orthogonal to global Y direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local y=t1 in the global YZ plane
     t1 = (/ 0d0, -t3(3), t3(2)  /)
     t2 = (/ lt, -t3(2)*t3(1), -t3(3)*t3(1) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (2)
   lt = (t3(3)*t3(3)+t3(1)*t3(1)) !component in th Z-X plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  Z-Y plane
     t2 = (/ t3(2), -t3(1), 0d0 /) !choose t2 orthogonal to global Z direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local z=t1 in the global ZX plane
     t1 = (/ t3(3), 0d0, -t3(1)  /)
     t2 = (/  -t3(1)*t3(2), lt, -t3(3)*t3(2) /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 CASE (3)
   lt = (t3(1)*t3(1)+t3(2)*t3(2)) !component in th X-Y plane
   IF( lt  < 1.0d-5) THEN         !If t3 is almost orthogonal to  X-Y plane
     t2 = (/ 0d0, t3(3), -t3(2) /) !choose t2 orthogonal to global X direction
     CALL vecuni(3,t2,lt)
     CALL vecpro(t2,t3,t1)
   ELSE       !         SELECT local x=t1 in the global xy plane
     t1 = (/ -t3(2), t3(1) , 0d0 /)
     t2 = (/ -t3(1)*t3(3), -t3(2)*t3(3), lt /)
     CALL vecuni(3,t1,lt)   !     normalizes t1 & t2
     CALL vecuni(3,t2,lt)
   END IF
 END SELECT

 cosa = COS(angle)                            !angle to compute
 sina = SIN(angle)                            !local X1 direction

 t3 = t1                                !auxiliar
 t1 = t3*cosa + t2*sina                 !local X1 direction
 t2 =-t3*sina + t2*cosa                 !local X2 direction

 !*** find the local coordinates

 a(1,0) = DOT_PRODUCT(l1,t1) ! l1 . t1
 a(2,0) = DOT_PRODUCT(l2,t1) ! l2 . t1
 a(3,0) = -a(1,0)-a(2,0)
 b(1,0) = DOT_PRODUCT(l1,t2) ! l1 . t2
 b(2,0) = DOT_PRODUCT(l2,t2) ! l2 . t2
 b(3,0) = -b(1,0)-b(2,0)

 ns = 0  !initializes number of simple supported or free sides

 DO i=1,3  !for each side
   IF(sides(i)) THEN                 !adjacent elememt exists
     l1 = x(:,kk(3,i)) - x(:,kk(2,i))              !side J-K  (I)
     l2 = x(:,kk(1,i)) - x(:,kk(3,i))              !side K-I  (J)

     ll1 = SQRT(l1(1)*l1(1)+l1(2)*l1(2)+l1(3)*l1(3)) !length of side (I)
     ll2 = SQRT(l2(1)*l2(1)+l2(2)*l2(2)+l2(3)*l2(3)) !length of side (J)

     CALL vecpro(l1,l2,t3)                          !normal vector
     area(i) = SQRT(DOT_PRODUCT(t3,t3))            !element area (twice)
     IF (area(i) == 0.0D0)THEN
       WRITE (55,*,ERR=9999) 'ERROR: AREA2 == 0 '
       WRITE(55,"(3e15.5)",ERR=9999)x(:,kk(:,i))
       CALL runen3('AXEP14: WRONG ELEMENT DEFINITION   ')
     END IF

     a(1,i) = -a(i,0)                              !projec. of side (I)
     b(1,i) = -b(i,0)
     cosa = a(1,i)/ll1                              !angle of side (I)
     sina = b(1,i)/ll1

     cosb = (l1(1)*l2(1)+l1(2)*l2(2)+l1(3)*l2(3))/ll1/ll2 !angle between (I-J)
     sinb = area(i)/ll1/ll2

     cosg = cosa*cosb - sina*sinb                    !angle of side (J) g = a+b
     sing = cosa*sinb + sina*cosb

     a(2,i) = ll2*cosg                              !projec. of side (J)
     b(2,i) = ll2*sing

     IF( nonrg )THEN      !for non-regular meshes
       t11 = b(i,0)*a(2,i)       - a(i,0)*b(2,i)        ! height on adjacent element
       t22 = b(i,0)*a(kk(2,i),0) - a(i,0)*b(kk(2,i),0)  ! height of central element
       t11  = 2d0*t11/(t11+t22)                          !coefficient for curvature computation
       t22  = 2d0-t11                                    !coefficient for curvature computation
       ci(i)= t11

       a(1:2,i) = a(1:2,i)/area(i)*t22             !shape function derivatives
       b(1:2,i) = b(1:2,i)/area(i)*t22

       a(3,i) = -a(1,i)-a(2,i)                      !side (K) shape function derivatives
       b(3,i) = -b(1,i)-b(2,i)

       a(2:3,i) = a(2:3,i) + a(kk(2:3,i),0)/area2*t11 !shape function derivatives
       b(2:3,i) = b(2:3,i) + b(kk(2:3,i),0)/area2*t11

     ELSE
       a(1:2,i) = a(1:2,i)/area(i)                 !shape function derivatives
       b(1:2,i) = b(1:2,i)/area(i)

       a(3,i) = -a(1,i)-a(2,i)                      !side (K) shape function derivatives
       b(3,i) = -b(1,i)-b(2,i)

     END IF


   ELSE IF( is(i)) THEN !clamped
     lt  = 1d0/SQRT(b(i,0)**2+a(i,0)**2)           ! side length (inverse)
     !local proyections of outer normal
     cosa =  b(i,0)*lt                              ! n1
     sina = -a(i,0)*lt                              ! n2
     !keep normal vector to plane of symmetry (outer)
     IF( b(1,i) == 0d0 ) a(:,i) = cosa*t1 + sina*t2 ! n1*t1 + n2*t2
     b(:,i) = (/ lt, cosa, sina /)                  ! side length inverse, n1 & n2
   ELSE                 !simple supported
     ns = ns+1
     ss(ns) = i
   END IF
 END DO
 SELECT CASE(ns)      !if SS or Free sides exist

 CASE (1)          !one side free or simple supported
   i = ss(1)                              !side number
   lt  = SQRT(b(i,0)**2+a(i,0)**2)        !side length
   cosa =  b(i,0)/lt                      !side normal
   sina = -a(i,0)/lt
   r11 = cosa*cosa                        !auxiliar factors
   r12 = cosa*sina
   r22 = sina*sina
   !mn =  1d0-r11*r11     -r22*r11     -2*r12*r11
   !         -r11*r22  1d0-r22*r22     -2*r12*r22
   !         -r11*r12     -r22*r12     1d0-2*r12
   a(:,i) = (/ 1d0-r11*r11,    -r22*r11,    -2*r12*r11 /) !first line of mn matrix
   b(:,i) = (/ 1d0-r22*r22,    -2*r12*r22,   1d0-2*r12 /) !(2,2),(2,3) & (3,3)
   !   mn = RESHAPE ((/ a(1,i),   a(2,i),   a(3,i),    &
   !                    a(2,i),   b(1,i),   b(2,i),    &
   !                  a(3,i)/2, b(2,i)/2,   b(3,i) /), &
   !                 (/3,3/))
 CASE (2)         !Two sides free or simple supported
   DO j=1,2       ! compute the normal vectors
     i = ss(j)    ! free side
     lt  = SQRT(b(i,0)**2+a(i,0)**2)         !side length
     no(1:2,j) = (/ b(i,0), -a(i,0) /)/lt    !side normal
   END DO
   ! second compute the dual base
   lt  = no(1,1)*no(2,2) - no(1,2)*no(2,1)   !determinant
   nd(1,1) =  no(2,2)/lt                     !inverse
   nd(2,1) = -no(1,2)/lt
   nd(1,2) = -no(2,1)/lt
   nd(2,2) =  no(1,1)/lt

   r11 = nd(1,1)*nd(1,2)                     !auxiliar factors
   r22 = nd(2,1)*nd(2,2)
   r12 = nd(1,1)*nd(2,2)+nd(2,1)*nd(1,2)
   t11 = 2d0*no(1,1)*no(1,2)
   t22 = 2d0*no(2,1)*no(2,2)
   t12 = no(1,1)*no(2,2)+no(2,1)*no(1,2)
   !mn =  r11*t11, r11*t22, r11*t12,
   !      r22*t11, r22*t22, r22*t12,
   !      r12*t11, r12*t22, r12*t12
   i = ss(1)                                 !first side
   j = ss(2)                                 !second side
   a(:,i) = (/  t11*r11, t11*r22, t11*r12 /) !first line
   b(:,i) = (/  t22*r11, t22*r22, t22*r12 /) !second line
   a(:,j) = (/  t12*r11, t12*r22, t12*r12 /) !third line
   !mn = RESHAPE ((/ a(1,i),   a(2,i),   a(3,i),
   !                 b(1,i),   b(2,i),   b(3,i),
   !                 a(1,j),   a(2,j),   a(3,j) /)  (/3,3/))

 END SELECT

 a(:,0) = a(:,0)/area2                            !shape function derivatives
 b(:,0) = b(:,0)/area2
 area2 = area2/2d0                                !keep the element area

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE axep14

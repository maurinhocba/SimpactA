 SUBROUTINE axep13(x,a,b,c,cd,area2,sides,angle,a0,locax)
 !***********************************************************************
 !
 !    this routine compute the element local axes system
 !    for the 3 node element, and for the adjacent elements
 !    (local x-axis is directed along fiber at an Angle with
 !    standard direction (intersection with X-Y plane)
 !
 !***********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8),INTENT(IN) :: angle      !angle between standard X1 and local X1
 REAL (kind=8),INTENT(IN OUT) :: x(3,6) !nodal coordinates
 REAL (kind=8),INTENT(OUT) :: a(3),b(3),c(3,3,2),cd(4,2,3),area2(4),a0(3)
 INTEGER(kind=4), INTENT(IN) :: sides(3),locax

 ! local variables
 INTEGER (kind=4) ii,jj,i,j,k
 REAL (kind=8) l1(3),l2(3),l3(3),t1(3),t2(3),t3(3),tb1(3),tb2(3),tb3(3), &
               ll1,lt,cosa,sina,t(3,3),ji(2,2),jac(2,2),cosb,sinb

 REAL (kind=8), PARAMETER :: nd(6,2,3) =            &
   (/ -0.5d0, 0.5d0,-0.5d0, 0.5d0, 0.0d0, 0.0d0 ,   &
      -0.5d0,-0.5d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0 ,   &
      -0.5d0, 0.5d0, 0.5d0, 0.0d0,-0.5d0, 0.0d0 ,   &
      -1.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ,   &
      -1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 ,   &
      -0.5d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0,-0.5d0 /)

 !*** evaluate the side vectors
 !
 l1 = x(1:3,3) - x(1:3,2)                             !side 1
 l2 = x(1:3,1) - x(1:3,3)                             !side 2
 l3 = x(1:3,2) - x(1:3,1)                             !side 3

 !*** evaluate the cross product => plane normal

 t3(1) = l1(2)*l2(3) - l1(3)*l2(2)                    !normal * area2
 t3(2) = l1(3)*l2(1) - l1(1)*l2(3)
 t3(3) = l1(1)*l2(2) - l1(2)*l2(1)

 area2(1) = SQRT(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3)) !computes twice area

 IF (area2(1) == 0.0D0)THEN
   WRITE (55,*,ERR=9999) 'ERROR: AREA2 == 0 '
   WRITE(55,"(3e15.5)",ERR=9999) x(1:3,1:3)
   CALL runen3('AXEP13: WRONG ELEMENT DEFINITION   ')
 END IF

 t3 = t3/area2(1)                                     !t3 (unit length)
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

 t(1:3,1) = t1*cosa + t2*sina                 !local X1 direction
 t(1:3,2) =-t1*sina + t2*cosa                 !local X2 direction
 t(1:3,3) = t3                                !Normal direction

 !*** find the local coordinates

 a(1) = l1(1)*t(1,1)+l1(2)*t(2,1)+l1(3)*t(3,1)      ! l1 . t1
 a(2) = l2(1)*t(1,1)+l2(2)*t(2,1)+l2(3)*t(3,1)      ! l2 . t1
 a(3) = -a(1)-a(2)
 b(1) = l1(1)*t(1,2)+l1(2)*t(2,2)+l1(3)*t(3,2)      ! l1 . t2
 b(2) = l2(1)*t(1,2)+l2(2)*t(2,2)+l2(3)*t(3,2)      ! l2 . t2
 b(3) = -b(1)-b(2)

 ! first side normal derivative
 lt = SQRT(a(1)*a(1)+b(1)*b(1))
 c(1,1,1) = -lt/area2(1)
 c(2,1,1) = -DOT_PRODUCT(l1,l2)/lt/area2(1)
 c(3,1,1) = -c(1,1,1)-c(2,1,1)
 ! second side normal derivative
 lt = SQRT(a(2)*a(2)+b(2)*b(2))
 c(2,2,1) = -lt/area2(1)
 c(1,2,1) = -DOT_PRODUCT(l1,l2)/lt/area2(1)
 c(3,2,1) = -c(1,2,1)-c(2,2,1)
 ! third side normal derivative
 lt = SQRT(a(3)*a(3)+b(3)*b(3))
 c(3,3,1) = -lt/area2(1)
 c(1,3,1) = -DOT_PRODUCT(l1,l3)/lt/area2(1)
 c(2,3,1) = -c(1,3,1)-c(3,3,1)

 DO ii=1,3     !for each side
   IF(sides(ii) == 1 )THEN                   !adjacent element exist
     jj = ii+1                               !position in AREA2
     i = hh(1,ii)                            !local first node
     j = hh(2,ii)                            !local second node
     k = hh(3,ii)                            !local third node
     l1 = x(1:3,k) - x(1:3,j)                !side J-K  (I)
     l2 = x(1:3,i) - x(1:3,k)                !side K-I  (J)

     ll1 = SQRT(b(ii)*b(ii)+a(ii)*a(ii))     !length of side (I)

     t3(1) = l1(2)*l2(3) - l1(3)*l2(2)       !normal of elem (I)
     t3(2) = l1(3)*l2(1) - l1(1)*l2(3)
     t3(3) = l1(1)*l2(2) - l1(2)*l2(1)
     area2(jj) = SQRT(t3(1)*t3(1)+t3(2)*t3(2)+t3(3)*t3(3))!Area of Elem (I)
     IF (area2(jj) == 0.0D0)THEN
       WRITE (55,*,ERR=9999) 'ERROR: AREA2 == 0, side element: ',jj
       WRITE(55,"(3e15.5)",ERR=9999),x(1:3,i),x(1:3,j),x(1:3,k)
       CALL runen3('AXEP13: WRONG ELEMENT DEFINITION   ')
     END IF
     t3 = t3/area2(jj)                       !unit normal of elem (I)
     ! normal derivative at side (outward)
     c(1,ii,2) = -ll1/area2(jj)              !extra node is always -
     c(2,ii,2) = -DOT_PRODUCT(l1,l2)/ll1/area2(jj)
     c(3,ii,2) = -c(1,ii,2)-c(2,ii,2)        !sum must be 0
     ! compute initial curvatures and angles
     cosb = DOT_PRODUCT(t(:,3),t3)            !angle between normal at elements
     l3 = (b(ii)*t(:,1) - a(ii)*t(:,2))/ll1   !local normal direction (outward)
     sinb = DOT_PRODUCT(l3,t3)                !angle between normal & side n
     a0(ii) = ATAN2(sinb,cosb)                !keep initial angle
     ! compute cartesian derivatives for membrane part
     l1 = MATMUL(x,nd(:,1,ii))      ! natural derivative (xi)
     l2 = MATMUL(x,nd(:,2,ii))      ! natural derivative (eta)
     CALL vecpro(l1,l2,tb3)         ! normal at gauss point
     CALL vecuni(3,tb3,lt)          ! unit vector
     CALL vecpro(t3,tb3,tb1)        ! intersection between both planes
     CALL vecuni(3,tb1,lt)          ! unit intersection vector

     IF( lt > 1.0D-4 )THEN         !if planes are not paralell
       CALL vecpro(tb3,tb1,tb2)     !normal to intersection on tangent pl
       cosb = DOT_PRODUCT(tb1,t1)   !angle between intersection and t1
       sinb = DOT_PRODUCT(tb1,t2)   !angle between intersection and t2
       tb3 = tb1                    !use as auxiliar
       tb1 = cosb*tb3 - sinb*tb2    !first direction on tangent plane at g-p
       tb2 = sinb*tb3 + cosb*tb2    !second direction on tangent plane at g-p
     ELSE IF( lt > 1.0D-12 )THEN   !if planes are almost paralell
       CALL vecpro(tb3,t1,tb2)      ! select T2 normal to t1
       CALL vecuni(3,tb2,lt)        ! unit vector
       CALL vecpro(tb2,tb3,tb1)     ! compute selected T1
     ELSE
       tb1 = t1                 !same directions if planes are paralell
       tb2 = t2
     END IF
     IF( angle /= 0d0 )THEN
       tb3 = tb1
       tb1 = tb3*cosa + tb2*sina        !local X1 direction
       tb2 =-tb3*sina + tb2*cosa        !local X2 direction
     END IF
     !     Jacobian Matrix
     jac(1,1) = DOT_PRODUCT(tb1,l1) !x1,xi
     jac(1,2) = DOT_PRODUCT(tb1,l2) !x1,eta
     jac(2,1) = DOT_PRODUCT(tb2,l1) !x2,xi
     jac(2,2) = DOT_PRODUCT(tb2,l2) !x2,eta
     !     Jacobian Inverse
     lt = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)  !determinant
     ji(1,1) = jac(2,2)/lt         !xi,x1
     ji(1,2) =-jac(1,2)/lt         !xi,x2
     ji(2,1) =-jac(2,1)/lt         !eta,x1
     ji(2,2) = jac(1,1)/lt         !eta,x2
     !compute cartesian derivatives
     DO j=1,3             !main nodes
       cd(j,1,ii) = ji(1,1)*nd(j,1,ii) + ji(2,1)*nd(j,2,ii)   !N(j),x1
       cd(j,2,ii) = ji(1,2)*nd(j,1,ii) + ji(2,2)*nd(j,2,ii)   !N(j),x2
     END DO
     k = ii+3              !extra node
     cd(4,1,ii) = ji(1,1)*nd(k,1,ii) + ji(2,1)*nd(k,2,ii)     !N(4),x1
     cd(4,2,ii) = ji(1,2)*nd(k,1,ii) + ji(2,2)*nd(k,2,ii)     !N(4),x2
   END IF
 END DO
 ! triangle area functions derivatives
 a(1:3) = a(1:3)/area2(1)
 b(1:3) = b(1:3)/area2(1)  !note that sign is not changed here

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE axep13

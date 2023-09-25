 SUBROUTINE axep20(x,cd,area,sides,angle,ntype,r)
 !***********************************************************************
 !
 !    this routine compute the element local axes system
 !    for the 3 node element, and for the adjacent elements
 !    (local x-axis is directed along fiber at an Angle with
 !    standard direction (intersection with X-Y plane)
 !
 !    Average cartesian derivatives are also computed
 !***********************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ntype
 REAL (kind=8),INTENT(IN) :: x(2,6),  & !nodal coordinates
                             angle      !angle between standard X1 and local X1
 REAL (kind=8),INTENT(OUT) :: cd(6,3),area
 REAL (kind=8),INTENT(IN) :: r(3,2)     !radius at sides and centers
 LOGICAL, INTENT(IN) :: sides(3)

 ! local variables
 INTEGER (kind=4) i,j,k

 REAL (kind=8) l1(2),l2(2),t1(2),t2(2),t3,jac(2,2),ji(2,2),cosa,sina, &
               c(4,2),rm,f,fs

 REAL (kind=8), PARAMETER :: nd(6,2,3) =            &
   (/ -0.5d0, 0.5d0,-0.5d0, 0.5d0, 0.0d0, 0.0d0 ,   &
      -0.5d0,-0.5d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0 ,   &
      -0.5d0, 0.5d0, 0.5d0, 0.0d0,-0.5d0, 0.0d0 ,   &
      -1.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ,   &
      -1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 ,   &
      -0.5d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0,-0.5d0 /), &
               u6 = 0.1666666666666667

 ! initializes average derivatives
 cd = 0d0
 IF( ntype == 3 )THEN  !initializes shape functions
   rm = SUM(x(1,1:3))  !radius at main triangle x 3
   cd(1:3,3) = 0.5d0/rm
   !factor for side element contribution to L3
   fs = 0d0       !count the side elements
   DO i=1,3
     IF( sides (i) )fs = fs + 1d0
   END DO
   fs = 3d0/fs
 END IF

 !*** evaluate the first two side vectors

 l1 = x(1:2,3) - x(1:2,2)                             !side 1
 l2 = x(1:2,1) - x(1:2,3)                             !side 2

 !*** evaluate the cross product => plane normal

 area = (l1(1)*l2(2) - l1(2)*l2(1))/2d0  ! area

 IF (area <= 0d0)THEN
   WRITE (55,*,ERR=9999) 'ERROR: AREA2 <= 0d0 '
   WRITE(55,"(6e15.5)",ERR=9999) x(1:2,1:3)
   CALL runen3('AXEP20: WRONG ELEMENT DEFINITION   ')
 END IF

 l1 = (/  1d0, 0d0 /)
 l2 = (/  0d0, 1d0 /)

 cosa = COS(angle)                            !angle to compute
 sina = SIN(angle)                            !local X1 direction

 t1 = l1*cosa + l2*sina                 !local X1 direction
 t2 =-l1*sina + l2*cosa                 !local X2 direction

 !*** find the local cartesian derivatives at the three gauss points
 f = 0d0
 DO i=1,3       !for each Gauss point
   l1 = MATMUL(x,nd(:,1,i))      !natural derivative (xi)
   l2 = MATMUL(x,nd(:,2,i))      !natural derivative (eta)
   jac(1,1) = DOT_PRODUCT(l1,t1) !x,xi
   jac(1,2) = DOT_PRODUCT(l2,t1) !x,eta
   jac(2,1) = DOT_PRODUCT(l1,t2) !y,xi
   jac(2,2) = DOT_PRODUCT(l2,t2) !y,eta
   t3 = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)  !determinant
   ji(1,1) = jac(2,2)/t3         !xi,x
   ji(1,2) =-jac(1,2)/t3         !xi,y
   ji(2,1) =-jac(2,1)/t3         !eta,x
   ji(2,2) = jac(1,1)/t3         !eta,y
   !compute cartesian derivatives
   DO j=1,3             !main nodes
     c(j,1) = ji(1,1)*nd(j,1,i) + ji(2,1)*nd(j,2,i)
     c(j,2) = ji(1,2)*nd(j,1,i) + ji(2,2)*nd(j,2,i)
   END DO
   k = i+3              !extra node
   c(4,1) = ji(1,1)*nd(k,1,i) + ji(2,1)*nd(k,2,i)
   c(4,2) = ji(1,2)*nd(k,1,i) + ji(2,2)*nd(k,2,i)

   IF( .NOT.sides(i) )THEN    ! correct derivatives for non-existent sides
     c(i,:) = c(i,:) - c(4,:)
     c(kk(1,i),:) = c(kk(1,i),:) + c(4,:)
     c(kk(2,i),:) = c(kk(2,i),:) + c(4,:)
     c(4,:) = 0d0
   END IF
   ! pass to element arrays
   IF( ntype == 3 )THEN
     ! derivatives
     cd(1:3,1:2) = cd(1:3,1:2) + c(1:3,1:2)*r(i,1)  !main nodes
     cd(k,1:2)   = cd(k,1:2)   + c(4,1:2)  *r(i,1)  !extra node
     f = f + r(i,1)                                 !weigth factor
     IF( sides(i) ) THEN     !shape functions for hoop strain
       cd(kk(1:2,i),3) = cd(kk(1:2,i),3) + u6/r(i,2)*fs  !main nodes
       cd(k,3)         =                   u6/r(i,2)*fs  !auxiliar node
     END IF
   ELSE
     cd(1:3,1:2) = cd(1:3,1:2) + c(1:3,1:2)         !main nodes
     cd(k,1:2)   = cd(k,1:2)   + c(4,1:2)           !extra node
     f = f + 1d0
   END IF
 END DO
 cd(:,1:2) = cd(:,1:2) / f

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE axep20

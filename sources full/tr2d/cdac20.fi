 SUBROUTINE cdac20(numel,lnods,xa,cd,area1,istop,thick,angdf,angle,ntype)
 !-------------------------------------------------------------------
 !--- actualizes derivatives in eulerian formulation for element 20
 !--------------------------------------------------------------------
 IMPLICIT NONE

 !dummy variables
 INTEGER (kind=4), INTENT(IN) :: numel        !element number
 REAL(kind=8),INTENT (IN) :: xa(2,6)          !present coordinates
 REAL (kind=8), INTENT(OUT) :: cd(6,3), &     !cartessian derivatives
                               area1          !element area
 INTEGER (kind=4), INTENT(IN) :: ntype, &    !problem type
                                 lnods(6)    !patch conectivities
 INTEGER (kind=4), INTENT(OUT) :: istop     !flag for error
 REAL (kind=8), INTENT(IN) :: thick, &      !thickness for Plane Stress
                              angdf, &      !default angle X_1 - Ort_1
                              angle         !
 !local variables
 INTEGER (kind=4) i,m,n
 REAL (kind=8) r1,r(3,2),x(2,6)
 REAL (kind=8) anglf,factor
 REAL (kind=8), PARAMETER :: twopi=6.283185307209586
 LOGICAL :: sides(3)

 !     ---------------------------------------------------------
 factor = twopi/360d0         !pi/180

   x(:,1:3) = xa(:,1:3) !present coordinates of vertex nodes
   DO i =4,6                       !for all nodes in the patch
     n = lnods(i)                  !node number
     m = i-3                       !associated side element
     IF(n /= 0) THEN               !if node exist
       x(:,i) = xa(:,i)         !get original coordinates
       sides(m) = .TRUE.           !side element exists
     ELSE
       sides(m) = .FALSE.          !side element does not exists
       x(:,i) = - x(:,m) + x(:,kk(1,m)) + x(:,kk(2,m))
     END IF
   END DO

   IF( ntype == 3 )THEN
     DO i=1,3    !compute original r at each mid side
       r(i,1) = ( x(1,kk(1,i)) + x(1,kk(2,i)) )/2d0  !mid side radius
       IF( sides(i) )r(i,2) =  x(1,i+3) + r(i,1)*2d0 !side elm radius x 3
     END DO
   END IF
   ! compute element local systems, side proyections and element areas
   anglf = ( angdf + angle )*factor
   CALL axep20(x(1,1),cd(1,1),area1,sides(1),anglf,ntype,r(1,1))

   IF( area1 <= 0d0 ) THEN         !check that area is possitive
     WRITE(*,1000) numel, area1
     istop = 1
   END IF
   IF( ntype == 1 ) THEN  !plane stress
     area1= area1*thick
   ELSE IF( ntype == 3 )THEN ! include thickness or radius in the jacobian (element volume)
     r1 = SUM( x(1,1:3) )/3d0           !average radius from nodal points
     area1 = area1*r1*twopi
   END IF

 1000 FORMAT(' NEGATIVE OR NUL JACOBIAN AT 2-D TRIANGLE  (20) ',I5,/ &
&       ' AREA = ', E12.4 )

 RETURN
 END SUBROUTINE cdac20

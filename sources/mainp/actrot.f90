 SUBROUTINE actrot(lambda,theta)
!*****************************************************************
!
!     UPDATES THE NODAL LOCAL AXES
!
!*****************************************************************
 IMPLICIT NONE
!     dummy arguments
 REAL (kind=8),INTENT(IN) :: theta(3)        !rotation increment
 REAL (kind=8),INTENT(IN OUT) :: lambda(3,3) !local system
!     local variables
 REAL (kind=8) :: t(3),dl(3,3),l(3,3),modulo,c1,c2,c3

 modulo=SQRT(DOT_PRODUCT(theta,theta))      !check if Theta /= 0
 IF(modulo == 0d0) RETURN
 t = theta/modulo                           !rotation axis
 l = lambda                                 !previous local axis

 c1 = SIN(modulo)                           !factors
 c2 = 2d0*(SIN(modulo/2d0))**2
 c3 = 1d0-c2

 dl(1,1) = c3           + c2 * t(1)*t(1)    !Delta Lambda
 dl(2,1) =      c1*t(3) + c2 * t(2)*t(1)
 dl(3,1) =    - c1*t(2) + c2 * t(3)*t(1)
 dl(1,2) =    - c1*t(3) + c2 * t(1)*t(2)
 dl(2,2) = c3           + c2 * t(2)*t(2)
 dl(3,2) =      c1*t(1) + c2 * t(3)*t(2)
 dl(1,3) =      c1*t(2) + c2 * t(1)*t(3)
 dl(2,3) =    - c1*t(1) + c2 * t(2)*t(3)
 dl(3,3) = c3           + c2 * t(3)*t(3)

 lambda = MATMUL(l,dl)                      !new local axis

 RETURN
 END SUBROUTINE actrot

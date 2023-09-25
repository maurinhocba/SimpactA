 SUBROUTINE projq4(y,x,xi,eta,shape,flag)
 !
 ! proyects a point over a quadrilateral
 ! first order approach
 ! initial guess is assumed good enough
 !
 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT(IN) :: y(3), & !point coordinates
                             x(3,4)  !quad nodal coordinates
 REAL(kind=8), INTENT(IN OUT) :: xi,eta, & !guess
                                 shape(4)  !shape functions at the guess
 LOGICAL, INTENT(IN) :: flag         !to recompute shape functions
 ! local variables
 REAL(kind=8) :: deriv(4,2), &
                 xx(3,2),xy(3,2),t(3),mx,xc(3)

 ! compute shape functions and derivatives at the point
 CALL shape1 (4 ,deriv ,shape ,xi,eta )
 ! distance from point to initial guress
 xc = y - MATMUL(x,shape)
 ! jacobian matrix at the guess
 xx = MATMUL(x,deriv)
 ! compute inverse jacobian matrix (dual base)
 CALL vecpro(xx(1,1),xx(1,2),t(1))   !normal vector
 CALL vecuni(3,t,mx)                 !unit normal vector and Jacobian
 CALL vecpro(xx(1,2),t(1),xy(1,1))   !first vector of the dual base
 CALL vecpro(t(1),xx(1,1),xy(1,2))   !second vector of the dual base
 xi = xi + DOT_PRODUCT(xc,xy(:,1))/mx   !compute improved position
 eta= eta+ DOT_PRODUCT(xc,xy(:,2))/mx
 IF( flag )THEN
   mx = xi*eta
   shape(1) = (1d0-xi-eta+mx)/4d0
   shape(2) = (1d0+xi-eta-mx)/4d0
   shape(3) = (1d0+xi+eta+mx)/4d0
   shape(4) = (1d0-xi+eta-mx)/4d0
 END IF

 RETURN
 END SUBROUTINE projq4

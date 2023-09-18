 SUBROUTINE visco20(ms,m2,is,np,curve,ppro,vpar,stra,dtime,g,ier)
 !
 ! compute viscoplastic consistency parameter defined by a power law
 ! using a Newton-Raphson schemme combined with bisection algorithm
 !
 USE ctrl_db, ONLY: vefac !factor for strain rate smoothing
 USE mat_dba, ONLY : inte_cr   
 IMPLICIT NONE
 !dummy variables
 REAL (kind=8), INTENT(IN)    :: ms,      & !mises trial stress 
                                 m2,      & !2 * mu                                 
                                 ppro(5), & !plastic properties (hardening)
                                 vpar(3), & !viscoplastic model parameters
                                 stra(2), & !efective plastic strain & evp strain rate
                                 dtime       !time increment
 INTEGER (kind=4), INTENT(IN) :: is, &
                                 np                                  
 REAL (kind=8), INTENT(INOUT) :: g    !viscoplas. consistency parameter
 INTEGER(kind=4), INTENT(OUT) :: ier  !if no convergence then is TRUE
 REAL (kind=8), POINTER :: curve(:,:) !(3,np) yield stress curve
 
 !local variables
 REAL (kind=8) :: gmax,  & !radial return plastic consist. parameter (maximum)
                  yi,    & !yield stress
                  ap,    & !hardenind modulus
                  fbar,  & !mises regularized function
                  fp,    & !1st derivative of mises regularized function
                  dg,    & !increment in viscoplastic consist parameter
                  evp,   & !viscoplastic effective strain
                  rvp      !viscoplastic effective strain rate
 INTEGER (kind=4) :: i

 REAL (kind=8)            :: a=0d0,b=1d0,c=1d0,d=0d0 ! deriv parameters
 REAL (kind=8), PARAMETER :: r23 = 0.816496580927726d0, & ! sqrt(2/3)
                             tol = 1d-6

 ! initalize some variables
 ier  = 0 ! error flag
 gmax = g ! elastic-plastic consistency parameter (maximum) from radial return
 evp  = stra(2)  ! old equivalent strain
 rvp  = stra(1)  ! old equivalent strain rate 
 g = 0d0  !
 
 DO
   !     setup initial yield FUNCTION radius
   IF( is == 5 ) THEN        !for points defined yield value
     i = 1   !begin at first interval
     yi = inte_cr (curve,np,evp,i)    !s_y
     ap = curve(3,i)                  !A'
   ELSE
     CALL isoha14(is,yi,ap,evp,ppro(1),ppro(2),ppro(3),ppro(4))  !compute s_y and A'
   END IF  
   
   ! regularized mises criteria
   fbar = ms - m2*g - r23*yi  - vpar(1)*(evp**vpar(2))*(rvp**vpar(3))
  
   IF( ABS( fbar / ms ) <= tol )EXIT
   
   ! first derivative of regularized mises function 
   IF( vpar(2) > 0d0 )THEN     ! derivative of hardening term
     a = 2d0*r23*vpar(2)*(evp**(vpar(2)-1d0))
     c = evp**vpar(2)
   END IF
   IF( vpar(3) > 0d0 )THEN     ! derivative of rate sensivity term
     b = rvp**vpar(3)
     d = (r23/dtime)*vpar(3)*(rvp**(vpar(3)-1d0))
   END IF
   ! evaluates first derivative
   fp = - m2 - (2d0/3d0)*ap - vpar(1)*(a*b+c*d) !tangent to curve fbar(vp)
   dg = fbar/fp  !increment in viscoplastic parameter
   g  = g - dg    !corrected viscoplastic parameter
   
   DO     ! bisection if g is greather gmax
     IF( g < gmax )EXIT 
     g = g / 2d0  
   END DO
      
   !update variables
   evp = stra(1) + r23 * g  !new equivalent strain
   rvp = (r23 * g) / dtime  !new equivalent strain rate   
   !rvp   = (rvp_i + vefac*rvp) / (1d0 + vefac) !new smoothed eq. strain rate
   
 END DO

 IF(g > gmax) g = gmax

 !IF (ier /= 0) WRITE(55,999) ncon,gmax,g
 !999 FORMAT('Viscoplastic parameter correction no convergence',/,  &
 !          'iter = ',     i3,/,                                   &
 !          'gmax = ', f16.12,/,                                   &
 !          'g_vp = ', f16.12,/,                                   &
 !          'CALLED FROM VISCO20',/)

 RETURN
 END SUBROUTINE visco20
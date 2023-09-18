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

 END SUBROUTINE visco20
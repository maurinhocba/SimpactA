  SUBROUTINE corr20_a (visc,t,lb3,ethm,gausv,is,np,curve,props, &
                      propv,gm,km,stres,pstr,dtime,ierr)
                    
   !
   ! plastic correction algorithm for UPDATED LAGRANGIAN FORMULATION (actual configuration)
   ! linear strain triangle (based on GARCIA GARINO thesis)
   !
   USE mat_dba, ONLY : inte_cr     
   IMPLICIT NONE
   LOGICAL, INTENT(IN) :: visc     ! TRUE if viscoplastic analysis
   REAL (kind=8), INTENT(IN) :: dtime,    & !time step size
                                t(2,2),   & !in-plane Deformation gradient
                                lb3,      & !Lambda in out of plane direction
                                ethm,     & !thermal strain trace * 2
                                props(5), & !plastic properties (hardening)
                                propv(3), & !viscoplastic properties
                                gm,km       !elastic properties (km = 3K)
   INTEGER(kind=4), INTENT(IN) :: is,       & !isotropic hardening model
                                  np          !number of points defining curve (is = 5)                                
   REAL (kind=8), INTENT(IN OUT) :: gausv(6), & !Internal variables
                                    stres(4)    !Kirchhoff stresses   
   REAL (kind=8), INTENT(OUT)    :: pstr(5)     !plastic strain tensor
   INTEGER(kind=4), INTENT(OUT)  :: ierr     !error flag
   REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress curve


 END SUBROUTINE corr20_a

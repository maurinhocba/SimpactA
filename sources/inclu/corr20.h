 SUBROUTINE corr20 (eulrf,t,lb3,gausv,is,props,propv,gm,km, &
                    stres,sigma,np,curve,ierr)
   !
   ! plastic correction algorithm for linear triangle in plane problems
   !
   USE lispa0
   USE ctrl_db, ONLY : dtime,therm
   USE mat_dba, ONLY : inte_cr
   IMPLICIT NONE
   !dummy variables for mechanical analysis
   LOGICAL, INTENT(IN) :: eulrf  !TRUE use spatial configuration else use intermediate configuration
   REAL (kind=8), INTENT(IN) :: t(2,2),   & !in-plane Deformation gradient
                                lb3,      & !Lambda in out of plane direction
                                props(5), & !plastic properties (hardening)
                                propv(3), & !viscoplastic properties
                                gm,km       !elastic properties (km = 3K)
   REAL (kind=8), INTENT(IN OUT) :: gausv(:), & !Internal variables TLF: (1:5)=Fp^-1 (6)=ep
                                                !                   ULF: (1:4)=be^-1 (5)=ep_dot (6)=ep
                                    stres(4)    !Kirchhoff stresses (for post-process)
   REAL (kind=8), INTENT(OUT)    :: sigma(4)    !for internal forces evaluation
   INTEGER(kind=4), INTENT(IN) :: is,       & !isotropic hardening model
                                  np          !number of points defining curve (is = 5)
   INTEGER(kind=4), INTENT(OUT) :: ierr     !error flag (1)
   REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress curve


 END SUBROUTINE corr20

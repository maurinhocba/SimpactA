 SUBROUTINE corr19 (t,lb3,gausv,is,props,propv,gm,km,stres,sigma, &
                    np,curve,ierr,tcont,strpl,dg,fac,alpha,dther)
   !
   ! plastic correction algorithm for linear triangle in plane problems
   ! coupled thermo-mechanical analysis
   !
   USE lispa0
   USE ctrl_db, ONLY : dtime
   USE mat_dba, ONLY : inte_cr
   IMPLICIT NONE
   !dummy variables
   REAL (kind=8), INTENT(IN) :: t(2,2),   & !in-plane Deformation gradient
                                lb3,      & !Lambda in out of plane direction
                                props(5), & !plastic properties (hardening)
                                propv(3), & !viscoplastic properties
                                gm,km,    & !elastic properties (km = 3K)
                                fac,      & !(1+nu) for plane strain and axilsymmetry
                                alpha,    & !thermal expansion coefficient
                                dther,    & !element temperature increment
                                dg          !derivative of shear modulus with respect to temperature                                 
   REAL (kind=8), INTENT(IN OUT) :: gausv(:), & !Internal variables ULF: (1:4)=be^-1 (5)=ep_dot (6)=ep
                                    stres(4), & !Kirchhoff stresses (for post-process)
                                    strpl(:)    !plastic strain tensor & elastic strain trace
   REAL (kind=8), INTENT(OUT)    :: sigma(4), & !for internal forces evaluation
                                    tcont(3)    !temperature dependant contributions for coupled thermo-mechanical 
   INTEGER(kind=4), INTENT(IN) :: is,       & !isotropic hardening model
                                  np          !number of points defining curve (is = 5)
   INTEGER(kind=4), INTENT(OUT) :: ierr     !error flag (1)
   REAL (kind=8), POINTER :: curve(:,:)     !(3,np) yield stress curve


 END SUBROUTINE corr19

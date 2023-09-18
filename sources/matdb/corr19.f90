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

   ! local variables
   REAL (kind=8) :: efps0,ethm,pstr0(4),pstr(4),etrc0,trac,estr(4),mises
   !LOGICAL :: visc !.TRUE. for elastic-visco-plastic material
   
   INTERFACE
     INCLUDE 'corr20_a.h'
   END INTERFACE   

   efps0 = gausv(6)          !initial (old) Equivalent Plastic Strain
   ethm = 2d0*fac*alpha*dther   !thermal strain trace * 2
   pstr0 = strpl(1:4) !old plastic strain tensor
   etrc0 = strpl(5)   !old elastic strain trace
   !visc = ( propv(1) > 0d0 .AND. ANY(propv(2:) > 0d0) ) !is TRUE when viscoplast material ... not yet
   !plastic corrector
   CALL corr20_a (.FALSE.,t,lb3,ethm,gausv,is,np,curve, &
                    props,propv,gm,km,stres,strpl,dtime,ierr)
   sigma = stres !for internal forces evaluation

   !evaluates thermal dissipation contributions
   trac = stres(1)+stres(2)+stres(4) ! kirchhoff stress trace
   estr(1) = stres(1)-trac/3d0 ! deviatoric kirchhoff stress
   estr(2) = stres(2)-trac/3d0 !
   estr(3) = stres(3)          !
   estr(4) = stres(4)-trac/3d0 !
   mises = SQRT(estr(1)*estr(1)+estr(2)*estr(2)+estr(4)*estr(4)+ &
                2d0*estr(3)*estr(3)) !von mises stress
   estr =  estr/gm/2d0  !new deviatoric elastic strain tensor
   !store coupled thermo-mechanical factors
   tcont    = 0d0                      !
   tcont(1) = mises*(gausv(6) - efps0) ! plastic work (s_ij : ep_dot_ij)
   tcont(2) = strpl(5) - etrc0         ! change in elastic strain trace
   pstr(1:4)= strpl(1:4) - pstr0  !change in plastic strain tensor
   tcont(3) = dther*dg*(estr(1)*pstr(1)+estr(2)*pstr(2)+ &   !change in plastic work due temperature
                        estr(4)*pstr(4)+2d0*estr(3)*pstr(3)) ! DT * dgm/dT * e_ij : ep_dot_ij

   RETURN
 END SUBROUTINE corr19



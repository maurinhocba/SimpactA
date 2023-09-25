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
   
   ! local variables
   INTEGER (kind=4) :: i
   REAL (kind=8)    :: efps0,reps0,se,devs(4),pres,trac,be(4),bo(4),fi(5), &
                       an(4),de,m2,m3,ftria,yield,aprim
   REAL (kind=8), PARAMETER :: toler = 1d-4,             &
                               r23 = 0.816496580927726d0
                               
   INTERFACE
     INCLUDE 'visco20.h'
   END INTERFACE                               

   efps0 = gausv(6) !initial (old) Equivalent Plastic Strain
   reps0 = gausv(5) !initial (old) equivalent plastic strain rate (viscoplastic  only)
   !     setup initial yield FUNCTION radius
   IF( is == 5 ) THEN        !for points defined yield value
     i = 1   !begin at first interval
     yield = inte_cr (curve,np,efps0,i)    !s_y
     aprim = curve(3,i)                    !A'
   ELSE
     CALL isoha14(is,yield,aprim,efps0,props(1),props(2),props(3),props(4))  !compute s_y and A'
   END IF   
   m2 = gm*2d0      !2*Mu to compute stresses
   m3 = gm*3d0      !3*Mu to use as auxiliar
   
   !inverse relative deformation gradient
   fi(1:2) = t(:,1)
   fi(3:4) = t(:,2)
   fi(  5) = lb3
   !compute elastic trial Finger tensor  (f^-T * b_old * f^-1)
   bo    = gausv(1:4)   !previous elastic Finger tensor
   be(1) = (fi(1)*bo(1)+fi(2)*bo(3))*fi(1) + (fi(1)*bo(3)+fi(2)*bo(2))*fi(2)
   be(2) = (fi(3)*bo(1)+fi(4)*bo(3))*fi(3) + (fi(3)*bo(3)+fi(4)*bo(2))*fi(4)
   be(3) = (fi(1)*bo(1)+fi(2)*bo(3))*fi(3) + (fi(1)*bo(3)+fi(2)*bo(2))*fi(4)
   be(4) =  fi(5)*bo(4)*fi(5)
   !elastic predictor state
   trac    = (be(1)+be(2)+be(4))/3d0 !Finger tensor trace / 3d0
   pres    =  km*(1d0-trac-ethm)/2d0 !p = 3 K * 1/3 tr(e_almansi)
   devs(1) = be(1)-trac  !deviatoric Finger tensor
   devs(2) = be(2)-trac  !
   devs(3) = be(3)       !
   devs(4) = be(4)-trac  !
   devs    = -gm*devs      !deviatoric Kirchhoff stresses
   de = 0d0                !initializes consistency parameter
   se = SQRT(devs(1)*devs(1)+devs(2)*devs(2)+devs(4)*devs(4)+ & !vMstress
         2d0*devs(3)*devs(3))                                   !
   ftria = se - r23*yield    !f_trial
   IF( ftria/se > toler )THEN  !plastic
     !radial return (no iteration) in terms of A' at eps old
     de = 1.5d0*ftria/(aprim+m3) !increment in eps
     IF( visc ) & !iteration to convergence in viscoplastic
       CALL visco20(se,m2,is,np,curve,props,propv,gausv(5:6),dtime,de,ierr)
       !de = (1d0/m2)*mises / (1d0 + (1d0/m3)*(aprim+propv(1)/dtime))
   END IF
   ! Updates internal variables
   IF( de > 0d0 ) THEN
     an   = devs / se        !plastic flow direction tensor
     devs = devs - m2*de*an  !dev Kirchhoff stress
     be   = be + 2d0*de*an   !elastic Finger tensor
     gausv(6) = efps0 + r23*de !efective plast strain
     IF( visc ) gausv(5) = r23*de / dtime !effective plast strain rate
   END IF
   !update internal variables
   gausv(1:4) = be(1:4)
   !Actualized Kirchhoff stress tensor
   stres(1) = pres + devs(1)  !tau(1,1)
   stres(2) = pres + devs(2)  !tau(2,2)
   stres(3) =        devs(3)  !tau(1,2) = tau(2,1)
   stres(4) = pres + devs(4)  !tau(3,3)

   pstr(1:4) = -de*an        !plastic deviatoric strain tensor 
   pstr(5)   = pres/(km/3d0) !elastic strain trace
 END SUBROUTINE corr20_a
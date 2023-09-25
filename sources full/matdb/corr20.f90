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

   ! local variables
   INTEGER(kind=4) :: i
   REAL (kind=8) :: yield,aprim,efps0,ethm,strpl(5)
   LOGICAL :: visc !.TRUE. for elastic-visco-plastic material

   INTERFACE
     INCLUDE 'corr20_a.h'
   END INTERFACE

   IF(eulrf)THEN !update lagrangian formulation (actual configuration) from Garcia Garino tesis
     ethm = 0d0
     visc = ( propv(1) > 0d0 .AND. ANY(propv(2:) > 0d0) ) !is TRUE when viscoplast material
     !plastic corrector
     CALL corr20_a (visc,t,lb3,ethm,gausv,is,np,curve, &
                    props,propv,gm,km,stres,strpl,dtime,ierr)
     sigma = stres !for internal forces evaluation

   ELSE !total lagrangian formulation (intermediate configuration) from Crisfield - not thermal yet
     efps0 = gausv(6) !initial (old) Equivalent Plastic Strain
     !     setup initial yield FUNCTION radius
     IF( is == 5 ) THEN        !for points defined yield value
       i = 1   !begin at first interval
       yield = inte_cr (curve,np,efps0,i)    !s_y
       aprim = curve(3,i)                    !A'
     ELSE
       CALL isoha14(is,yield,aprim,efps0,props(1),props(2),props(3),props(4))  !compute s_y and A'
     END IF      
     
     CALL corr20_i (t,lb3,gausv,aprim,yield,gm,km,stres,sigma,ierr)

   END IF

   RETURN
 END SUBROUTINE corr20
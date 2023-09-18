  SUBROUTINE ud_solid_3(stran,delta,stres,gausv,ttime,ielem,ierr, &
                       newmt,submo,props,chead)
 !-----------------------------------------------------------------------
 ! stress computation
 ! USER DEFINED MATERIAL
 !-----------------------------------------------------------------------
 USE lispa0 , ONLY : lures    !Fortran UNIT to write results (ASCII)
 USE mat_dba , ONLY : curve   !Gives access to TYPE(curve) definition
 IMPLICIT NONE                !advisable
 LOGICAL,INTENT(IN OUT):: newmt !flag to indicate that a different material
                                !than previous one is used.
                                !This variable can be used to initialize
                                !constant, and must be changed to .FALSE.
 ! The next 5 lines give access to the data read in arrays
 INTEGER (kind=4),INTENT(IN):: submo   !material model
 REAL(kind=8), POINTER :: props(:)     ! array of properties
 TYPE (curve), POINTER :: chead  !pointer to first curve
 ! The next 2 lines include control data
 INTEGER(kind=4),INTENT(IN):: ielem !element label, to print error messages
 REAL(kind=8),INTENT(IN):: ttime    !Process time
 ! The next 3 lines include the most relevant variables
 REAL(kind=8),INTENT(IN):: stran(4) !twice the total shear strains in local system
 REAL(kind=8),INTENT(IN):: delta    !total volumentric strain
 REAL(kind=8),INTENT(IN OUT):: gausv(*) !Internal variables
 REAL(kind=8),INTENT(OUT):: stres(4) !total stresses in local system
 ! The next line include an error flag (output)
 INTEGER(kind=4),INTENT(OUT):: ierr !error flag, default value is .FALSE.
                                    !if changed to .TRUE. it will stop
                                    !the execution
 !Local variables
 REAL(kind=8) :: strpl(4),str(4),efpst,press
 REAL(kind=8), SAVE  :: gm,kv,death
 REAL(kind=8), POINTER  :: val(:,:)
 LOGICAL, SAVE :: yiel,pflag
 INTEGER (kind=4), SAVE :: isoh,kinh,np

 INTERFACE
   INCLUDE 'corr17.h'
 END INTERFACE

 IF( newmt )THEN
   yiel = props(3) > 0
   IF( yiel )THEN

     ! isotropic hardening
     IF( props(6) > 0 )THEN
       isoh = 4                   !exponential + saturation
     ELSE IF( props(5) > 0 )THEN
       isoh = 3                   !Ludwik Nadai
     ELSE IF( props(4) > 0 )THEN
       isoh = 2                   !linear
     ELSE
       isoh = 1                   !no hardening
     END IF

     ! kinematic hardening
     IF(props(8) > 0 )THEN
       kinh = 2         !linear
       IF(props(9) > 0 ) kinh = 3 !non-linear
     ELSE
       kinh = 1         !no kinematic hardening
     END IF
   END IF
   ! elastic properties
   gm = props(1)
   kv = 3d0*props(2)
   death = props(7)
   np = 0          !not a curve
   NULLIFY(val)
 END IF

 IF (kinh == 1 )THEN  ! no kinematic hardening
   !elastic (trial) strains
   strpl(1:4) = gausv(1:4)               !previous (twice) plastic strains
   str(1:4) = stran(1:4) - strpl(1:4)      !trial Elastic log strains

   stres= gm*str                       !Trial elastic shear stresses
   efpst = gausv(5)                    !effect plastic strain
   CALL corr17(stres(1),stres(2),stres(3),stres(4),efpst,gm, &
               props(3:6),props(10:16),ierr,strpl,pflag,isoh,ielem,np,val)
   IF(ierr == 1) RETURN              !no convergence in plasticity
   IF( pflag )THEN                   !if plastic flow
     gausv(1:4) = gausv(1:4) + strpl(1:4)  !total plastic shear strains
     gausv(5)   = gausv(5)   + efpst       !Total effect. plastic strain
   END IF
   press = kv*delta
   stres(1) = stres(1) + press
   stres(2) = stres(2) + press
   stres(4) = stres(4) + press
 ELSE
   !no yet
 END IF
 RETURN
 END SUBROUTINE ud_solid_3

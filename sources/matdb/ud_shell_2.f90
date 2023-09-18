  SUBROUTINE ud_shell_2(stran,stres,gausv,ttime,ielem,ierr, &
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
 REAL(kind=8),INTENT(IN):: stran(3) !total strains in local system
 REAL(kind=8),INTENT(IN OUT):: gausv(:) !Internal variables
 REAL(kind=8),INTENT(OUT):: stres(3) !total stresses in local system
 ! The next line include an error flag (output)
 INTEGER(kind=4),INTENT(OUT):: ierr !error flag, default value is .FALSE.
                                    !if changed to .TRUE. it will stop
                                    !the execution

   !! Derived type for a curve
   !TYPE curve
   !  INTEGER (kind=4) :: np            !number of points
   !  CHARACTER (len=mnam) :: name        !Associated label
   !  REAL(kind=8), POINTER :: val(:,:) !x-y values
   !  TYPE(curve), POINTER :: next      !pointer for next curve
   !END TYPE curve

   !Local variables
   REAL(kind=8):: efpst, khard(3),st(3)
   LOGICAL, SAVE :: yiel
   INTEGER (kind=4), SAVE :: isoh,kinh
   REAL (kind=8), SAVE :: death

   IF( newmt )THEN
     yiel = props(5) > 0
     IF( yiel )THEN

       ! isotropic hardening
       IF( props(8) > 0 )THEN
         isoh = 4                   !exponential + saturation
       ELSE IF( props(7) > 0 )THEN
         isoh = 3                   !Ludwik Nadai
       ELSE IF( props(6) > 0 )THEN
         isoh = 2                   !linear
       ELSE
         isoh = 1                   !no hardening
       END IF

       ! kinematic hardening
       IF(props(10) > 0 )THEN
         kinh = 2         !linear
         IF(props(11) > 0 ) kinh = 3 !non-linear
       ELSE
         kinh = 1         !no kinematic hardening
       END IF
     END IF
     death = props(9)
     newmt = .FALSE.
   END IF

   IF( yiel )THEN
     st = stran - gausv(1:3)
   ELSE
     st = stran
   END IF

   stres(1) =  props(1)* st(1) +  props(2)* st(2)
   stres(2) =  props(2)* st(1) +  props(3)* st(2)
   stres(3) =  props(4)* st(3)

   IF( yiel )THEN  !if plastic then
     IF( death >= ttime) THEN !if plasticity active
       efpst = gausv(4)       !equivalent plastic strain
       IF( kinh > 1 ) THEN    !if kinematic hardening
         khard(1:3) = gausv(5:7)      !back stress (in local coordinates)
         !CALL corr31(str(1),khard(1)),efpst,props(1),props(5),props(10), &
         !            props(11),propp(12),propp(18),ierr,isoh,kinh,st(1))
         gausv(5:7) = khard(1:3)      !effective stress
       ELSE
         CALL ud_corr_30(stres(1),stres(2),stres(3),efpst,props(1), &
                         props(5),props(18),props(12),ierr,isoh,st(1))
       END IF
       IF(ierr == 1) THEN
         WRITE(55,   "(' Error in the constitutive equation at element',i7)",ERR=9999) ielem
         WRITE(lures,"(' Error in the constitutive equation at element',i7)",ERR=9999) ielem
         WRITE(*,    "(' Error in the constitutive equation at element',i7)")ielem
         RETURN
       END IF
       IF( efpst > 0 )THEN
         gausv(4)   = gausv(4) + efpst    !stores new equivalent plastic strain
         gausv(1:3) = gausv(1:3) + st     !stores new plastic strain
       END IF
     END IF
   END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE ud_shell_2

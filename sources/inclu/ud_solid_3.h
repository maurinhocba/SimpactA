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
  END SUBROUTINE ud_solid_3

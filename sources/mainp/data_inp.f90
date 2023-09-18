 SUBROUTINE data_inp(actio)

 USE param_db,ONLY : milb,mich				!Global parameters
 USE ctrl_db, ONLY : nstra,text,therm, aero	!global control parameters
 USE outp_db, ONLY : iwrit					!output asked
 USE export_db, ONLY : exp_data				!export data
 USE damp_db, ONLY : rddamp
 IMPLICIT NONE
   !dummy arguments
   CHARACTER(len=milb), INTENT(IN OUT) :: actio
   !Functions
   CHARACTER(len=mich):: inttoch   !Puts an integer into a string

    !           Control parameters input

    CALL timing(3,1)                     !init time for Control data input
    CALL contol(actio)                   !read Main control parameters
    CALL timing(3,2)                     !end time for Control data input
    WRITE(6,"(/,5X,'***  BEGIN STAGE ',A,': ',A,/)") TRIM(inttoch(nstra,0)),text
    WRITE(55,"(/,5X,'***  BEGIN STAGE ',A,': ',A,/)",err=9999)TRIM(inttoch(nstra,0)),text

    ! modify ACTIO to major changes (NSTRA2) if CUTTING or and added set

    !           Finite element mesh data input
    CALL timing(6,1)                     !init time for mesh data input
    CALL inpdat(actio )                  !read Mesh data and constraints
    CALL timing(6,2)                     !end time for mesh data input
    
    !			Aerodynamical data input
    IF(aero .AND. TRIM(actio)=='NEW') THEN
      CALL timing(20,1)                   !init time for aero data input
      CALL aero_data(actio)               !read Aerodynamical Mesh and conditions data
      CALL timing(20,2)                   !end time for aero data input
    END IF
    !           Initial conditions and time history data input
    CALL timing(7,1)                     !init time for initial conditions data input
    CALL histor( )                       !read data for history curves
    CALL rddamp( iwrit )                 !read damping DATA
    CALL timing(7,2)                     !end time for initial conditions data input

    !           Load data input
    CALL timing(9,1)                     !init time for load data input
    CALL loadpl( )                       !read load data
    IF( therm )CALL heatpl ( actio )     !read external heat data (thermal problem)
    CALL timing(9,2)                     !end time for load data input

    !          Contact data input
    CALL timing(4,1)
    CALL conini(actio)
    CALL timing(4,2)

    !           Initial conditions  !!  this must be revised
    CALL timing(7,1)
    CALL intime( )   !read initial conditions (disp, vel, stresses)
    CALL timing(7,2)

    !          Export data input
    CALL timing(3,1)
    CALL exp_data ( )                    !read EXPort DATA for future use
    CALL timing(3,2)

    CALL endstr( )      ! end of input data for one strategy
    RETURN
    9999 CALL runen2('ERROR WHILE WRITING "DEB" FILE.')
 END SUBROUTINE data_inp

 SUBROUTINE initial_comp(actio,istop)

 USE param_db,  ONLY : milb              !Global parameters
 USE ctrl_db,   ONLY : numct,ttime,dtime,memo,therm,tmscal,lumped,aero  !global control parameters
 USE npo_db,    ONLY : resid,tmass,tresi       !point information
 USE outp_db,   ONLY : iwrit             !output asked
 USE gvar_db,   ONLY : actchk
 USE inter_sr,  ONLY : inter_ini
 USE introut_sr,ONLY : introut_mng
 !USE gener_db,  ONLY : include_damp
 
 IMPLICIT NONE
   !dummy arguments
   CHARACTER(len=milb), INTENT(IN OUT) :: actio
   INTEGER(kind=4), INTENT(IN OUT) :: istop

   INTERFACE
     INCLUDE 'contac.h'
     INCLUDE 'elemnt.h'
   END INTERFACE

    CALL timing(6,1)
    CALL eqnums(actio)   !compute active DOF numbers
    CALL timing(6,2)

    !           Evaluating gauss point variables

    CALL timing(5,1)
    CALL elemnt('GAUSSC', istop=istop, flag1=.TRUE.)
    IF (istop == 1) THEN !istop=1 -> error in jacobian calculation
       WRITE(6 ,901)
       WRITE(55,901,ERR=9999)
       STOP
    END IF
    CALL timing(5,2)

    !           Evaluating lumped/consistent mass matrix

    CALL timing(10,1)
    CALL masmtx( )
    CALL wrtpos(0)  ! write mesh data for post-process
    CALL timing(10,2)

    !           Reference Conservative Load Calculation

    CALL timing(9,1)
    CALL loaini( )
    CALL timing(9,2)

    IF (therm) THEN
      CALL tlmass ( )              !compute capacity matrix
      tmass = tmass*tmscal         !scale capacity matrix
      CALL heaini (actio)          !reference heating Calculation
    END IF
	
	!			Aeroelastic Interaction
	
	IF( aero ) THEN
      CALL timing(1,1)
	  IF( TRIM(actio)=='NEW' )THEN  ! new problem
        CALL inter_ini              ! compute master nodes and initial distancies
        CALL introut_mng('initial')
      END IF
      CALL introut_mng('openFiles',actio)
      !CALL include_damp
      CALL timing(1,2)
	END IF

    !           Calculation of internal forces due to initial stresses

    IF( .NOT. actchk )THEN
      CALL timing(14,1)
      resid = 0d0
      CALL elemnt('RESVPL', deltc=dtime, istop=istop, ttime= ttime)
      IF(therm) THEN
        tresi = 0d0
        CALL elemnt ('TRESID', deltc=dtime, ttime= ttime)
      END IF
      CALL timing(14,2)
      IF (istop == 1) THEN !istop=1 -> there was no convergence
         WRITE(6 ,900)
         WRITE(55,900,err=9999)
         STOP
      END IF

      IF (numct > 0) THEN
         CALL contac('INITIA', iwrit=iwrit, maxve=memo, ttime=ttime)

      END IF


      WRITE(6,"(/,18x,'*** PRELIMINARY CALCULATIONS COMPLETED ***'/ &
              & /,23x,'*** TIME INTEGRATION STARTED ***',/)")
    END IF


    RETURN
    9999 CALL runen2('ERROR WHILE WRITING "REP" FILE.')
    900 FORMAT(/,12x,' DIVERGENCE IN THE RETURN ALGORITHM : ', &
                     'PROGRAM STOPPED',/,17x,' FINAL GEOMETRY ', &
                     'SAVED FOR POSTPROCESSING!')
    901 FORMAT(' PROGRAM STOPPED DUE TO A DETECTED ERROR'// &
               ' IN THE JACOBIAN CALCULATION. MAKE A GEOMETRY REVIEW !')
  END SUBROUTINE initial_comp

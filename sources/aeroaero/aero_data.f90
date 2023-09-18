SUBROUTINE aero_data(actio)

  USE inter_db,     ONLY: bEType
  !USE wind_sr,      ONLY: rdwind
    USE wind_db, ONLY: windFlag
  !USE control_sr,   ONLY: rdctrl

  IMPLICIT NONE

  CHARACTER(*) :: actio
  
  
  CALL Input_simulation     ! subroutine not included in a module
  CALL Input_extras         ! subroutine not included in a module
  CALL Input_geometry       ! subroutine not included in a module
  
  IF(TRIM(bEType)=='NBST')THEN ! shell model for blades
    CALL Input_bed          ! subroutine not included in a module
  ENDIF
  
  !CALL rdwind              NOT IN USE, IN ORDER TO PREVENT INCLUDING ERRORS
    windFlag = .FALSE.
  !CALL rdctrl              NOT IN USE, IN ORDER TO PREVENT INCLUDING ERRORS
  !CALL gener_inp           NOT IN USE, IN ORDER TO PREVENT INCLUDING ERRORS

END SUBROUTINE aero_data

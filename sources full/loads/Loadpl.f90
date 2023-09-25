SUBROUTINE loadpl ( )

  !*** load input routine

  USE ctrl_db, ONLY: ndime, ndofn, nload
  USE outp_db, ONLY: iwrit
  USE c_input, ONLY: lures, exists, listen, backs
  IMPLICIT NONE

  INTERFACE
    INCLUDE 'rdload.h'
  END INTERFACE

  CALL listen('LOADPL')              !read first card
  IF (.NOT.exists('LOADDA')) THEN      !if key-word LOAD_DATA not found
    backs = .TRUE.                   !one card back
    IF ( nload > 0 ) &
      WRITE (lures,"('Load from the previous strategy used.')",ERR=9999)
  ELSE                               !LOAD_DATA present
                                       !Header
    IF (iwrit == 1) WRITE (lures, &
       "(/,'  L O A D   R E A D   I N   T H E   ', &
       &   'P R E S E N T   S T R A T E G Y ',/)",ERR=9999)
    CALL rdload (iwrit,ndime,ndofn,nload)

    CALL listen('LOADPL')
    IF (.NOT.exists('ENDLOA'))CALL runend('LOADPL: END_LOAD_DATA CARD EXPECTED')

  END IF

  RETURN
 9999 CALL runen2('')
END SUBROUTINE loadpl

SUBROUTINE heatpl (actio)

  !*** heat input routine

  USE curv_db, ONLY: curnam,getcun
  USE ctrl_db, ONLY: ndime, ndoft, nheat, npoin, npoio, itemp
  USE outp_db, ONLY: iwrit
  USE c_input, ONLY: ludat, lures, exists, param, listen, get_name, backs
  USE heat_db
  USE npo_db, ONLY: label,oldlb,tempe,dtemp
  IMPLICIT NONE

  CHARACTER(len=* ), INTENT(IN) :: actio

  INTERFACE
    INCLUDE 'rdheat.h'
  END INTERFACE

  CALL listen('HEATPL')              !read first card
  IF (.NOT.exists('HEATDA')) THEN      !if key-word HEAT_DATA not found
    backs = .TRUE.    !one card back
    IF (actio == 'NEW') THEN           !for a new problem
      nheat = 0                      !no external heats considered
    ELSE                               !for a new strategy
      WRITE (lures,"('Heat from the previous strategy used.')",ERR=9999)
    END IF

  ELSE                               !HEAT_DATA present
    IF (actio == 'NEW') THEN           !for a new problem, initializes counters
      nheat = 0                        !number of different heat sets
    END IF
    CALL rdheat (iwrit,ndoft,nheat,actio) !input or delete heat set data
    CALL listen('HEATPL')
    IF (.NOT.exists('ENDHEA'))CALL runend('HEATPL: END_HEAT_DATA CARD EXPECTED')

  END IF

  RETURN
9999 CALL runen2('')
END SUBROUTINE heatpl

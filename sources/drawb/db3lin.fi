SUBROUTINE db3lin (ndbpoin,lnd, xl, coora)

  ! compute db line connectivies (internal) and coordinates (x-y)
  ! checks coincident nodes

  !   m o d u l e s
  USE lispa0
  USE line_db
  USE outp_db, ONLY : iwrit
  USE npo_db, ONLY : label

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s
  INTEGER, INTENT (IN) :: ndbpoin  !number of nodes in the db line
  REAL (kind=8), INTENT (IN) :: coora (:,:)    !present coordinates
  INTEGER (kind=4), INTENT (OUT) :: lnd (:,:)  !db line connectivities
  REAL (kind=8), INTENT (OUT) :: xl (:,:)      !db line nodes coords

  !   l o c a l   v a r i a b l e s
  INTEGER (kind=4) :: i, n, chnode, ipoin, ndbpoin1
  TYPE (lsg), POINTER :: seg
  !-----------------------------------------------------------------

  ndbpoin1 = ndbpoin - 1                   !number of segments
  WRITE (lures,260,ERR=9999) ndbpoin       !echo
  IF(iwrit == 1)WRITE (lures,261,ERR=9999) !echo heading
  seg => headli                            !point to first segment
  DO n = 1,ndbpoin1                     !loop over segments
    ipoin = chnode (seg%nods(1))           !first node
    xl (1:2, n) = coora (1:2, ipoin)       !x-y coordinates
    lnd(1,n) = ipoin                       !keep node in array
    ipoin = chnode (seg%nods(2))           !second node
    lnd(2,n) = ipoin                       !keep node in array
    IF(iwrit == 1)WRITE (lures,"(2i8)",ERR=9999)seg%nods(1:2)  !echo
    seg => seg%next
  END DO
  xl (1:2, ndbpoin) = coora (1:2, ipoin)   !keep last node coords
  ! check coincident nodes
  DO i = 1, ndbpoin1       !for each segment
     IF (xl(1, i) == xl(1, i+1)) THEN    !compare first and last node
        IF (xl(2, i) == xl(2, i+1)) THEN  !write message and stop
           WRITE (0, 320,ERR=9999)     i, i + 1
           WRITE (lures, 320,ERR=9999) i, i + 1
           STOP
        END IF
     END IF
  END DO

  RETURN
260 FORMAT (/ ' Drawbead-sheet pair definition', /, /, &
             ' Number of drawbead line nodes ...', i8)
261 FORMAT (/' List of drawbead line nodes ' /)

320 FORMAT (/, /, '        **********  Fatal error ********** ', /, &
                 ' Nodes', i8, ' and', i8,  &
             ' have the same coordinates.', /, /)
 9999 CALL runen2('')
END SUBROUTINE db3lin

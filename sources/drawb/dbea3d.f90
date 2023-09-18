SUBROUTINE dbea3d (itask, ttime, flag)

  !  3d drawbead main routine

  USE ctrl_db, ONLY: nconp !,ndofn, npoin, npoio
  USE db3_db !, ONLY : ndbea,db3inp,db3for,db3dmp,db3res,db3udl,db3ini,prdbln,db3pair_db
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: itask  !task to perform
  REAL  (KIND=8),    INTENT(IN) :: ttime
  INTEGER (kind=4), INTENT(IN) :: flag
  !  local variables

  INTEGER(kind=4), SAVE :: ndold = 0
  TYPE(db3pair_db),POINTER:: pair
!======================================================================

  SELECT CASE (TRIM(itask))

  CASE ('NEW', 'NSTRA0', 'NSTRA1', 'NSTRA2')

    CALL db3inp (ndold)

  CASE ('FORCES')

    CALL db3for ( )

  CASE ('DUMPIN')

    CALL db3dmp ( itask )

  CASE ('RESTAR')
    ndbea = nconp(4)
    CALL db3res( itask )

  CASE ('UPDLON')
    CALL db3udl ( )
    nconp(4) = ndbea

  CASE ('INITIA')

    CALL db3ini (ndold)

  CASE ('WRTSUR')
    ! writes drawbead data for history
    CALL prdbln(flag)

  CASE ('OUTDY1')
    !  writes drawbead FORCES
    IF( numtdf > 0) THEN
      WRITE (42,ERR=9999) ttime
      pair => headp
      DO
        IF( .NOT.ASSOCIATED (pair)) EXIT
        IF( pair%iparm(6) == 1 )WRITE (42,ERR=9999) pair%tdbf
        pair => pair%next
      END DO
    END IF

  CASE ('OUTDY2')
    !  writes drawbead Affected elements
    WRITE (43,ERR=9999) ttime
    pair => headp
    DO
      IF( .NOT.ASSOCIATED (pair)) EXIT
      WRITE(43,err=9999) pair%daz
      pair => pair%next
    END DO

  END SELECT

  RETURN
 9999 CALL runen2('')
END SUBROUTINE dbea3d

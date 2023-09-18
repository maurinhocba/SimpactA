      SUBROUTINE wrtfl1 (ttime)

      USE  loa_db
      IMPLICIT NONE
      REAL (kind=8) :: ttime

      !   l o c a l   v a r i a b l e s

      LOGICAL :: wrtn
      TYPE (loa_set), POINTER :: loas
      !---------------------------------------------

      IF (ifloa == 0) RETURN

      wrtn = .FALSE.
      loas => headls
      DO 
        IF(.NOT.ASSOCIATED(loas)) EXIT
        IF ( loas%numfl > 0 .AND. loas%flpar%nodcen /= 0) THEN
           IF (.NOT.wrtn)THEN
             WRITE (52,ERR=9999) ttime
             wrtn = .TRUE.
           ENDIF
           WRITE (52,ERR=9999) loas%flpar%vol, loas%flpar%press
        END IF
        loas => loas%next
      END DO

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE wrtfl1

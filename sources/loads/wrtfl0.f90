      SUBROUTINE wrtfl0

      !writes information for post-processing associated to follower loads

      USE  loa_db
      USE c_input, ONLY : openfi
      IMPLICIT NONE

      !   l o c a l   v a r i a b l e s

      INTEGER ::  nvd
      TYPE (loa_set), POINTER :: loas
      !---------------------------------------------
      nvd = 0

      loas => headls
      DO ! counts volume dependent follower loads
        IF(.NOT.ASSOCIATED(loas)) EXIT
        IF ( loas%numfl > 0 .AND. loas%flpar%nodcen /= 0) nvd = nvd + 1
        loas => loas%next
      END DO

      IF (nvd > 0) THEN
         CALL openfi(52)
         WRITE (52,ERR=9999) nvd, 2 ! volume + pressure = 2
         loas => headls
         DO
           IF(.NOT.ASSOCIATED(loas)) EXIT
           IF ( loas%numfl > 0 .AND. loas%flpar%nodcen /= 0) THEN
              WRITE (52,ERR=9999) loas%lbl
           END IF
           loas => loas%next
         END DO
      END IF

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE wrtfl0

      SUBROUTINE endstr
!********************************************************************
!
!*** END STRATEGY card READ
!
!********************************************************************
      USE c_input, ONLY : listen,exists
      IMPLICIT NONE

      CALL listen('ENDSTR')
      IF (.NOT.exists('ENDSTR'))  CALL runend('ENDSTR: END_STRATEGY CARD EXPECTED ')

      RETURN
      END SUBROUTINE endstr

      SUBROUTINE explit(azero, dtend, dtime, dtold, dtrec, ttime, flag)
!********************************************************************
!
! *** time stepping routine
!
!********************************************************************
      IMPLICIT NONE

      LOGICAL, INTENT (IN OUT) :: flag
      REAL (kind=8),INTENT(IN) :: azero, dtrec, dtime, dtend
      REAL (kind=8),INTENT(IN OUT) :: ttime,dtold

      END SUBROUTINE explit

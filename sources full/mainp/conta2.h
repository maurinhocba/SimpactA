      SUBROUTINE conta2(itask,dtcal,ttime,iwrit,velnp,maxve)

!     main contac routine (ALGORITHM 2)

      IMPLICIT NONE

!        Dummy arguments
                                     !task to perform
      CHARACTER(len=*),INTENT(IN) :: itask
      INTEGER (kind=4),INTENT(IN) :: iwrit,maxve
      REAL (kind=8),INTENT(IN), OPTIONAL ::dtcal,ttime,velnp(:,:)
      END SUBROUTINE conta2

      SUBROUTINE conta1(itask,dtcal,ttime,iwrit,      &
     &                  velnp,maxve)

!     main contac routine (TUBE IN A SKIN)

      IMPLICIT NONE

!        Dummy arguments
                                     !task to perform
      CHARACTER(len=*),INTENT(IN) :: itask
      INTEGER (kind=4),INTENT(IN) :: iwrit,maxve
      REAL (kind=8),INTENT(IN), OPTIONAL ::dtcal,ttime,velnp(:,:)
      END SUBROUTINE conta1

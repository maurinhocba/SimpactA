      SUBROUTINE contac(task,iwrit,dtcal,ttime,velnp,maxve)

  !     main contac routine

      IMPLICIT NONE

  !        Dummy arguments
      CHARACTER(len=*),INTENT(IN) :: task  !task to perform
      INTEGER (kind=4),INTENT(IN) :: iwrit
      INTEGER (kind=4),INTENT(IN), OPTIONAL :: maxve
      REAL (kind=8),INTENT(IN), OPTIONAL :: dtcal,ttime,velnp(:,:)
      END SUBROUTINE contac

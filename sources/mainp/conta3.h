      SUBROUTINE conta3(itask,dtcal,ttime,iwrit,     &
     &                  velnp,maxve)

      IMPLICIT NONE

!        Dummy arguments
                                     !task to perform
      CHARACTER(len=*),INTENT(IN):: itask
      INTEGER (kind=4),INTENT(IN), OPTIONAL :: iwrit,maxve
      REAL (kind=8),INTENT(IN), OPTIONAL ::dtcal,ttime,velnp(:,:)

      END SUBROUTINE conta3

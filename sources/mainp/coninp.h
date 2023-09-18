      SUBROUTINE coninp(itask,ncont,iwrit,dtcal,ttime,  &
     &                  velnp,maxve)

!     main contac input routine

      IMPLICIT NONE

!        Dummy arguments
      CHARACTER(len=*),INTENT(IN) :: itask  !task to perform
      INTEGER (kind=4),INTENT(IN) :: ncont
      INTEGER (kind=4),INTENT(IN) :: iwrit
      INTEGER (kind=4),INTENT(IN), OPTIONAL :: maxve
      REAL (kind=8),INTENT(IN), OPTIONAL :: dtcal,ttime,velnp(:,:)
      END SUBROUTINE coninp

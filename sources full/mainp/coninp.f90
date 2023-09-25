      SUBROUTINE coninp(itask,ncont,iwrit,dtcal,ttime,  &
     &                  velnp,maxve)

!     main contac input routine

      IMPLICIT NONE

!        Dummy arguments
      CHARACTER(len=*),INTENT(IN):: itask  !task to perform
      INTEGER(kind=4),INTENT(IN):: ncont
      INTEGER(kind=4),INTENT(IN):: iwrit
      INTEGER(kind=4),INTENT(IN),OPTIONAL:: maxve
      REAL(kind=8),INTENT(IN),OPTIONAL:: dtcal,ttime,velnp(:,:)
!        Local variables

      INTERFACE
        INCLUDE 'conta1.h'
        INCLUDE 'conta2.h'
        INCLUDE 'conta3.h'
        INCLUDE 'dbea3d.h'
      END INTERFACE

      SELECT CASE (ncont)

      CASE (1) ! special contact algorithm - TUBE IN A SKIN
        CALL conta1(itask,dtcal,ttime,iwrit,velnp,maxve)
      CASE (2) ! general contact in 2-D
        CALL conta2(itask,dtcal,ttime,iwrit,velnp,maxve)
      CASE (3) ! general contact in 3-D
        CALL conta3(itask,dtcal,ttime,iwrit,velnp,maxve)
      CASE (4) ! special contact algorithm - drawbeads
        CALL dbea3d(itask,ttime,iwrit)
      CASE DEFAULT
        CALL runend (' Wrong contact algorithm code')

      END SELECT

      RETURN
      END SUBROUTINE coninp

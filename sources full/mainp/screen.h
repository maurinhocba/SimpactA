      SUBROUTINE screen(istep,dtime,ttime,cpui,endtm,tflag)
      IMPLICIT NONE
      INTEGER (kind=4),INTENT(IN) :: istep
      REAL (kind=8),INTENT(IN) :: dtime,ttime,cpui,endtm
      LOGICAL,INTENT(INOUT):: tflag
      END SUBROUTINE screen

      SUBROUTINE wrtrpt(task,narg,params,text,newnum,oldnum,vint,nrea,blflg)
      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: task
      INTEGER(kind=4),INTENT(IN),OPTIONAL:: narg, newnum, oldnum,        &
     &                vint(:)
      REAL(kind=8),INTENT(IN),OPTIONAL:: nrea
      CHARACTER(len=*),INTENT(IN),OPTIONAL:: params(:), text
      LOGICAL,INTENT(IN),OPTIONAL:: blflg

      END SUBROUTINE wrtrpt

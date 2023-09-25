      SUBROUTINE elmdat(task,nelem,elsnam,itype)

!     read element sets

      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: elsnam
      CHARACTER(len=*),INTENT(IN):: task
      INTEGER (kind=4),INTENT(IN):: itype
      INTEGER (kind=4),INTENT(IN OUT):: nelem

      END SUBROUTINE elmdat

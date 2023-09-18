      SUBROUTINE closef(nreqd,nreqv,nreqa,nreqc,ncont)

      IMPLICIT NONE
      INTEGER (kind=4),INTENT(IN) :: nreqd,nreqv,nreqa,nreqc,ncont

      INTERFACE
        INCLUDE 'elemnt.h'
      END INTERFACE

      IF(nreqd > 0)CLOSE(11,STATUS='KEEP')
      IF(nreqc > 0)CLOSE(14,STATUS='KEEP')
      CLOSE(16,STATUS='KEEP')
      CLOSE(17,STATUS='KEEP')
      IF(nreqv > 0)CLOSE(19,STATUS='KEEP')
      IF(nreqa > 0)CLOSE(20,STATUS='KEEP')
      IF(ncont ==1)CLOSE(30,STATUS='KEEP')

!      CALL elemnt ('CLOSEF',eset,esets)

      RETURN
      END SUBROUTINE closef

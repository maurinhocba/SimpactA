      SUBROUTINE close1(nreqs,narch)

      IMPLICIT NONE
      INTEGER(kind=4) nreqs,narch

      IF(nreqs > 0)CLOSE(narch,STATUS='KEEP')
      RETURN
      END SUBROUTINE close1

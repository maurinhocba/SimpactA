      SUBROUTINE elem10(TASK, nelms, elsnam, istop, flag1, flag2)

      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters

      LOGICAL, OPTIONAL :: flag1, flag2
      CHARACTER (len=*), OPTIONAL  :: elsnam
      INTEGER (kind=4), OPTIONAL  :: istop, nelms(:)

      END SUBROUTINE elem10

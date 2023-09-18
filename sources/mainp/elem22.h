      SUBROUTINE elem22 (TASK, elsnam, dtime,                    &
     &                   ttime, istop, flag2 )

      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters

      LOGICAL, OPTIONAL :: flag2
      CHARACTER (len=*), OPTIONAL :: elsnam
      INTEGER (kind=4), OPTIONAL :: istop
      REAL (kind=8), OPTIONAL :: dtime,ttime

      END SUBROUTINE elem22


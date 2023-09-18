      SUBROUTINE elemt2(TASK, nelms, elsnam, dtime, ttime, istop, &
     &                  flag1, flag2)


      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters

      LOGICAL, OPTIONAL :: flag1,flag2
      CHARACTER (len=*), OPTIONAL :: elsnam
      INTEGER (kind=4), OPTIONAL :: istop, nelms(:)
      REAL (kind=8), OPTIONAL :: dtime,ttime

      END SUBROUTINE elemt2

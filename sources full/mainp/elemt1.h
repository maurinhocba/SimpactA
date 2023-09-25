      SUBROUTINE elemt1(TASK, nelms, elsnam, dtime, ttime, flag1, flag2)

      !USE ctrl_db, ONLY: iwrit, ndime, neulr, ndofn, npoin
      !USE outp_db, ONLY: sumat
      !USE ele01_db
      !USE lispa0
      !USE npo_db

      IMPLICIT NONE
      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters

      LOGICAL, OPTIONAL :: flag1,flag2
      CHARACTER (len=*), OPTIONAL :: elsnam
      INTEGER (kind=4), OPTIONAL :: nelms(:)
      REAL (kind=8), OPTIONAL :: dtime,ttime
      END SUBROUTINE elemt1

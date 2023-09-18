      SUBROUTINE elem18(TASK, nelms, elsnam, dtime, ttime, istop,       &
     &                  ivect, lnod, flag1, flag2)

      !master routine for element 18 (TLF) 2-D solid element

      USE param_db,ONLY: mnam
      !USE ctrl_db, ONLY: iwrit, ndofn, npoin
      !USE outp_db, ONLY: sumat
      !USE ele18_db
      !USE npo_db

      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters

      LOGICAL, OPTIONAL :: flag1,flag2
      CHARACTER (len=*), OPTIONAL :: elsnam
      INTEGER (kind=4), OPTIONAL :: istop, ivect(:), nelms(:)
      INTEGER (kind=4), POINTER, OPTIONAL :: lnod(:,:)
      REAL (kind=8), OPTIONAL :: dtime,ttime

      END SUBROUTINE elem18

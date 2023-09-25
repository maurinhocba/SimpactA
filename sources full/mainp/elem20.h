      SUBROUTINE elem20(TASK, elsnam, dtime, ttime, istop,             &
     &                  ivect, lnod, flag1, flag2)

      !master routine for element 20 (TLF) 2-D solid triangular element

      USE param_db,ONLY: mnam
      USE ctrl_db, ONLY: ndofn, npoin
      USE outp_db, ONLY: sumat, iwrit
      USE ele20_db
      USE c_input
      !USE meshmo_db
      USE npo_db

      IMPLICIT NONE

      CHARACTER(len=*),INTENT(IN):: TASK

      ! optional parameters

      LOGICAL, OPTIONAL :: flag1,flag2
      CHARACTER (len=*), OPTIONAL :: elsnam
      INTEGER (kind=4), OPTIONAL :: istop, ivect(:)
      INTEGER (kind=4), POINTER, OPTIONAL :: lnod(:,:)
      REAL (kind=8), OPTIONAL :: dtime, ttime

      END SUBROUTINE elem20

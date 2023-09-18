SUBROUTINE dbea3d (itask, ttime, flag)

  !  3d drawbead main routine

  USE ctrl_db, ONLY: ndofn, npoin, npoio
  USE npo_db
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: itask  !task to perform
  REAL  (KIND=8),    INTENT(IN) :: ttime
  INTEGER(KIND=4),   INTENT(IN) :: flag

END SUBROUTINE dbea3d

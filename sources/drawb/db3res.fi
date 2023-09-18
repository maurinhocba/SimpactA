SUBROUTINE db3res (itask )

!     restore 3d drawbead database

  USE c_input, ONLY : openfi
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT (IN) :: itask

  ! local
  INTEGER (KIND=4) :: n

  TYPE (db3pair_db), POINTER :: pair

  READ (51) numtdf, puncht

  DO n = 1, ndbea
     CALL new_db3(pair) !ALLOCATE (pair)
     CALL rsdbe3 (pair, itask)
     CALL add_db3pair (pair,headp,tailp)
  END DO

  RETURN

END SUBROUTINE db3res

SUBROUTINE db3dmp (itask)

!     dump 3d drawbead database

  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT (IN) :: itask

  ! local
  INTEGER (KIND=4) :: n

  TYPE (db3pair_db), POINTER :: pair


  WRITE (50,ERR=9999) numtdf, puncht

  pair => headp

  DO n = 1, ndbea
     CALL rsdbe3 (pair, itask)
     pair => pair%next
  END DO

  RETURN
 9999 CALL runen2('')
END SUBROUTINE db3dmp

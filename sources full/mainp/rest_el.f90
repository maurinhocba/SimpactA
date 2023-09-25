      SUBROUTINE rest_el

      !  restores element databases at restart

      USE param_db,ONLY: midn, mnam
      USE lispa0, ONLY : lures
      USE esets_db,ONLY: msets, elty
      IMPLICIT NONE

      INTEGER (kind=4) :: i, nelem
      CHARACTER(len=midn) :: elmty
      CHARACTER(len=midn) :: heads
      CHARACTER(len=mnam) :: elsnam

      INTERFACE
        INCLUDE 'elmdat.h'
        INCLUDE 'elemnt.h'
      END INTERFACE

      DO
        READ (51) heads
        IF (heads == 'ENDSET' )EXIT

        READ (51) elmty, elsnam
        READ (51) nelem
        i = 1
        DO
          IF( TRIM(elmty) == TRIM(elty(i)) ) EXIT
          i = i+1
          IF( i <= msets)CYCLE
          WRITE(lures,*,ERR=9999)' Unrecognizable element type ',TRIM(elmty)
          CALL runen3 ('REST_EL: NON-VALID ELEMNT TYPE READ AT RESTART')
        END DO

        CALL elmdat('RESTAR',nelem,elsnam, i)

      END DO

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE rest_el


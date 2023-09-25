SUBROUTINE remesh(nn)
!-----------------------------------------------------------------------
! main remeshing routine
!-----------------------------------------------------------------------
USE lispa0, ONLY : lures
USE meshmo_db  !remeshing database
IMPLICIT NONE

  !--- Dummy variables
  CHARACTER(len=*),INTENT(IN) :: nn
  !--- Local variables
  LOGICAL          :: found
  INTEGER (kind=4) :: i, etype

  INTERFACE
    INCLUDE 'elemnt.h'
  END INTERFACE

  DO i=1,r_nrfset
    IF (.NOT.r_meshmo(i)) CYCLE
    CALL elemnt ('SEARCH', name=r_refset(i), igrav=etype, flag2=found)
    ! igrav - elemnt type
    IF (found) THEN
      WRITE (lures,'(/" Set ",A," will be remeshed.")',ERR=9999) TRIM(r_refset(i))
    ELSE
      WRITE (lures, '(" Set ",A," to be remeshed not found")',ERR=9999) TRIM(r_refset(i))
      CALL runen3 ('Remesh: Required set not found.')
    END IF
    IF ((etype == 20)) THEN
      ! remeshing for 2D solid elements
      CALL solrem(r_refset(i))
    ELSE
      CALL runen3('Remesh: Required set can not be remeshed.')
    END IF
  END DO
  WRITE(6,"(5X,'***  END OF REMESH NUMBER ',A,/)",err=9999) TRIM(nn)
  RETURN
 9999 CALL runen2('')
END SUBROUTINE remesh

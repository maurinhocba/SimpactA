SUBROUTINE inc_vcha(vcha,ninc,dflt)
!-----------------------------------------------------------------------
! Resize a vector of characters maintaining old data
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: ninc  !Number of new elements
  CHARACTER(len=*),INTENT(IN),OPTIONAL:: dflt   !Initialize default value
  CHARACTER(len=*),POINTER:: vcha(:)  !Vector to resize
  !--- Local variables
  INTEGER(kind=4):: ndim, nn
  CHARACTER(len=LEN(vcha)):: dft
  CHARACTER(len=LEN(vcha)),ALLOCATABLE:: work(:)

  nn = ninc
  IF (PRESENT(dflt)) THEN
    dft = dflt
  ELSE
    dft = ''
  END IF

  IF (ASSOCIATED(vcha)) THEN
    ndim = SIZE(vcha)      !Determine vector dimension
    ALLOCATE(work(ndim))   !Allocate auxiliar vector

    work(1:ndim) = vcha(1:ndim)   !Store old data in the auxiliar vector

    DEALLOCATE(vcha)              !Vector resize
    ALLOCATE(vcha(ndim+nn))       !

    IF (nn > 0) vcha(ndim+1:ndim+nn)=dft    !Initializes
    ndim = MIN0(ndim,ndim+nn)
    vcha(1:ndim) = work(1:ndim)   !Restore old data in the vector

    DEALLOCATE(work)       !Deallocate auxiliar vector

  ELSE
    IF (nn == 0) nn=1
    ALLOCATE(vcha(nn))     !Allocate a new vector
    vcha(1:nn) = dft       !Initializes

  END IF

RETURN
END SUBROUTINE inc_vcha

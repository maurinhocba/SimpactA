SUBROUTINE inc_vint(vint,ninc,dflt)
!-----------------------------------------------------------------------
! Resize a vector of integers maintaining old data
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: ninc  !Number of new elements
  INTEGER(kind=4),INTENT(IN),OPTIONAL:: dflt   !Initialize default value
  INTEGER(kind=4),POINTER:: vint(:)  !Vector to resize
  !--- Local variables
  INTEGER(kind=4):: ndim, nn, dft
  INTEGER(kind=4),ALLOCATABLE:: work(:)

  nn = ninc
  IF (PRESENT(dflt)) THEN
    dft = dflt
  ELSE
    dft = 0
  END IF

  IF (ASSOCIATED(vint)) THEN
    ndim = SIZE(vint)      !Determine vector dimension
    ALLOCATE(work(ndim))   !Allocate auxiliar vector

    work(1:ndim) = vint(1:ndim)   !Store old data in the auxiliar vector

    DEALLOCATE(vint)              !Vector resize
    ALLOCATE(vint(ndim+nn))       !

    IF (nn > 0) vint(ndim+1:ndim+nn)=dft    !Initializes
    ndim = MIN0(ndim,ndim+nn)
    vint(1:ndim) = work(1:ndim)   !Restore old data in the vector

    DEALLOCATE(work)       !Deallocate auxiliar vector

  ELSE
    IF (nn == 0) nn=1
    ALLOCATE(vint(nn))     !Allocate a new vector
    vint(1:nn) = dft       !Initializes

  END IF

RETURN
END SUBROUTINE inc_vint

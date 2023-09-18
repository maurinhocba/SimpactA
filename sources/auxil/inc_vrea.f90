SUBROUTINE inc_vrea(vrea,ninc,dflt)
!-----------------------------------------------------------------------
! Resize a vector of reals maintaining old data
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: ninc  !Number of new elements
  REAL(kind=8),INTENT(IN),OPTIONAL:: dflt   !Initialize default value
  REAL(kind=8),POINTER:: vrea(:)  !Vector to resize
  !--- Local variables
  INTEGER(kind=4):: ndim, nn
  REAL(kind=8):: dft
  REAL(kind=8),ALLOCATABLE:: work(:)

  nn = ninc
  IF (PRESENT(dflt)) THEN
    dft = dflt
  ELSE
    dft = 0
  END IF

  IF (ASSOCIATED(vrea)) THEN
    ndim = SIZE(vrea)      !Determine vector dimension
    ALLOCATE(work(ndim))   !Allocate auxiliar vector

    work(1:ndim) = vrea(1:ndim)   !Store old data in the auxiliar vector

    DEALLOCATE(vrea)              !Vector resize
    ALLOCATE(vrea(ndim+nn))       !

    IF (nn > 0) vrea(ndim+1:ndim+nn)=dft    !Initializes
    ndim = MIN0(ndim,ndim+nn)
    vrea(1:ndim) = work(1:ndim)   !Restore old data in the vector

    DEALLOCATE(work)       !Deallocate auxiliar vector

  ELSE
    IF (nn == 0) nn=1
    ALLOCATE(vrea(nn))     !Allocate a new vector
    vrea(1:nn) = dft       !Initializes

  END IF

RETURN
END SUBROUTINE inc_vrea

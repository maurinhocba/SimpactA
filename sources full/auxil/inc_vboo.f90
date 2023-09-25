SUBROUTINE inc_vboo(vboo,ninc,dflt)
!-----------------------------------------------------------------------
! Resize a vector of integers maintaining old data
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: ninc  !Number of new elements
  LOGICAL,INTENT(IN),OPTIONAL:: dflt   !Initialize default value
  LOGICAL,POINTER:: vboo(:)  !Vector to resize
  !--- Local variables
  INTEGER(kind=4):: ndim, nn
  LOGICAL:: dft
  LOGICAL,ALLOCATABLE:: work(:)

  nn = ninc
  IF (PRESENT(dflt)) THEN
    dft = dflt
  ELSE
    dft = .FALSE.
  END IF

  IF (ASSOCIATED(vboo)) THEN
    ndim = SIZE(vboo)      !Determine vector dimension
    ALLOCATE(work(ndim))   !Allocate auxiliar vector

    work(1:ndim) = vboo(1:ndim)   !Store old data in the auxiliar vector

    DEALLOCATE(vboo)              !Vector resize
    ALLOCATE(vboo(ndim+nn))       !

    IF (nn > 0) vboo(ndim+1:ndim+nn)=dft  !Initializes
    ndim = MIN0(ndim,ndim+nn)
    vboo(1:ndim) = work(1:ndim)     !Restore old data in the vector

    DEALLOCATE(work)     !Deallocate auxiliar vector

  ELSE
    IF (nn == 0) nn=1
    ALLOCATE(vboo(nn))   !Allocate a new vector
    vboo(1:nn) = dft     !Initializes

  END IF

RETURN
END SUBROUTINE inc_vboo

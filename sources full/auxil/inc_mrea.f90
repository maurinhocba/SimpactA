SUBROUTINE inc_mrea(mrea,nfil,ncol,dflt)
!-----------------------------------------------------------------------
! Resize an array of integers maintaining old data
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: nfil,   & !Number of new raws
                               ncol      !Number of new columns
  REAL(kind=8),INTENT(IN),OPTIONAL:: dflt    !Initialize default value
  REAL(kind=8),POINTER:: mrea(:,:)  !Array to resize
  !--- Local variables
  INTEGER(kind=4):: ndf, ndc, nf, nc
  REAL(kind=8):: dft
  REAL(kind=8),ALLOCATABLE:: work(:,:)

  nf = nfil
  nc = ncol

  IF (PRESENT(dflt)) THEN
    dft = dflt
  ELSE
    dft = 0d0
  END IF

  IF (ASSOCIATED(mrea)) THEN
    ndf = SIZE(mrea,DIM=1)      !Determine number of raws
    ndc = SIZE(mrea,DIM=2)      !Determine number of columns
    ALLOCATE(work(ndf,ndc))     !Allocate auxiliar array

    work(1:ndf,1:ndc) = mrea(1:ndf,1:ndc)   !Store old data in the auxiliar vector

    DEALLOCATE(mrea)                 !Array resize
    ALLOCATE(mrea(ndf+nf,ndc+nc))    !

    IF (nf > 0) mrea(ndf+1:ndf+nf,1:ndc)=dft            !Initializes
    IF (nc > 0) THEN
      mrea(1:ndf,ndc+1:ndc+nc) = dft                    !Initializes
      IF (nf > 0) mrea(ndf+1:ndf+nf,ndc+1:ndc+nc)=dft   !Initializes
    END IF

    ndf = MIN0(ndf,ndf+nf)
    ndc = MIN0(ndc,ndc+nc)
    mrea(1:ndf,1:ndc) = work(1:ndf,1:ndc)   !Restore old data in the vector

    DEALLOCATE(work)       !Deallocate auxiliar vector

  ELSE
    IF (nf == 0) nf=1
    IF (nc == 0) nc=1
    ALLOCATE(mrea(nf,nc))     !Allocate a new array
    mrea(1:nf,1:nc) = dft     !Initializes

  END IF

RETURN
END SUBROUTINE inc_mrea

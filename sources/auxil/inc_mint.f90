SUBROUTINE inc_mint(mint,nfil,ncol,dflt)
!-----------------------------------------------------------------------
! Resize an array of integers maintaining old data
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: nfil,   & !Number of new raws
                               ncol      !Number of new columns
  INTEGER(kind=4),INTENT(IN),OPTIONAL:: dflt    !Initialize default value
  INTEGER(kind=4),POINTER:: mint(:,:)  !Array to resize
  !--- Local variables
  INTEGER(kind=4):: ndf, ndc, nf, nc, dft
  INTEGER(kind=4),ALLOCATABLE:: work(:,:)

  nf = nfil
  nc = ncol

  IF (PRESENT(dflt)) THEN
    dft = dflt
  ELSE
    dft = 0
  END IF

  IF (ASSOCIATED(mint)) THEN
    ndf = SIZE(mint,DIM=1)      !Determine number of rows
    ndc = SIZE(mint,DIM=2)      !Determine number of columns
    ALLOCATE(work(ndf,ndc))     !Allocate auxiliar array

    work(1:ndf,1:ndc) = mint(1:ndf,1:ndc)   !Store old data in the auxiliar vector

    DEALLOCATE(mint)                 !Array resize
    ALLOCATE(mint(ndf+nf,ndc+nc))    !

    IF (nf > 0) mint(ndf+1:ndf+nf,1:ndc)=dft            !Initializes
    IF (nc > 0) THEN
      mint(1:ndf,ndc+1:ndc+nc) = dft                    !Initializes
      IF (nf > 0) mint(ndf+1:ndf+nf,ndc+1:ndc+nc)=dft   !Initializes
    END IF

    ndf = MIN0(ndf,ndf+nf)
    ndc = MIN0(ndc,ndc+nc)
    mint(1:ndf,1:ndc) = work(1:ndf,1:ndc)   !Restore old data in the vector

    DEALLOCATE(work)       !Deallocate auxiliar vector

  ELSE
    IF (nf == 0) nf=1
    IF (nc == 0) nc=1
    ALLOCATE(mint(nf,nc))     !Allocate a new array
    mint(1:nf,1:nc) = dft     !Initializes

  END IF

RETURN
END SUBROUTINE inc_mint

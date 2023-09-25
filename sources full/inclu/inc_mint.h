SUBROUTINE inc_mint(mint,nfil,ncol,dflt)
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN):: nfil,   & !Number of new raws
                               ncol      !Number of new columns
  INTEGER(kind=4),INTENT(IN),OPTIONAL:: dflt    !Initialize default value
  INTEGER(kind=4),POINTER:: mint(:,:)  !Array to resize

END SUBROUTINE inc_mint

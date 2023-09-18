SUBROUTINE inc_mrea(mrea,nfil,ncol,dflt)
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN):: nfil,   &
                               ncol
  REAL(kind=8),INTENT(IN),OPTIONAL:: dflt
  REAL(kind=8),POINTER:: mrea(:,:)

END SUBROUTINE inc_mrea

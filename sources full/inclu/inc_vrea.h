SUBROUTINE inc_vrea(vrea,ninc,dflt)
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN),OPTIONAL:: ninc
  REAL(kind=8),INTENT(IN),OPTIONAL:: dflt
  REAL(kind=8),POINTER:: vrea(:)

END SUBROUTINE inc_vrea

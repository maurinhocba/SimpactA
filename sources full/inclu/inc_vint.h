SUBROUTINE inc_vint(vint,ninc,dflt)
IMPLICIT NONE

  INTEGER(kind=4),INTENT(IN),OPTIONAL:: ninc
  INTEGER(kind=4),INTENT(IN),OPTIONAL:: dflt
  INTEGER(kind=4),POINTER:: vint(:)

END SUBROUTINE inc_vint

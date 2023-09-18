SUBROUTINE vecuni(n,v,modul)
!*****************************************************************************
! This routine compute de length of the vector v and converts it to a unit one
!*****************************************************************************
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: n     !Vector dimension
  REAL(kind=8),INTENT(INOUT):: v(n)  !Vector
  REAL(kind=8),INTENT(OUT):: modul   !Module of input vector

  modul = SQRT(DOT_PRODUCT(v(1:n),v(1:n)))
  IF (modul > 0d0) v(1:n)=v(1:n)/modul

RETURN
END SUBROUTINE vecuni

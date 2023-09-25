SUBROUTINE get_ifxb(nelbnd,lnbnd,ifx_b)
! get ifx codes and coordinates for boundary nodes in a remeshed set
IMPLICIT NONE

  !--- Dummy variables (INPUT)
  INTEGER(kind=4):: nelbnd               !number of segments defining contour
  INTEGER(kind=4),POINTER:: lnbnd(:)     !Conectivities of boundary segments
  !--- Dummy variables (OUTPUT)
  INTEGER(kind=4),POINTER:: ifx_b(:,:)   !(nd1,nelem)   boundary conditions

END SUBROUTINE get_ifxb

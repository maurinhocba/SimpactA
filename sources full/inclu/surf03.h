SUBROUTINE surf03 ( lnods, nnode, nelem)
!******************************************************************
! Get the boundary definition from the element set segments are
! oriented (outward normal)
!******************************************************************
!USE surf_db
IMPLICIT NONE

  !--- Dummy arguments
  INTEGER(kind=4),INTENT(IN):: nelem      !number of elements
  INTEGER(kind=4),INTENT(IN):: nnode      !number of nodes/element
  INTEGER(kind=4),INTENT(IN):: lnods(:,:) !(nnode,nelem) connectivities
  !TYPE(cont_srf),POINTER:: surfa   !INTENT(OUT) surface data
END SUBROUTINE surf03

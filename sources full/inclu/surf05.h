      SUBROUTINE surf05 ( lnods, nelem, nnode)
!******************************************************************
!
!*** get the boundary definition from the element set
!    segments are oriented (outward normal)
!
!******************************************************************
      USE surf_db
      IMPLICIT NONE
      !dummy arguments
      INTEGER (kind=4), INTENT(IN) :: nelem,nnode      !no of elements & nodes
      INTEGER (kind=4), INTENT(IN) :: lnods(:,:) !(nnode,nelem) connectivities
      !TYPE (cont_srf), POINTER :: surfa   !INTENT(OUT) surface data

      END SUBROUTINE surf05

      SUBROUTINE rdsegs (nnode, surface)

      !reads surface definitions

      USE surf_db
      IMPLICIT NONE
      INTEGER(kind=4) :: nelem,nnode
      TYPE (cont_srf), POINTER :: surface
      END SUBROUTINE rdsegs

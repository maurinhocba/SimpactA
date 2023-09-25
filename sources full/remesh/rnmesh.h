SUBROUTINE rnmesh(scalef,f_axsym)   !=> remesh_db   coorn(ndime,numnp)  l_new(nnode,ne_new)
! READ new mesh: nodal coordinates and element connectivities
IMPLICIT NONE

  !--- Dummy variables
  REAL(kind=8),INTENT(IN):: scalef
  LOGICAL,OPTIONAL:: f_axsym  ! TRUE => axisym. problem, introduced to correct negative x

END SUBROUTINE rnmesh

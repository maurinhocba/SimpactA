SUBROUTINE rdvels (iwrit, headn, tailn, nrve, ndofn)

  !reads a set of prescribed velocities

  USE rve_db
  IMPLICIT NONE
  INTEGER (kind=4) :: iwrit,nrve,ndofn
  TYPE (rve_nod), POINTER :: headn, tailn

END SUBROUTINE rdvels

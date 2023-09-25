SUBROUTINE rdpret (iwrit, headn, tailn, nrpt)

  !reads a set of prescribed temperatures

  USE c_input
  USE nsets_db
  USE ift_db
  IMPLICIT NONE
  INTEGER (kind=4) :: iwrit,nrpt
  TYPE (rpt_nod), POINTER :: headn, tailn

END SUBROUTINE rdpret

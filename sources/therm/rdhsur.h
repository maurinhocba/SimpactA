 SUBROUTINE rdhsur(nsurf,iwrit,heads,tails)

 !     Reads surface distributed flows

 USE c_input
 USE heat_db
 IMPLICIT NONE
 INTEGER (kind=4) :: nsurf,iwrit
 TYPE (srf_nod),POINTER :: heads, tails

 END SUBROUTINE rdhsur

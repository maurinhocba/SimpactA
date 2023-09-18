 SUBROUTINE rdhedg(nedge,iwrit,heads,tails)

 !     Reads edge distributed flows (solid 2D elements)

 USE c_input
 USE heat_db
 IMPLICIT NONE
 INTEGER (kind=4) :: nedge,iwrit
 TYPE (edg_nod),POINTER :: heads, tails

 END SUBROUTINE rdhedg

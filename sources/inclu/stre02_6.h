 SUBROUTINE stre02_6 (mat,l,stres,j)

 ! computes the internal nodal forces  1D (truss elements)
 ! for rubbers

 USE mat_dba, ONLY : mater
 IMPLICIT NONE
 REAL (kind=8), INTENT(IN) :: l  !lambda
 REAL (kind=8), INTENT(IN OUT) :: stres(:)
 REAL (kind=8), INTENT(OUT) :: j
 TYPE (mater), POINTER :: mat

 END SUBROUTINE stre02_6

 SUBROUTINE reduce_mass(mass,maxav,neq,maxa)
 !
 !  reduce mass matrix size checking null values
 !
 IMPLICIT NONE
 INTEGER(kind=4), INTENT(IN) :: neq
 INTEGER(kind=4), INTENT(IN OUT) :: maxav(:)
 INTEGER(kind=4), INTENT(OUT) :: maxa
 REAL(kind=8), INTENT(IN OUT) :: mass(:)

 END SUBROUTINE reduce_mass

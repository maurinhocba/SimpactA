 SUBROUTINE lump_mass(mass,maxav,neq)
 !
 !  reduce mass matrix size checking null values
 !
 IMPLICIT NONE
 INTEGER(kind=4), INTENT(IN) :: neq
 INTEGER(kind=4), INTENT(IN OUT) :: maxav(:)
 REAL(kind=8), INTENT(IN OUT) :: mass(:)

 REAL(kind=8) :: val
 INTEGER(kind=4) :: i,j,pi,pj,h

 DO i=2,neq
   pi = maxav(i)
   h = maxav(i+1) - pi - 1
   DO j=1,h
     val = mass(pi+j)
     pj = maxav(i-j)
     mass(pi) = mass(pi) + val
     mass(pj) = mass(pj) + val
     mass(pi+j) = 0d0
   END DO
 END DO
 RETURN
 END SUBROUTINE lump_mass

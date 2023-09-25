 SUBROUTINE rubber2d(pr,lb,model,stre,r1,r2,vec)

 ! calculates stiffness and stresses for rubber and elastomeric
 ! foam materials

 IMPLICIT NONE
 ! DUMMY arguments
 REAL(kind=8), INTENT(IN):: pr(12)    !(12) material constants
 REAL(kind=8), INTENT(IN OUT):: lb(3) !(3)  Principal stretches Eigenvalues of U
 INTEGER(kind=4), INTENT(IN OUT) :: model   !Model
 REAL(kind=8), INTENT(OUT) :: stre(:)           !  2nd Piola-Kirchhoff stresses S11 S22 S12 S33
 REAL(kind=8), INTENT(IN):: r1,r2     ! first eigenvector in-plane components
 ! optional arguments
 REAL(kind=8), OPTIONAL, INTENT(OUT) :: vec(:)  !elasticity (stiffness) matrix in vector form
 END SUBROUTINE rubber2d

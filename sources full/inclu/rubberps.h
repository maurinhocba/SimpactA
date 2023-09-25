 SUBROUTINE rubberps(pr,lb,model,stre,vec,r1,r2)

 ! calculates stiffness and stresses for rubber and elastomeric
 ! foam materials

 IMPLICIT NONE
 ! DUMMY arguments
 REAL(kind=8), INTENT(IN):: pr(12)    !(12) material constants
 REAL(kind=8), INTENT(IN OUT):: lb(3) !(3)  Principal stretches Eigenvalues of U
 INTEGER(kind=4), INTENT(IN) :: model   !Model
                                !1:      Arruda-Boyce
                                !2:      Mooney-Rivlin
                                !3:      Neo-Hooke
                                !4:      Ogden (N=1)
                                !5:      Ogden (N=2)
                                !6:      Ogden (N=3)
                                !7:      Polynomial (N=1)
                                !8:      Polynomial (N=2)
                                !9:      Polynomial (N=3)
                                !10:     Reduced Polynomial (N=1)
                                !11:     Reduced Polynomial (N=2)
                                !12:     Reduced Polynomial (N=3)
                                !13:     Van der Waals (not implemented yet)
                                !14:     Yeoh
                                !15:     Hyperfoam (N=1)
                                !16:     Hyperfoam (N=2)
                                !17:     Hyperfoam (N=3)
 REAL(kind=8), INTENT(OUT) :: stre(:)           !  2nd Piola-Kirchhoff stresses S11 S22 S12
 ! optional arguments
 REAL(kind=8), OPTIONAL, INTENT(OUT) :: vec(:)  ! elasticity (stiffness) matrix in vector form
 REAL(kind=8), OPTIONAL, INTENT(IN):: r1,r2     ! first eigenvector in-plane components
 END SUBROUTINE rubberps

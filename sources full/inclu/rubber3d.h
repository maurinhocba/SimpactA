 SUBROUTINE rubber3d(pr,g,model,stre,mat)

 ! calculates stiffness and stresses for rubber
 ! and elastomeric foam materials (not yet)

 IMPLICIT NONE
 ! DUMMY arguments
 REAL(kind=8), INTENT(IN):: pr(12)    !(12) material constants
 REAL(kind=8), INTENT(IN OUT):: g(6)  !(6)  Metric tensor (symmetric)
 INTEGER(kind=4), INTENT(IN OUT) :: model   !Model
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
 ! optional arguments
 REAL(kind=8), OPTIONAL, INTENT(OUT) :: &
           stre(6), &   !  2nd Piola-Kirchhoff stresses S11 S22 S33 S12 S13 S23
           mat(21)      !elasticity (stiffness) matrix in vector form

 END SUBROUTINE rubber3d

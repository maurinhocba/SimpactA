 SUBROUTINE bmat16q(cd,bb,t)
 !*** B matrix for a 3-D  Prism
 !    Quadratic approach for in plane components
 !    computed at center of triangular faces
 !    Green Lagrange Strains (TLF)
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8), INTENT(IN) :: t(3,2,3,2), & !deformation gradient at mid-sides
                              cd(4,2,3,2)   !cartesyan derivatives at faces
 REAL (kind=8), INTENT(OUT) :: bb(3,6,3,2)  !B-matrix at each face
 ! local variables
 INTEGER (kind=4) :: i,j,k

 bb = 0d0    !initializes
 DO i=1,2     !for each face
   DO k=1,3         !for each mid side point
     DO j=1,3         !for each main node
       bb(1:3,j,1,i) = bb(1:3,j,1,i) + cd(j,1,k,i)*t(1:3,1,k,i)
       bb(1:3,j,2,i) = bb(1:3,j,2,i) + cd(j,2,k,i)*t(1:3,2,k,i)
       bb(1:3,j,3,i) = bb(1:3,j,3,i) + cd(j,2,k,i)*t(1:3,1,k,i)+cd(j,1,k,i)*t(1:3,2,k,i)
     END DO
     j = k+3          !extra node
     bb(1:3,j,1,i) = bb(1:3,j,1,i) + cd(4,1,k,i)*t(1:3,1,k,i)
     bb(1:3,j,2,i) = bb(1:3,j,2,i) + cd(4,2,k,i)*t(1:3,2,k,i)
     bb(1:3,j,3,i) = bb(1:3,j,3,i) + cd(4,2,k,i)*t(1:3,1,k,i)+cd(4,1,k,i)*t(1:3,2,k,i)
   END DO
 END DO
 bb = bb/3d0     !average

 RETURN
 END SUBROUTINE bmat16q

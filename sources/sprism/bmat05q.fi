 SUBROUTINE bmat05q(cd,bb,t,quad)
 !*** B matrix for a 3-D  Prism
 !    Quadratic approach for in plane components
 !    computed at center of triangular faces
 !    Green Lagrange Strains (TLF)
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8), INTENT(IN) :: t(:,:,:,:), & !deformation gradient at mid-sides
                              cd(:,:,:,:)   !cartesyan derivatives at faces
 REAL (kind=8), INTENT(OUT) :: bb(:,:,:,:)  !B-matrix at each face
 LOGICAL, INTENT(IN) :: quad
 ! local variables
 INTEGER (kind=4) :: i,j,k

 IF( quad )THEN
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
 ELSE     !standard Constant Strain Triangle
   DO i=1,2     !for each face
     DO j=1,3         !for each main node
       bb(1:3,j,1,i) =  cd(j,1,1,i)*t(1:3,1,1,i)
       bb(1:3,j,2,i) =  cd(j,2,1,i)*t(1:3,2,1,i)
       bb(1:3,j,3,i) =  cd(j,2,1,i)*t(1:3,1,1,i)+cd(j,1,1,i)*t(1:3,2,1,i)
     END DO
   END DO
 END IF

 RETURN
 END SUBROUTINE bmat05q

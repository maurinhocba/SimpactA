 SUBROUTINE bmat04q(cd,bb,t,bm,area)
 !*** B matrix for a 3-D  Solid Shell Brick
 !    computed at mid-point sides
 !    Green Lagrange Strains (TLF)
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8), INTENT(IN) :: t(:,:,:,:),  &   !deformation gradient at mid-sides
                              cd(:,:,:,:), &   !cartesyan derivatives at faces
                              area(:,:)        !relative jacobian at mid-side points
 REAL (kind=8), INTENT(OUT) :: bb(:,:,:,:), & !(3,4,3,2)B-matrix at each face
                               bm(:,:,:)      !(12,4,2) B12-matrix at each mid-point side
 ! local variables
 INTEGER (kind=4) :: i,j,k
 REAL(kind=8) ::  bl(3,4,3,4) !B-matrix at each mid-point side

 bb = 0d0    !initializes
 DO i=1,2     !for each face
   DO k=1,4         !for each mid side point
     DO j=1,4         !for each main node
       bl(1:3,j,1,k) =  cd(j,1,k,i)*t(1:3,1,k,i)
       bl(1:3,j,2,k) =  cd(j,2,k,i)*t(1:3,2,k,i)
       bl(1:3,j,3,k) = (cd(j,2,k,i)*t(1:3,1,k,i)+cd(j,1,k,i)*t(1:3,2,k,i))
     END DO
     bb(:,:,:,i) = bb(:,:,:,i) + bl(:,:,:,k)*area(k,i)
     bm(:,k,i)   = RESHAPE(bl(:,:,3,k),(/12/))
   END DO
   DO k=1,4         !for each mid side point (Stabilization)
     bm(:,k,i) = bm(:,k,i) - RESHAPE(bb(:,:,3,i),(/12/))
   END DO
 END DO

 RETURN
 END SUBROUTINE bmat04q

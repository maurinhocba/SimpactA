 SUBROUTINE bmat04(bmat,bq,sh,bs,b3,f)
 !
 !*** B matrix  for a 3-D 8-node Brick solid-shell element
 !    Green Lagrange Strains (TLF)
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8), INTENT(IN) :: bq(:,:,:,:),  & !in-plane matrix at faces
                              bs(:,:,:,:),  & !assumed transverse shear matrix at faces
                              b3(:,:),      & !transverse strain at center
                              sh(:),        & !shape functions in zeta direction
                              f               !EXP(2*alpha*zeta)

 REAL (kind=8), INTENT(OUT) :: bmat(:,:,:)    !B-matrix


 ! local variables
 INTEGER (kind=4) :: j,j4                !indexes

 ! in-plane components
 DO j=1,4   !for each main node in lower face
   j4= j+4  !         main node in upper face
   bmat(1:3,j ,1) = sh(1)*bq(:, j,1,1)
   bmat(1:3,j ,2) = sh(1)*bq(:, j,2,1)
   bmat(1:3,j ,4) = sh(1)*bq(:, j,3,1)
   bmat(1:3,j4,1) = sh(2)*bq(:, j,1,2)
   bmat(1:3,j4,2) = sh(2)*bq(:, j,2,2)
   bmat(1:3,j4,4) = sh(2)*bq(:, j,3,2)
 END DO

 ! Assumed shear strains
 bmat(:,1:nnode,5:6) = sh(1)*bs(:,:,:,1) + sh(2)*bs(:,:,:,2)
 ! EAS for transverse strain
 bmat(:,1:nnode,3) =  b3*f

 RETURN
 END SUBROUTINE bmat04

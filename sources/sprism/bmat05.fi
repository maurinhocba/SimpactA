 SUBROUTINE bmat05(bmat,bq,sh,bs,b3,f,quad)
 !
 !*** B matrix  for a 3-D 6-node Prism solid-shell element
 !    Green Lagrange Strains (TLF)
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8), INTENT(IN) :: bq(:,:,:,:),  & !in-plane matrix at faces
                              bs(:,:,:,:),  & !assumed transverse shear matrix at faces
                              b3(:,:),      & !transverse strain at center
                              sh(:),        & !shape functions in zeta direction
                              f               !EXP(2*alpha*zeta)
 LOGICAL, INTENT(IN) :: quad                  !.TRUE. if quad approach

 REAL (kind=8), INTENT(OUT) :: bmat(:,:,:)    !B-matrix


 ! local variables
 INTEGER (kind=4) :: i,j,i3,j3                !indexes

 ! in-plane components
 DO j=1,3   !for each main node in lower face
   j3= j+3  !         main node in upper face
   bmat(1:3,j ,1) = sh(1)*bq(:, j,1,1)
   bmat(1:3,j ,2) = sh(1)*bq(:, j,2,1)
   bmat(1:3,j ,4) = sh(1)*bq(:, j,3,1)
   bmat(1:3,j3,1) = sh(2)*bq(:, j,1,2)
   bmat(1:3,j3,2) = sh(2)*bq(:, j,2,2)
   bmat(1:3,j3,4) = sh(2)*bq(:, j,3,2)
   IF( quad ) THEN
     i = j+6  !        extra node in lower face
     i3= i+3  !        extra node in upper face
     bmat(1:3,i ,1) = sh(1)*bq(:,j3,1,1)
     bmat(1:3,i ,2) = sh(1)*bq(:,j3,2,1)
     bmat(1:3,i ,4) = sh(1)*bq(:,j3,3,1)
     bmat(1:3,i3,1) = sh(2)*bq(:,j3,1,2)
     bmat(1:3,i3,2) = sh(2)*bq(:,j3,2,2)
     bmat(1:3,i3,4) = sh(2)*bq(:,j3,3,2)
   END IF
 END DO

 ! Assumed shear strains
 bmat(:,1:nnb,5:6) = sh(1)*bs(:,:,:,1) + sh(2)*bs(:,:,:,2)
 ! EAS for transverse strain
 bmat(:,1:nnb,3) =  b3*f

 RETURN
 END SUBROUTINE bmat05

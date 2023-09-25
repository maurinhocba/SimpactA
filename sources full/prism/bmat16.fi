 SUBROUTINE bmat16(cartd,bmat,t,quad,bb,sh,shell,bs,jacin)
 !
 !*** B matrix  for a 3-D 6-node Prism
 !    Green Lagrange Strains (TLF)
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8), INTENT(IN) :: cartd(nnode,3), & !cartesyan derivatives
                              t(3,3),         & !deformation gradient
                              sh(2)             !shape functions in zeta direction
 REAL (kind=8), INTENT(OUT) :: bmat(:,:,:)      !B-matrix
 REAL (kind=8) :: bb(:,:,:,:)                   !in-plane quadratic approach matrix
 REAL (kind=8), OPTIONAL :: bs(:,:,:)                     !assumed transverse shear approach matrix
 REAL (kind=8), OPTIONAL :: jacin(:,:)                    !in-plane inverse jacobian (assumed transverse shear approach)
 LOGICAL :: quad                                !flag for quad approach
 LOGICAL :: shell                               !flag for transverse shear approach
 ! local variables
 INTEGER (kind=4) :: i,j,i3,j3
 REAL(kind=8) :: baux(3,3)

 ! in-plane components
 IF( quad )THEN    !for quad approach
   DO j=1,3   !for each main node in lower face
     j3= j+3  !         main node in upper face
     i = j+6  !        extra node in lower face
     i3= i+3  !        extra node in upper face
     bmat(1:3,j ,1) = sh(1)*bb(:, j,1,1)
     bmat(1:3,j3,1) = sh(2)*bb(:, j,1,2)
     bmat(1:3,j ,2) = sh(1)*bb(:, j,2,1)
     bmat(1:3,j3,2) = sh(2)*bb(:, j,2,2)
     bmat(1:3,j ,4) = sh(1)*bb(:, j,3,1)
     bmat(1:3,j3,4) = sh(2)*bb(:, j,3,2)
     bmat(1:3,i ,1) = sh(1)*bb(:,j3,1,1)
     bmat(1:3,i3,1) = sh(2)*bb(:,j3,1,2)
     bmat(1:3,i ,2) = sh(1)*bb(:,j3,2,1)
     bmat(1:3,i3,2) = sh(2)*bb(:,j3,2,2)
     bmat(1:3,i ,4) = sh(1)*bb(:,j3,3,1)
     bmat(1:3,i3,4) = sh(2)*bb(:,j3,3,2)
   END DO
 ELSE               !standard
   DO j=1,nnode
     bmat(1:3,j,1) =  cartd(j,1)*t(1:3,1)
     bmat(1:3,j,2) =  cartd(j,2)*t(1:3,2)
     bmat(1:3,j,4) =  cartd(j,2)*t(1:3,1)+cartd(j,1)*t(1:3,2)
   END DO
 END IF

 ! out-of-plane components
 DO j=1,nnode
   bmat(1:3,j,3) =  cartd(j,3)*t(1:3,3)
   IF( shell )THEN
     baux(:,3) = (bs(:,j,1) + bs(:,j,2) + bs(:,j,3) )/3d0
     baux(:,1) =  -baux(:,3) + bs(:,j,3)
     baux(:,2) =  +baux(:,3) - bs(:,j,2)
     bmat(1:3,j,5) = jacin(1,1)*baux(:,1)+ jacin(2,1)*baux(:,2)
     bmat(1:3,j,6) = jacin(1,2)*baux(:,1)+ jacin(2,2)*baux(:,2)
   ELSE !standard
     bmat(1:3,j,5) =  cartd(j,3)*t(1:3,1)+cartd(j,1)*t(1:3,3)
     bmat(1:3,j,6) =  cartd(j,3)*t(1:3,2)+cartd(j,2)*t(1:3,3)
   END IF
 END DO

 RETURN
 END SUBROUTINE bmat16

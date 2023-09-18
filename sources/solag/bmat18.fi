 SUBROUTINE bmat18(cartd,bmat,t,nnode,shell)
 !
 !*** B matrix  for a 3-D hexaedral
 !    Green Lagrange Strains (TLF)
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: nnode
 REAL (kind=8), INTENT(IN) :: cartd(nnode,3),t(3,3)
 REAL (kind=8), INTENT(OUT) :: bmat(3,nnode,6)
 LOGICAL, INTENT(IN) :: shell

 INTEGER (kind=4) :: j

 DO j=1,nnode
   bmat(1:3,j,1) =  cartd(j,1)*t(1:3,1)
   bmat(1:3,j,2) =  cartd(j,2)*t(1:3,2)
   bmat(1:3,j,3) =  cartd(j,3)*t(1:3,3)
   bmat(1:3,j,4) =  cartd(j,2)*t(1:3,1)+cartd(j,1)*t(1:3,2)
   IF( shell ) CYCLE
   bmat(1:3,j,5) =  cartd(j,3)*t(1:3,1)+cartd(j,1)*t(1:3,3)
   bmat(1:3,j,6) =  cartd(j,3)*t(1:3,2)+cartd(j,2)*t(1:3,3)
 END DO

 RETURN
 END SUBROUTINE bmat18

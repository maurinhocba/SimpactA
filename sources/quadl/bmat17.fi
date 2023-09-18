 SUBROUTINE bmat17(cartd,bmat,t,nnode,rf,shape,ntype)

 !*** membrane matrix  for a 2-D quadrilateral

 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: nnode,ntype
 REAL (kind=8), INTENT(IN) :: cartd(nnode,2),t(2,2),rf,shape(nnode)
 REAL (kind=8), INTENT(OUT) :: bmat(2,nnode,4)
 INTEGER (kind=4) :: j

 DO j=1,nnode
   bmat(1:2,j,1) =  cartd(j,1)*t(1:2,1)
   bmat(1:2,j,2) =  cartd(j,2)*t(1:2,2)
   bmat(1:2,j,3) =  cartd(j,2)*t(1:2,1)+cartd(j,1)*t(1:2,2)
   IF( ntype == 3 ) bmat(1,j,4) = rf*shape(j)
 END DO

 RETURN
 END SUBROUTINE bmat17

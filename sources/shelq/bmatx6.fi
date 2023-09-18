 SUBROUTINE bmatx6(nnode,cartd,lb,bmatx,dx,dt,ndofe,nangl)
 !******************************************************************************
 !
 !     computes the tangent deformation matrix for shell elements
 !
 !******************************************************************************
 IMPLICIT NONE
 !***  routine parameters
 INTEGER (kind=4), INTENT(IN) :: nnode,ndofe
 REAL (kind=8), INTENT(IN) :: cartd(:,:),lb(:,:,:),dx(:,:),dt(:,:)
 REAL (kind=8), INTENT(IN), POINTER :: nangl(:,:)
 REAL (kind=8), INTENT(OUT) :: bmatx(:,:,:)
 !***  local variables
 INTEGER (kind=4) :: n
 REAL    (kind=8) :: coef1,coef2,coef3,coef4

 DO n = 1,nnode
   bmatx(1:3,n,1) = cartd(n,1)*dx(:,1)
   bmatx(1:3,n,2) = cartd(n,2)*dx(:,2)
   bmatx(1:3,n,3) = cartd(n,1)*dx(:,2) + cartd(n,2)*dx(:,1)
   bmatx(1:3,n,4) = cartd(n,1)*dt(:,1)
   bmatx(1:3,n,5) = cartd(n,2)*dt(:,2)
   bmatx(1:3,n,6) = cartd(n,1)*dt(:,2) + cartd(n,2)*dt(:,1)
   coef1 = DOT_PRODUCT(dx(:,1),lb(:,1,n))
   coef2 = DOT_PRODUCT(dx(:,1),lb(:,2,n))
   coef3 = DOT_PRODUCT(dx(:,2),lb(:,1,n))
   coef4 = DOT_PRODUCT(dx(:,2),lb(:,2,n))
   bmatx(4,n,4) =  - cartd(n,1)*coef2
   bmatx(5,n,4) =  + cartd(n,1)*coef1
   bmatx(4,n,5) =  - cartd(n,2)*coef4
   bmatx(5,n,5) =  + cartd(n,2)*coef3
   bmatx(4,n,6) =  - cartd(n,1)*coef4 - cartd(n,2)*coef2
   bmatx(5,n,6) =  + cartd(n,1)*coef3 + cartd(n,2)*coef1
   IF( ndofe == 7 )THEN
     bmatx(6:7,n, 9) = cartd(n,1)*nangl(:,n)
     bmatx(6:7,n,10) = cartd(n,2)*(/-nangl(2,n),nangl(1,n)/)
     bmatx(6:7,n,11) = cartd(n,2)*nangl(:,n)
     bmatx(6:7,n,12) = cartd(n,1)*(/-nangl(2,n),nangl(1,n)/)
   END IF
 END DO
 RETURN
 END SUBROUTINE bmatx6

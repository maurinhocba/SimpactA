 SUBROUTINE bmatx7(stype,cartd,lb,bmatx,dx,dt,ndofe,nangl)
 !******************************************************************************
 !
 !***** this routine compute the tangent deformation matrix for element 7
 !      membrane and bending and shear for stype = 1
 !******************************************************************************
 IMPLICIT NONE
 !                      routine parameters

 INTEGER (kind=4) stype,ndofe

 REAL (kind=8), INTENT(IN) :: cartd(:,:),lb(:,:,:),dx(:,:),dt(:,:)
 REAL (kind=8), INTENT(IN), POINTER :: nangl(:,:)
 REAL (kind=8), INTENT(OUT) :: bmatx(:,:,:)

 !                      local variables

 INTEGER (kind=4) n

 REAL    (kind=8) c11,c12,c21,c22

 DO n=1,6
   IF(stype <= 1) THEN
     ! membrane contribution (mid-surface displacement)
     bmatx(1:3,n,1) =  cartd(n,1)*dx(:,1)
     bmatx(1:3,n,2) =  cartd(n,2)*dx(:,2)
     bmatx(1:3,n,3) =  cartd(n,1)*dx(:,2) + cartd(n,2)*dx(:,1)
   END IF

   ! bending contribution to mid-surface displacement
   bmatx(1:3,n,4) =  cartd(n,1)*dt(:,1)
   bmatx(1:3,n,5) =  cartd(n,2)*dt(:,2)
   bmatx(1:3,n,6) =  cartd(n,1)*dt(:,2) + cartd(n,2)*dt(:,1)

   c11 = DOT_PRODUCT(dx(:,1),lb(:,1,n))
   c12 = DOT_PRODUCT(dx(:,1),lb(:,2,n))
   c21 = DOT_PRODUCT(dx(:,2),lb(:,1,n))
   c22 = DOT_PRODUCT(dx(:,2),lb(:,2,n))
   ! bending contribution to directors
   bmatx(4,n,4) = - cartd(n,1)*c12
   bmatx(5,n,4) = + cartd(n,1)*c11
   bmatx(4,n,5) = - cartd(n,2)*c22
   bmatx(5,n,5) = + cartd(n,2)*c21
   bmatx(4,n,6) =(- cartd(n,1)*c22 - cartd(n,2)*c12)
   bmatx(5,n,6) =(+ cartd(n,1)*c21 + cartd(n,2)*c11)
   IF( ndofe == 7 )THEN
     bmatx(6:7,n, 9) = cartd(n,1)*nangl(:,n)
     bmatx(6:7,n,10) = cartd(n,2)*(/-nangl(2,n),nangl(1,n)/)
     bmatx(6:7,n,11) = cartd(n,2)*nangl(:,n)
     bmatx(6:7,n,12) = cartd(n,1)*(/-nangl(2,n),nangl(1,n)/)
   END IF

 END DO
 RETURN
 END SUBROUTINE bmatx7

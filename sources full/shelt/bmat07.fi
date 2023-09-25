 SUBROUTINE bmat07(cartd,lb,bb,bmem,dx,dt,quad,mems,tm,cd,ndofe,nangl,bf)
 !******************************************************************************
 !
 !***** this routine compute the tangent deformation matrix for element 03
 !      membrane and bending
 !******************************************************************************
 IMPLICIT NONE
 !                      routine parameters

 INTEGER(kind=4), INTENT(IN) :: ndofe
 REAL (kind=8), INTENT(IN) :: cartd(:,:),lb(:,:,:),dx(:,:),dt(:,:),tm(:,:,:)
 REAL (kind=8), POINTER :: cd(:,:,:)
 REAL (kind=8), INTENT(IN), POINTER :: nangl(:,:)
 LOGICAL, INTENT(IN) :: quad,mems(:)
 REAL (kind=8), INTENT(OUT) :: bb(:,:,:),bmem(:,:,:),bf(:,:,:)

 !                      local variables

 INTEGER (kind=4) n,i,j

 REAL    (kind=8) c11,c12,c21,c22

!membrane part
 IF(quad)THEN
   bmem = 0d0        !initializes
   DO i=1,3          !for each side
     IF( mems(i) )THEN  !linear strain approach
       !d1 = ABS(a0(i)+2d0*gamma(i))    !angle at boundary
       !IF( d1 <= alp1 )THEN            !smooth  use side (cuadratic) metric tensor
          DO j=1,3          !node
            bmem(:,j,1) = bmem(:,j,1) + cd(j,1,i)*tm(:,1,i)
            bmem(:,j,2) = bmem(:,j,2) + cd(j,2,i)*tm(:,2,i)
            bmem(:,j,3) = bmem(:,j,3) + cd(j,2,i)*tm(:,1,i) + cd(j,1,i)*tm(:,2,i)
          END DO
          j = i+3          ! extra node
          bmem(:,j,1) = bmem(:,j,1) + cd(4,1,i)*tm(:,1,i)
          bmem(:,j,2) = bmem(:,j,2) + cd(4,2,i)*tm(:,2,i)
          bmem(:,j,3) = bmem(:,j,3) + cd(4,2,i)*tm(:,1,i) + cd(4,1,i)*tm(:,2,i)

       !ELSE IF( d1 >= alp2 )THEN       !folded, use main triangle metric tensor
       !  DO j=1,3          !node
       !    bmem(:,j,1) = bmem(:,j,1) - b(j)*t(:,1)
       !    bmem(:,j,2) = bmem(:,j,2) + a(j)*t(:,2)
       !    bmem(:,j,3) = bmem(:,j,3) + a(j)*t(:,1) - b(j)*t(:,2)
       !  END DO

       !ELSE                            !interpolate
       !  d1 = (d1-alp1)/alp3           !function for main triangle
       !  d2 = 1d0 - d1                 !function for cuadratic interpolation
       !  DO j=1,3          !node
       !    bmem(:,j,1) = bmem(:,j,1)- d1*b(j)*t(:,1) + d2*cd(j,1,i)*t(:,n)
       !    bmem(:,j,2) = bmem(:,j,2)+ d1*a(j)*t(:,2) + d2*cd(j,2,i)*t(:,n1)
       !    bmem(:,j,3) = bmem(:,j,3)+ d1*(a(j)*t(:,1) - b(j)*t(:,2)) &
       !                             + d2*(cd(j,2,i)*t(:,n) + cd(j,1,i)*t(:,n1))
       !  END DO
       !  j = i+3          ! extra node
       !  bmem(:,j,1) = bmem(:,j,1) + d2*cd(4,1,i)*t(:,n)
       !  bmem(:,j,2) = bmem(:,j,2) + d2*cd(4,2,i)*t(:,n1)
       !  bmem(:,j,3) = bmem(:,j,3) + d2*(cd(4,2,i)*t(:,n) + cd(4,1,i)*t(:,n1))

       !END IF
     ELSE   ! constant strain triangle
       DO j=1,3          !node
         bmem(:,j,1) = bmem(:,j,1) + cartd(j,1)*dx(:,1)
         bmem(:,j,2) = bmem(:,j,2) + cartd(j,2)*dx(:,2)
         bmem(:,j,3) = bmem(:,j,3) + cartd(j,2)*dx(:,1) + cartd(j,1)*dx(:,2)
       END DO
     END IF
   END DO
   bmem =bmem/3d0
 ELSE   ! constant strain triangle
   DO j=1,3          !node
     bmem(:,j,1) = cartd(j,1)*dx(:,1)
     bmem(:,j,2) = cartd(j,2)*dx(:,2)
     bmem(:,j,3) = cartd(j,2)*dx(:,1) + cartd(j,1)*dx(:,2)
   END DO
 END IF

!bending part

 DO n=1,3

   ! bending contribution to mid-surface displacement
   bb(1:3,n,1) =  cartd(n,1)*dt(:,1)
   bb(1:3,n,2) =  cartd(n,2)*dt(:,2)
   bb(1:3,n,3) =  cartd(n,1)*dt(:,2) + cartd(n,2)*dt(:,1)

   c11 = DOT_PRODUCT(dx(:,1),lb(:,1,n))
   c12 = DOT_PRODUCT(dx(:,1),lb(:,2,n))
   c21 = DOT_PRODUCT(dx(:,2),lb(:,1,n))
   c22 = DOT_PRODUCT(dx(:,2),lb(:,2,n))
   ! bending contribution to directors
   bb(4,n,1) = +2d0*   cartd(n,1)*c12
   bb(5,n,1) = -2d0*   cartd(n,1)*c11
   bb(4,n,2) = +2d0*   cartd(n,2)*c22
   bb(5,n,2) = -2d0*   cartd(n,2)*c21
   bb(4,n,3) =  2d0*(+ cartd(n,1)*c22 + cartd(n,2)*c12)
   bb(5,n,3) =  2d0*(- cartd(n,1)*c21 - cartd(n,2)*c11)
   IF( ndofe == 7 )THEN
     bf(1:2,n, 1) = cartd(n,1)*nangl(:,n)
     bf(1:2,n, 2) = cartd(n,2)*(/-nangl(2,n),nangl(1,n)/)
     bf(1:2,n, 3) = cartd(n,2)*nangl(:,n)
     bf(1:2,n, 4) = cartd(n,1)*(/-nangl(2,n),nangl(1,n)/)
   END IF
 END DO
 RETURN
 END SUBROUTINE bmat07

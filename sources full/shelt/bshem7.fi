 SUBROUTINE bshem7(x,t,lb,bs,bf,jas,nangl,ngamm)
 !*****************************************************************************
 !
 !*****evaluates b matrix for shear at the assumed strain points
 !
 !****************************************************************************
 IMPLICIT NONE
 !                        routine parameters

 REAL (kind=8), INTENT(IN) :: x(:,:),t(:,:),lb(:,:,:)
 REAL (kind=8), INTENT(IN), POINTER :: jas(:,:),nangl(:,:)
 REAL (kind=8), INTENT(OUT) :: bs(:,:,:),bf(:,:,:)
 INTEGER (kind=4), INTENT(IN) :: ngamm

 !                        local variables

 REAL    (kind=8) tt(3),dx(3) !bc(2,6,6)
 INTEGER (kind=8) :: i,j

 DO i=1,ngamm
   tt = nf(1)*t(:,kk(1,i))+ nf(2)*t(:,kk(2,i))+ nf(3)*t(:,kk(3,i))
   dx = nd(1)*x(:,kk(1,i))+ nd(2)*x(:,kk(2,i))+ nd(3)*x(:,kk(3,i))
   IF( MOD(i,2) == 0 )THEN
     dx = -dx  !this change is because nd() changes the sign
     tt = -tt  !this change is because is lately multiplied by -nd()
   END IF
   DO j=1,3
     bs(1:3,kk(j,i),i) = nd(j)*tt
     bs(  4,kk(j,i),i) = -DOT_PRODUCT(dx,lb(:,2,kk(j,i)))*nf(j)
     bs(  5,kk(j,i),i) =  DOT_PRODUCT(dx,lb(:,1,kk(j,i)))*nf(j)
   END DO
 END DO

! ! A  first node side 1 - 4 - 2
! tt = fxa*t(:,1)+fxc*t(:,4)+fxb*t(:,2) ; dx =-fxd*x(:,1)+fxf*x(:,4)-fxe*x(:,2)
! bs(1:3,1,1) = -fxd*tt ; bs(1:3,4,1) =  fxf*tt ; bs(1:3,2,1) = -fxe*tt
! bs(4,1,1) = -DOT_PRODUCT(dx,lb(:,2,1))*fxa ; bs(5,1,1) =  DOT_PRODUCT(dx,lb(:,1,1))*fxa
! bs(4,4,1) = -DOT_PRODUCT(dx,lb(:,2,4))*fxc ; bs(5,4,1) =  DOT_PRODUCT(dx,lb(:,1,4))*fxc
! bs(4,2,1) = -DOT_PRODUCT(dx,lb(:,2,2))*fxb ; bs(5,2,1) =  DOT_PRODUCT(dx,lb(:,1,2))*fxb
! ! B  second node side 1 - 4 - 2
! tt = fxb*t(:,1)+fxc*t(:,4)+fxa*t(:,2) ; dx = fxe*x(:,1)-fxf*x(:,4)+fxd*x(:,2)
! bs(1:3,1,2) =  fxe*tt ; bs(1:3,4,2) = -fxf*tt ; bs(1:3,2,2) =  fxd*tt
! bs(4,1,2) = -DOT_PRODUCT(dx,lb(:,2,1))*fxb ; bs(5,1,2) =  DOT_PRODUCT(dx,lb(:,1,1))*fxb
! bs(4,4,2) = -DOT_PRODUCT(dx,lb(:,2,4))*fxc ; bs(5,4,2) =  DOT_PRODUCT(dx,lb(:,1,4))*fxc
! bs(4,2,2) = -DOT_PRODUCT(dx,lb(:,2,2))*fxa ; bs(5,2,2) =  DOT_PRODUCT(dx,lb(:,1,2))*fxa
! ! C  first node side 2 - 5 - 3
! tt = fxa*t(:,2)+fxc*t(:,5)+fxb*t(:,3) ; dx =-fxd*x(:,2)+fxf*x(:,5)-fxe*x(:,3)
! bs(1:3,2,3) = -fxd*tt ; bs(1:3,5,3) =  fxf*tt ; bs(1:3,3,3) = -fxe*tt
! bs(4,2,3) = -DOT_PRODUCT(dx,lb(:,2,2))*fxa ; bs(5,2,3) =  DOT_PRODUCT(dx,lb(:,1,2))*fxa
! bs(4,5,3) = -DOT_PRODUCT(dx,lb(:,2,5))*fxc ; bs(5,5,3) =  DOT_PRODUCT(dx,lb(:,1,5))*fxc
! bs(4,3,3) = -DOT_PRODUCT(dx,lb(:,2,3))*fxb ; bs(5,3,3) =  DOT_PRODUCT(dx,lb(:,1,3))*fxb
! ! D  second node side 2 - 5 - 3
! tt = fxb*t(:,2)+fxc*t(:,5)+fxa*t(:,3) ; dx = fxe*x(:,2)-fxf*x(:,5)+fxd*x(:,3)
! bs(1:3,2,4) =  fxe*tt ; bs(1:3,5,4) = -fxf*tt ; bs(1:3,3,4) =  fxd*tt
! bs(4,2,4) = -DOT_PRODUCT(dx,lb(:,2,2))*fxb ; bs(5,2,4) =  DOT_PRODUCT(dx,lb(:,1,2))*fxb
! bs(4,5,4) = -DOT_PRODUCT(dx,lb(:,2,5))*fxc ; bs(5,5,4) =  DOT_PRODUCT(dx,lb(:,1,5))*fxc
! bs(4,3,4) = -DOT_PRODUCT(dx,lb(:,2,3))*fxa ; bs(5,3,4) =  DOT_PRODUCT(dx,lb(:,1,3))*fxa
! ! E  first node side 3 - 6 - 1
! tt = fxa*t(:,3)+fxc*t(:,6)+fxb*t(:,1) ; dx =-fxd*x(:,3)+fxf*x(:,6)-fxe*x(:,1)
! bs(1:3,3,5) = -fxd*tt ; bs(1:3,6,5) =  fxf*tt ; bs(1:3,1,5) = -fxe*tt
! bs(4,3,5) = -DOT_PRODUCT(dx,lb(:,2,3))*fxa ; bs(5,3,5) =  DOT_PRODUCT(dx,lb(:,1,3))*fxa
! bs(4,6,5) = -DOT_PRODUCT(dx,lb(:,2,6))*fxc ; bs(5,6,5) =  DOT_PRODUCT(dx,lb(:,1,6))*fxc
! bs(4,1,5) = -DOT_PRODUCT(dx,lb(:,2,1))*fxb ; bs(5,1,5) =  DOT_PRODUCT(dx,lb(:,1,1))*fxb
! ! F  second node side 3 - 6 - 1
! tt = fxb*t(:,3)+fxc*t(:,6)+fxa*t(:,1) ; dx = fxe*x(:,3)-fxf*x(:,6)+fxd*x(:,1)
! bs(1:3,3,6) =  fxe*tt ; bs(1:3,6,6) = -fxf*tt ; bs(1:3,1,6) =  fxd*tt
! bs(4,3,6) = -DOT_PRODUCT(dx,lb(:,2,3))*fxb ; bs(5,3,6) =  DOT_PRODUCT(dx,lb(:,1,3))*fxb
! bs(4,6,6) = -DOT_PRODUCT(dx,lb(:,2,6))*fxc ; bs(5,6,6) =  DOT_PRODUCT(dx,lb(:,1,6))*fxc
! bs(4,1,6) = -DOT_PRODUCT(dx,lb(:,2,1))*fxa ; bs(5,1,6) =  DOT_PRODUCT(dx,lb(:,1,1))*fxa

 IF( ASSOCIATED(jas) )THEN
   !bc = 0d0
   DO i=1,ngamm
     DO j=1,3
       bf(1,kk(j,i),i) = (nangl(1,kk(j,i))*jas(1,i)-nangl(2,kk(j,i))*jas(2,i))*nf(j)
       bf(2,kk(j,i),i) = (nangl(1,kk(j,i))*jas(2,i)+nangl(2,kk(j,i))*jas(1,i))*nf(j)
     END DO
   END DO
!   bf = 0d0
!   ! A  first node side 1 - 4 - 2
!   bf(1,1,1) = (nangl(1,1)*jas(1,1)-nangl(2,1)*jas(2,1))*fxa
!   bf(2,1,1) = (nangl(1,1)*jas(2,1)+nangl(2,1)*jas(1,1))*fxa
!   bf(1,4,1) = (nangl(1,4)*jas(1,1)-nangl(2,4)*jas(2,1))*fxc
!   bf(2,4,1) = (nangl(1,4)*jas(2,1)+nangl(2,4)*jas(1,1))*fxc
!   bf(1,2,1) = (nangl(1,2)*jas(1,1)-nangl(2,2)*jas(2,1))*fxb
!   bf(2,2,1) = (nangl(1,2)*jas(2,1)+nangl(2,2)*jas(1,1))*fxb
!   ! B  second node side 1 - 4 - 2
!   bf(1,1,2) = (nangl(1,1)*jas(1,2)-nangl(2,1)*jas(2,2))*fxb
!   bf(2,1,2) = (nangl(1,1)*jas(2,2)+nangl(2,1)*jas(1,2))*fxb
!   bf(1,4,2) = (nangl(1,4)*jas(1,2)-nangl(2,4)*jas(2,2))*fxc
!   bf(2,4,2) = (nangl(1,4)*jas(2,2)+nangl(2,4)*jas(1,2))*fxc
!   bf(1,2,2) = (nangl(1,2)*jas(1,2)-nangl(2,2)*jas(2,2))*fxa
!   bf(2,2,2) = (nangl(1,2)*jas(2,2)+nangl(2,2)*jas(1,2))*fxa
!   ! C  first node side 2 - 5 - 3
!   bf(1,2,3) = (nangl(1,2)*jas(1,3)-nangl(2,2)*jas(2,3))*fxa
!   bf(2,2,3) = (nangl(1,2)*jas(2,3)+nangl(2,2)*jas(1,3))*fxa
!   bf(1,5,3) = (nangl(1,5)*jas(1,3)-nangl(2,5)*jas(2,3))*fxc
!   bf(2,5,3) = (nangl(1,5)*jas(2,3)+nangl(2,5)*jas(1,3))*fxc
!   bf(1,3,3) = (nangl(1,3)*jas(1,3)-nangl(2,3)*jas(2,3))*fxb
!   bf(2,3,3) = (nangl(1,3)*jas(2,3)+nangl(2,3)*jas(1,3))*fxb
!   ! D  second node side 2 - 5 - 3
!   bf(1,2,4) = (nangl(1,2)*jas(1,4)-nangl(2,2)*jas(2,4))*fxb
!   bf(2,2,4) = (nangl(1,2)*jas(2,4)+nangl(2,2)*jas(1,4))*fxb
!   bf(1,5,4) = (nangl(1,5)*jas(1,4)-nangl(2,5)*jas(2,4))*fxc
!   bf(2,5,4) = (nangl(1,5)*jas(2,4)+nangl(2,5)*jas(1,4))*fxc
!   bf(1,3,4) = (nangl(1,3)*jas(1,4)-nangl(2,3)*jas(2,4))*fxa
!   bf(2,3,4) = (nangl(1,3)*jas(2,4)+nangl(2,3)*jas(1,4))*fxa
!   ! E  first node side 3 - 6 - 1
!   bf(1,3,5) = (nangl(1,3)*jas(1,5)-nangl(2,3)*jas(2,5))*fxa
!   bf(2,3,5) = (nangl(1,3)*jas(2,5)+nangl(2,3)*jas(1,5))*fxa
!   bf(1,6,5) = (nangl(1,6)*jas(1,5)-nangl(2,6)*jas(2,5))*fxc
!   bf(2,6,5) = (nangl(1,6)*jas(2,5)+nangl(2,6)*jas(1,5))*fxc
!   bf(1,1,5) = (nangl(1,1)*jas(1,5)-nangl(2,1)*jas(2,5))*fxb
!   bf(2,1,5) = (nangl(1,1)*jas(2,5)+nangl(2,1)*jas(1,5))*fxb
!   ! F  second node side 3 - 6 - 1
!   bf(1,3,6) = (nangl(1,3)*jas(1,6)-nangl(2,3)*jas(2,6))*fxb
!   bf(2,3,6) = (nangl(1,3)*jas(2,6)+nangl(2,3)*jas(1,6))*fxb
!   bf(1,6,6) = (nangl(1,6)*jas(1,6)-nangl(2,6)*jas(2,6))*fxc
!   bf(2,6,6) = (nangl(1,6)*jas(2,6)+nangl(2,6)*jas(1,6))*fxc
!   bf(1,1,6) = (nangl(1,1)*jas(1,6)-nangl(2,1)*jas(2,6))*fxa
!   bf(2,1,6) = (nangl(1,1)*jas(2,6)+nangl(2,1)*jas(1,6))*fxa
!   write(58,"(4e15.6,2f10.5)")((bf(:,i,j),bc(:,i,j),bf(:,i,j)-bc(:,i,j),i=1,6),j=1,6)
!   bc = bc-bf
 END IF
 RETURN
 END SUBROUTINE bshem7

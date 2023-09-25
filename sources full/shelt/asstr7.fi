 SUBROUTINE asstr7(x,t,gamma,psi,psib,jas,ngamm)
 !*****************************************************************************
 !
 !***** natural coordinates shear deformation at assumed strain points
 !
 !*****************************************************************************
 IMPLICIT NONE
 !                        routine parameters


 INTEGER(kind=4), INTENT(IN) :: ngamm
 REAL (kind=8), INTENT(IN) :: x(:,:),t(:,:),psi(:,:)
 REAL (kind=8), INTENT(OUT) :: gamma(:),psib(:)
 REAL (kind=8), POINTER :: jas(:,:)

 !                        local variables

 REAL (kind=8) dx(3),tt(3),p(2) !,gammb(6)
 INTEGER (kind=8) :: i

 DO i=1,ngamm
   tt = nf(1)*t(:,kk(1,i))+ nf(2)*t(:,kk(2,i))+ nf(3)*t(:,kk(3,i))
   dx = nd(1)*x(:,kk(1,i))+ nd(2)*x(:,kk(2,i))+ nd(3)*x(:,kk(3,i))
   gamma(i) = DOT_PRODUCT(tt,dx)
   IF( MOD(i,2) == 0 )gamma(i) = -gamma(i)
 END DO

! ! A  first node side 1 - 4 - 2
! tt = fxa*t(:,1)+fxc*t(:,4)+fxb*t(:,2) ; dx =-fxd*x(:,1)+fxf*x(:,4)-fxe*x(:,2)
! gamma(1) = DOT_PRODUCT(dx,tt)
! ! B  second node side 1 - 4 - 2
! tt = fxb*t(:,1)+fxc*t(:,4)+fxa*t(:,2) ; dx = fxe*x(:,1)-fxf*x(:,4)+fxd*x(:,2)
! gamma(2) = DOT_PRODUCT(dx,tt)
! ! C  first node side 2 - 5 - 3
! tt = fxa*t(:,2)+fxc*t(:,5)+fxb*t(:,3) ; dx =-fxd*x(:,2)+fxf*x(:,5)-fxe*x(:,3)
! gamma(3) = DOT_PRODUCT(dx,tt)
! ! D  second node side 2 - 5 - 3
! tt = fxb*t(:,2)+fxc*t(:,5)+fxa*t(:,3) ; dx = fxe*x(:,2)-fxf*x(:,5)+fxd*x(:,3)
! gamma(4) = DOT_PRODUCT(dx,tt)
! ! E  first node side 3 - 6 - 1
! tt = fxa*t(:,3)+fxc*t(:,6)+fxb*t(:,1) ; dx =-fxd*x(:,3)+fxf*x(:,6)-fxe*x(:,1)
! gamma(5) = DOT_PRODUCT(dx,tt)
! ! F  second node side 3 - 6 - 1
! tt = fxb*t(:,3)+fxc*t(:,6)+fxa*t(:,1) ; dx = fxe*x(:,3)-fxf*x(:,6)+fxd*x(:,1)
! gamma(6) = DOT_PRODUCT(dx,tt)

 IF ( ASSOCIATED(jas) )THEN
   DO i=1,ngamm
     p = nf(1)*psi(:,kk(1,i))+nf(2)*psi(:,kk(1,i))+nf(3)*psi(:,kk(1,i))
     psib(i) = DOT_PRODUCT(jas(:,i),p)
   END DO
   !! A  first node side 1 - 4 - 2
   !p = fxa*psi(:,1)+fxc*psi(:,4)+fxb*psi(:,2)
   !psib(1) = DOT_PRODUCT(jas(:,1),p)
   !! B  second node side 1 - 4 - 2
   !p = fxb*psi(:,1)+fxc*psi(:,4)+fxa*psi(:,2)
   !psib(2) = DOT_PRODUCT(jas(:,2),p)
   !! C  first node side 2 - 5 - 3
   !p = fxa*psi(:,2)+fxc*psi(:,5)+fxb*psi(:,3)
   !psib(3) = DOT_PRODUCT(jas(:,3),p)
   !! D  second node side 2 - 5 - 3
   !p = fxb*psi(:,2)+fxc*psi(:,5)+fxa*psi(:,3)
   !psib(4) = DOT_PRODUCT(jas(:,4),p)
   !! E  first node side 3 - 6 - 1
   !p = fxa*psi(:,3)+fxc*psi(:,6)+fxb*psi(:,1)
   !psib(5) = DOT_PRODUCT(jas(:,5),p)
   !! F  second node side 3 - 6 - 1
   !p = fxb*psi(:,3)+fxc*psi(:,6)+fxa*psi(:,1)
   !psib(6) = DOT_PRODUCT(jas(:,6),p)
   !gammb = gammb - psib
 END IF

 RETURN
 END SUBROUTINE asstr7

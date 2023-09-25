 SUBROUTINE stra07(cartd, x,t,stran,dx,dt,gamma,quad,sides,tm,cd,psi,psib,jas)
 !***********************************************************************
 !
 !*****this routine computes 1st and 2nd fundamental forms and distorsions
 !            for element 3
 !***********************************************************************
 IMPLICIT NONE
 !              routine parameters

 REAL (kind=8), INTENT(IN) :: cartd(:,:),x(:,:),t(:,:)
 REAL (kind=8), INTENT(OUT) :: stran(:),dx(:,:),dt(:,:),gamma(:),tm(:,:,:)
 LOGICAL, INTENT(IN) :: quad,sides(:)
 REAL (kind=8), POINTER :: cd(:,:,:)
 REAL (kind=8), INTENT(IN), OPTIONAL :: psi(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL :: psib(:)
 REAL (kind=8), INTENT(IN), POINTER, OPTIONAL :: jas(:,:)
!              local variables

 INTEGER (kind=4) i

 !     evaluates director & configuration & director derivatives

 dx(:,:) = MATMUL(x(:,1:3),cartd(:,:))
 dt(:,:) = -2d0*MATMUL(t(:,1:3),cartd(:,:))

 !   correct director derivatives (may be unnecessary)
 !DO l=1,2
 !  aux = DOT_PRODUCT(tt,dt(:,l))
 !  dt(:,l) = dt(:,l) - aux*tt
 !END DO

  !                                     evaluates 1st fundamental form
   IF(quad) THEN
     stran(1:3) = 0D0
     DO i=1,3        !for each triangle in the patch
       IF(sides(i))THEN                !if side exists
         ! tangent vectors at mid-side
         tm(:,:,i) = MATMUL(x(:,1:3),cd(1:3,:,i))   !x(1:2)(i) both derivatives
         tm(:,1,i)   = tm(:,1,i) + x(:,i+3)*cd(4,1,i)     !x(1)(i)  add fourth node
         tm(:,2,i)   = tm(:,2,i) + x(:,i+3)*cd(4,2,i)     !x(2)(i)
         stran(1) = stran(1) + DOT_PRODUCT(tm(:,1,i),tm(:,1,i))  !g11
         stran(2) = stran(2) + DOT_PRODUCT(tm(:,2,i),tm(:,2,i))  !g22
         stran(3) = stran(3) + DOT_PRODUCT(tm(:,1,i),tm(:,2,i))  !g12
       ELSE
       !  t(:,:,i) = dx
         stran(1) = stran(1) + DOT_PRODUCT(dx(:,1),dx(:,1))   !g11
         stran(2) = stran(2) + DOT_PRODUCT(dx(:,2),dx(:,2))   !g22
         stran(3) = stran(3) + DOT_PRODUCT(dx(:,1),dx(:,2))   !g12
       END IF
     END DO
     stran(1:3) = stran(1:3)/3d0
   ELSE
     stran(1) = DOT_PRODUCT(dx(:,1),dx(:,1))
     stran(2) = DOT_PRODUCT(dx(:,2),dx(:,2))
     stran(3) = DOT_PRODUCT(dx(:,1),dx(:,2))
   END IF
   !                                      evaluates 2nd fundamental form
   stran(4) = DOT_PRODUCT(dt(:,1),dx(:,1))
   stran(5) = DOT_PRODUCT(dt(:,2),dx(:,2))
   stran(6) = DOT_PRODUCT(dt(:,1),dx(:,2)) + DOT_PRODUCT(dt(:,2),dx(:,1))

   gamma(1) = DOT_PRODUCT(x(:,3)-x(:,2),t(:,1))
   gamma(2) = DOT_PRODUCT(x(:,1)-x(:,3),t(:,2))
   gamma(3) = DOT_PRODUCT(x(:,2)-x(:,1),t(:,3))

   IF ( ASSOCIATED(jas) )THEN
     psib(1) = DOT_PRODUCT(psi(:,3)+psi(:,2),jas(:,1))/2d0
     psib(2) = DOT_PRODUCT(psi(:,1)+psi(:,3),jas(:,2))/2d0
     psib(3) = DOT_PRODUCT(psi(:,2)+psi(:,1),jas(:,3))/2d0
   END IF

 RETURN
 END SUBROUTINE stra07

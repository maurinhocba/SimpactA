 SUBROUTINE stran7(stype,cartd,dn,x,t,stran,dx,dt,gamma,tt,shape,ngaus,nnode,psi,psib,jas,ngamm)
 !***********************************************************************
 !
 !*****this routine computes 1st and 2nd fundamental forms and distorsions
 !            for element 7
 !***********************************************************************
 IMPLICIT NONE
 !              routine parameters
 INTEGER (kind=4), INTENT(IN) :: stype,ngaus,nnode,ngamm

 REAL (kind=8), INTENT(IN) :: cartd(:,:,:),x(:,:),t(:,:),dn(:,:,:),shape(:,:)
 REAL (kind=8), INTENT(IN), OPTIONAL :: psi(:,:)
 REAL (kind=8), INTENT(OUT) :: stran(:,:),dx(:,:,:),dt(:,:,:),gamma(:),tt(:,:)
 REAL (kind=8), INTENT(OUT), OPTIONAL :: psib(:)
 REAL (kind=8), INTENT(IN), POINTER, OPTIONAL :: jas(:,:)

 !              local variables

 INTEGER (kind=4) i,g

 REAL    (kind=8) dfi(3,2,nnode)
 INTERFACE
   INCLUDE 'deriv6.h'
 END INTERFACE

 !     director field, mid-surface & director derivatives at gauss points
 CALL deriv6(ngaus,shape,cartd,x,t,tt,dx,dt)

 DO g = 1,ngaus
   IF (stype < 2) THEN
     !                                     evaluates 1st fundamental form
     stran(1,g) = DOT_PRODUCT(dx(:,1,g),dx(:,1,g))
     stran(2,g) = DOT_PRODUCT(dx(:,2,g),dx(:,2,g))
     stran(3,g) = DOT_PRODUCT(dx(:,1,g),dx(:,2,g))
   END IF
   !                                      evaluates 2nd fundamental form
   stran(4,g) = DOT_PRODUCT(dt(:,1,g),dx(:,1,g))
   stran(5,g) = DOT_PRODUCT(dt(:,2,g),dx(:,2,g))
   stran(6,g) = DOT_PRODUCT(dt(:,1,g),dx(:,2,g))       &
               +DOT_PRODUCT(dt(:,2,g),dx(:,1,g))
 END DO
 !            evaluates 1st fundamental form at assumed strain points
 IF(stype == 2) THEN
   CALL proma1(dfi(1,1,1),x(1,1),dn(1,1,1),3,6,nnode)
   !       dfi = RESHAPE(MULMAT(x,dn),(/ndime,2,3/))
   DO i = 1,3
     stran(1,i) = DOT_PRODUCT(dfi(:,1,i),dfi(:,1,i))
     stran(2,i) = DOT_PRODUCT(dfi(:,2,i),dfi(:,2,i))
     stran(3,i) = DOT_PRODUCT(dfi(:,1,i),dfi(:,2,i))
   END DO
 ELSE IF(stype == 3) THEN
   CALL proma1(dfi(1,1,1),x(1,1),nfdas(1,1,1),3,12,nnode)
   !     dfi = RESHAPE(MULMAT(x,nfdas),(/ndime,2,nnass/))
   stran(1,1) = DOT_PRODUCT(dfi(:,1,1),dfi(:,1,1))
   stran(2,1) = DOT_PRODUCT(dfi(:,1,2),dfi(:,1,2))

   stran(1,3) = DOT_PRODUCT(dfi(:,1,3),dfi(:,1,3))
   stran(2,3) = DOT_PRODUCT(dfi(:,2,3),dfi(:,2,3))
   stran(3,3) = DOT_PRODUCT(dfi(:,1,3),dfi(:,2,3))
   stran(3,1) = (stran(1,3)+stran(2,3)-2d0*stran(3,3))/2d0

   stran(1,3) = DOT_PRODUCT(dfi(:,1,4),dfi(:,1,4))
   stran(2,3) = DOT_PRODUCT(dfi(:,2,4),dfi(:,2,4))
   stran(3,3) = DOT_PRODUCT(dfi(:,1,4),dfi(:,2,4))
   stran(1,2) = (stran(1,3)+stran(2,3)-2d0*stran(3,3))/2d0

   stran(2,2) = DOT_PRODUCT(dfi(:,2,5),dfi(:,2,5))
   stran(3,2) = DOT_PRODUCT(dfi(:,2,6),dfi(:,2,6))

   dfi(:,1:2,1) = MATMUL(x,dn(:,:,1))
   stran(1,3) = DOT_PRODUCT(dfi(:,1,1),dfi(:,1,1))
   stran(2,3) = DOT_PRODUCT(dfi(:,2,1),dfi(:,2,1))
   stran(3,3) = 2d0*DOT_PRODUCT(dfi(:,1,1),dfi(:,2,1))
 END IF
 CALL asstr7(x,t,gamma,psi,psib,jas,ngamm)

 RETURN
 END SUBROUTINE stran7

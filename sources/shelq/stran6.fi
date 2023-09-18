 SUBROUTINE stran6(ngaus,cartd,x,t,stran,dx,dt)
 !***********************************************************************
 !
 !*****this routine computes 1st and 2nd fundamental forms
 !     for shell element
 !***********************************************************************
 IMPLICIT NONE

 !              routine parameters

 INTEGER (kind=4),INTENT(IN) :: ngaus

 REAL (kind=8),INTENT(IN) :: cartd(:,:,:),x(:,:),t(:,:)
 REAL (kind=8),INTENT(OUT) :: stran(:,:),dx(:,:,:),dt(:,:,:)
 !              local variables

 INTEGER (kind=4) :: g

 !     mid-surface & director derivatives at gauss points

 DO g=1,ngaus
   dt(:,:,g) = MATMUL(t,cartd(:,:,g))
   dx(:,:,g) = MATMUL(x,cartd(:,:,g))
   !                                     evaluates 1st fundamental form
   stran(1,g)= DOT_PRODUCT(dx(:,1,g),dx(:,1,g))
   stran(2,g)= DOT_PRODUCT(dx(:,2,g),dx(:,2,g))
   stran(3,g)= DOT_PRODUCT(dx(:,1,g),dx(:,2,g))
   !                                      evaluates 2nd fundamental form
   stran(4,g)= DOT_PRODUCT(dt(:,1,g),dx(:,1,g))
   stran(5,g)= DOT_PRODUCT(dt(:,2,g),dx(:,2,g))
   stran(6,g)= DOT_PRODUCT(dt(:,1,g),dx(:,2,g)) & !twice the torsion
              +DOT_PRODUCT(dt(:,2,g),dx(:,1,g))
 END DO
 RETURN
 END SUBROUTINE stran6

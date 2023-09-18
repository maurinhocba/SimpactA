 SUBROUTINE deriv6(ngaus,shape,cartd,x,t,tt,dx,dt)
 !***********************************************************************
 !
 !*****this routine computes configuration derivatives at gauss points
 !           for shell elements
 !***********************************************************************
 IMPLICIT NONE
 !                   routine parameters

 INTEGER (kind=4),INTENT(IN) :: ngaus

 REAL (kind=8),INTENT(IN) :: shape(:,:),cartd(:,:,:),x(:,:),t(:,:)
 REAL (kind=8),INTENT(OUT) :: tt(:,:),dx(:,:,:),dt(:,:,:)
 END SUBROUTINE deriv6

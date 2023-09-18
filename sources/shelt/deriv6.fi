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

 !                   local variables
 INTEGER (kind=4) :: g,l
 REAL (kind=8) :: norma,aux

 INTERFACE
   INCLUDE 'vecuni.h'
 END INTERFACE

 !     evaluates director & configuration & director derivatives
 tt = MATMUL(t,shape)

 DO g = 1,ngaus
   dx(:,:,g) = MATMUL(x,cartd(:,:,g))
   dt(:,:,g) = MATMUL(t,cartd(:,:,g))
   !  normalizes director field
   CALL vecuni(3,tt(:,g),norma)

   !   correct director derivat.
   DO l=1,2
     aux = DOT_PRODUCT(tt(:,g),dt(:,l,g))
     dt(:,l,g) = (dt(:,l,g)-aux*tt(1:,g))/norma
   END DO
 END DO

 RETURN
 END SUBROUTINE deriv6

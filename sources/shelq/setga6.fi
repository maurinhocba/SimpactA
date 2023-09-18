 SUBROUTINE setga6(nnode,ngaus,ngamm,    x,    t,                   &
&                  shape,deriv,weigp,dvolu,                         &
&                  cartd,jacin,stra0,gamm0,ambda,angle,locax)
 !***********************************************************************
 !
 !****this routine sets up some needed constant DATA for shell element
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters

 INTEGER (kind=4), INTENT(IN) :: nnode,ngaus,ngamm,locax

 REAL (kind=8), INTENT(IN) :: x(:,:),t(:,:),              &
&                 shape(:,:),deriv(:,:,:),angle,          &
&                 weigp(:)
 REAL (kind=8), INTENT(OUT) :: dvolu(:),cartd(:,:,:),         &
&                 jacin(:,:,:),stra0(:,:),gamm0(:),           &
&                 ambda(:,:)

 !                   local variables

 INTEGER (kind=4) :: g,n,l
 REAL (kind=8) :: tgp(3),norml(3),modul,dx(3,2,ngaus), &
                  dt(3,2,ngaus)

 INTERFACE
   INCLUDE 'lcasy6.h'
   INCLUDE 'vecuni.h'
 END INTERFACE

 DO g=1,ngaus
   !       compute the local cartesian system
   tgp = MATMUL(t,shape(:,g))
   CALL vecuni(3,tgp,modul)
   CALL lcasy6(deriv(:,:,g),x,jacin(:,:,g),dvolu(g),norml,angle,locax)
   ambda(1,g) = DOT_PRODUCT(norml,tgp)
   ambda(2,g) = 1d0
   dvolu(g) = dvolu(g)*weigp(g)
   DO n=1,nnode
     DO l=1,2
       cartd(n,l,g)=jacin(l,1,g)*deriv(n,1,g)+jacin(l,2,g)*deriv(n,2,g)
     END DO
   END DO
 END DO
 !     compute initial 1st and 2nd fundamental forms (& angles with t)
 CALL stran6(ngaus,cartd,x,t,stra0,dx,dt)
 !                  evaluates distorsion at the assumed strain points
 !     here is only valid for 4-noded element because of the assumed points
 !     we could use some dummy variables to use it for triangles also
 CALL asstr6(ngamm,x,t,gamm0)

 RETURN
 END SUBROUTINE setga6

 SUBROUTINE setga7(stype,dn,x,t,deriv,dvolu,cartd,jacin,stra0,gamm0,ambda, &
                   angle,locax,shape,weigp,ngaus,nnode,ngamm)
 !***********************************************************************
 !
 !****this routine sets up some needed constant DATA for shell element
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters

 INTEGER (kind=4), INTENT(IN) :: stype,locax,ngaus,nnode,ngamm
 REAL (kind=8), INTENT(IN) :: x(:,:),t(:,:),dn(:,:,:),deriv(:,:,:), angle, &
                              shape(:,:),weigp(:)
 REAL (kind=8), INTENT(OUT) :: dvolu(:),cartd(:,:,:),jacin(:,:,:), &
                               stra0(:,:),gamm0(:),ambda(:,:)

 !                   local variables
 INTEGER (kind=4) g,n,l
 REAL (kind=8) tgp(3),norml(3),dx(3,2,ngaus),dt(3,2,ngaus),modul,tt(3,ngaus)
 REAL (kind=8),POINTER :: nada(:,:)
 

 INTERFACE
   INCLUDE 'lcasy6.h'
   INCLUDE 'vecuni.h'
 END INTERFACE

 NULLIFY(nada)
 DO g=1,ngaus
   tgp = MATMUL(t,shape(:,g))
   CALL vecuni(3,tgp,modul)
   !       compute the local cartesian system
   CALL lcasy6(deriv(:,:,g),x,jacin(:,:,g),dvolu(g),norml,angle,locax)
   dvolu(g) = dvolu(g)*weigp(g)
   ambda(1,g) = DOT_PRODUCT(norml,tgp)
   ambda(2,g) = 1d0
   !cartd(:,:,g) = MATMUL( deriv(:,:,g), TRANSPOSE(jacin(:,:,g)) )
   DO n = 1,nnode
     DO l = 1,2
       cartd(n,l,g) = jacin(l,1,g)*deriv(n,1,g)+jacin(l,2,g)*deriv(n,2,g)
     END DO
   END DO
 END DO
 !     compute initial 1st and 2nd fundamental forms (& angles with t)
 CALL stran7(stype,cartd,dn,x,t,stra0,dx,dt,gamm0,tt,shape,ngaus,nnode,tt,tgp,nada,ngamm)
 RETURN
 END SUBROUTINE setga7

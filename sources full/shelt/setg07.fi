 SUBROUTINE setg07(x,t,dvolu,cartd,jacin,stra0,gamm0,ambda,angle,locax,quad,mems,cd)
 !***********************************************************************
 !
 !****this routine sets up some needed constant DATA for shell element
 !
 !***********************************************************************
 IMPLICIT NONE

 !                   routine parameters

 INTEGER (kind=4), INTENT(IN) :: locax
 REAL (kind=8), INTENT(IN) :: x(:,:),t(:,:),angle
 REAL (kind=8), INTENT(OUT) :: dvolu,cartd(:,:),jacin(:,:), &
                               stra0(:),gamm0(:),ambda(:)
 LOGICAL, INTENT(IN) :: quad,mems(:)
 REAL(kind=8), POINTER :: cd(:,:,:)
 REAL (kind=8),POINTER :: nada(:,:)

 !                   local variables

 REAL (kind=8) tgp(3),norml(3),dx(3,2),dt(3,2),modul,tm(3,2,3)

 INTERFACE
   INCLUDE 'vecuni.h'
 END INTERFACE

   NULLIFY(nada)
   tgp = ( t(:,1) + t(:,2) + t(:,3) )/3d0  !director at element center
   CALL vecuni(3,tgp,modul)                !unit vector
   !       compute the local cartesian system
   CALL axep07(x,cartd,cd,dvolu,mems,angle,locax,quad,norml)
   dvolu = dvolu/2d0
   jacin = RESHAPE( (/ cartd(2,1), cartd(2,2), cartd(3,1), cartd(3,2) /),(/2,2/))
   ambda(1) = DOT_PRODUCT(norml,tgp)
   ambda(2) = 1d0
   !     compute initial 1st and 2nd fundamental forms (& angles with t)
   CALL stra07(cartd,x,t,stra0,dx,dt,gamm0,quad,mems,tm,cd,t,tgp,nada)

 RETURN
 END SUBROUTINE setg07

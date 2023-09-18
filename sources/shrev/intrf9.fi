 SUBROUTINE intrf9(nstre,ntype,plast,nnode,ngaus,matty,axesc,     &
                   weigh,shape,deriv,dmatx,props,prop1,           &
                   jac,r0,stra0,stres,sede0,sedef,                &
                   ehist,strap,bts,x,phig,phil,                   &
                   energ,istop,tn,thick)
 !***********************************************************************
 !
 !****this routine compute the internal nodal forces for element no. 9
 !    2-3 node shell/beam  element
 !***********************************************************************
 IMPLICIT NONE

 !***  routine parameters
 LOGICAL, INTENT(IN) :: plast
 INTEGER(kind=4), INTENT(IN) :: ngaus,nstre,nnode,matty,axesc,ntype
 INTEGER(kind=4), INTENT(OUT) :: istop
 REAL(kind=8), INTENT(IN) :: dmatx(:),props(:),prop1(:),             &
                             shape(:,:),deriv(:,:),weigh(:),         &
                             jac(:),r0(:),stra0(:,:),sede0(:),       &
                             x(:,:),phig(:),phil(:),thick
 REAL(kind=8), INTENT(OUT) :: stres(:,:),sedef(:),bts(:,:),energ(:),tn(:,:)
 REAL (kind=8), POINTER :: ehist(:,:),strap(:,:)
 !***  local variables
 INTEGER (kind=4) n,g
 REAL    (kind=8) phin(nnode),facto,tg(2,ngaus),dx(2,ngaus),    &
                  dstra(nstre),stran(nstre,ngaus)
 REAL (kind=8),PARAMETER :: twopi=6.283185307179586d0

 !     compute the local cartesian system
 IF(axesc == 0) THEN
   phin = phig
 ELSE
   DO n=1,nnode
     phin(n) = phig(n) + phil(MIN(axesc,n))
   END DO
 END IF
 IF(ABS(phin(nnode)-phin(1)) > 3.2d0) THEN
   DO n=1,nnode
      IF(phin(n) < 0) phin(n) = phin(n)+twopi
   END DO
 END IF
 !     compute actual angles and curvatures
 CALL stran9(nnode,ngaus,nstre,ntype,x,dx,phin,stran,shape,deriv,jac,r0,tg)
 !     evaluates incremental and total stresses at gauss points
 CALL istgp9(matty,ngaus,nstre,plast,stres,stra0,sede0,sedef,ehist, &
             strap,dmatx,props,prop1,stran,energ,jac,weigh,istop)
 !     sums contribution to residual forces
 bts = 0d0
 DO g=1,ngaus
   facto = jac(g)*weigh(g)
   IF(ntype == 3) facto = facto*r0(g)*twopi
   dstra = facto*stres(:,g)
   CALL btstr9(nnode,ntype,deriv(:,g),shape(:,g),phin,x, &
               r0(g),bts,dx(:,g),jac(g),dstra)
 END DO
 facto = thick*SUM(sedef)/ngaus
 DO n=1,nnode
   tn(1,n) = -SIN(phin(n))*facto
   tn(2,n) =  COS(phin(n))*facto
 END DO
 RETURN
 END SUBROUTINE intrf9

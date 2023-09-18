 SUBROUTINE intrf8(ndime,ndofe,ng,nstre,nn,il,prop1,                &
                   weigh,shape,deriv,jac,stran,stres,epdef,sedef,   &
                   bmsig,x,v,lbd,llbd,istop,ielem,isotro)
 !***********************************************************************
 !
 !****this routine compute the internal nodal forces for element no. 8
 !    2-3 node beam element
 !***********************************************************************
 IMPLICIT NONE
 !***  routine parameters
 INTEGER (kind=4) ndime,ndofe,ng,nstre,nn,il,istop,ielem
 REAL (kind=8) prop1(:),shape(nn,ng),deriv(nn,ng),weigh(ng),        &
               jac(ng),stran(nstre,ng),stres(nstre,ng),             &
               epdef(ng),sedef(ng),bmsig(ndofe,nn),x(3,nn),         &
               v(3,nn),lbd(3,3,nn),llbd(3,3,nn)
 LOGICAL :: isotro
 !***  local variables
 INTEGER (kind=4) n,g,k
 REAL (kind=8) iro(3),vg(3),t(6),st(6),fu,dn,dstra(nstre,ng),       &
               bts(ndofe,nn),lb(3,3,nn),dx(3,ng),vb(3,nn)
 REAL (kind=8),SAVE :: h(3) = (/0d0,0d0,0d0/)

 !     compute the local cartesian system
 IF(il == 0) THEN  !same systems
   lb = lbd             !local system = global system
   vb = v               !local velocities = global velocities
 ELSE              !relative systems
   DO n=1,nn       !for each node
     k = MIN(il,n)
     lb(1:3,1:3,n) = MATMUL(lbd(1:3,1:3,n),llbd(1:3,1:3,k))
     vb(1:3,n) = MATMUL(TRANSPOSE(llbd(1:3,1:3,k)),v(1:3,n))
   END DO
 END IF
 !     SAVE former angles and curvatures
 dstra = stran  !original strains
 !     compute actual angles and curvatures
 CALL stran8(nn,ng,ndime,nstre,x,dx,lb,stran,shape,deriv,jac,lbd)
 !     evaluates incremental strain measures
 dstra = stran - dstra
 !     evaluates incremental and total stresses at gauss points
 IF( isotro )THEN
   CALL istg18(ng,nstre,stres,dstra,epdef,sedef,prop1,istop,ielem)
   IF(istop == 1) RETURN
   iro(1) = (prop1(6)-prop1(5))*prop1(15)/prop1(16)
   iro(2) =  prop1(5)*prop1(15)/prop1(16)
   iro(3) = -prop1(6)*prop1(15)/prop1(16)
 ELSE
   CALL istg28(ng,nstre,stres,dstra,sedef,prop1)
   iro(1) =  prop1(4)-prop1(3)
   iro(2) =  prop1(3)
   iro(3) = -prop1(4)
 END IF
 bts = 0d0
 !     sums contribution to residual forces
 DO g=1,ng
   vg = MATMUL(vb,shape(1:nn,g))
   h(1) = iro(1)*vg(2)*vg(3)
   h(2) = iro(2)*vg(3)*vg(1)
   h(3) = iro(3)*vg(1)*vg(2)
   t  = stran(1:6,g)
   st = jac(g)*weigh(g)*stres(1:6,g)
   DO n=1,nn
     fu = shape(n,g)
     dn = deriv(n,g)/jac(g)
     bts(1:3,n) = bts(1:3,n) + dn*MATMUL(lbd(1:3,1:3,g),st(1:3))
     bts(4,n)= bts(4,n) + dn*st(4) + fu*(h(1) + t(3)*st(2)          &
                        - t(2)*st(3) + t(6)*st(5) - t(5)*st(6))
     bts(5,n)= bts(5,n) + dn*st(5) + fu*(h(2) + t(1)*st(3)          &
                        - t(3)*st(1) + t(4)*st(6) - t(6)*st(4))
     bts(6,n)= bts(6,n) + dn*st(6) + fu*(h(3) + t(2)*st(1)          &
                        -t(1)*st(2) + t(5)*st(4) - t(4)*st(5))
   END DO
 END DO
 !     transform to nodal system from beam local system
 DO n=1,nn
   bmsig(1:3,n) = bts(1:3,n)
   IF(il == 0) THEN
     bmsig(4:6,n) = bts(4:6,n)
   ELSE
     k = MIN(n,il)
     bmsig(4:6,n) = MATMUL(llbd(1:3,1:3,k),bts(4:6,n))
   END IF
 END DO
 RETURN
 END SUBROUTINE intrf8

 SUBROUTINE intrf8(ndime,ndofe,ng,nstre,nn,il,rprop,weigh,shape,  &
                   deriv,jac,stra0,stran,stres,epdef,sedef,       &
                   bmsig,x,v,lbd,llbd,istop,ielem,isotro,visco,plast)
 !***********************************************************************
 !
 !****this routine compute the internal nodal forces for element no. 8
 !    2-3 node beam element
 !***********************************************************************
 IMPLICIT NONE
 !***  routine parameters
 INTEGER (kind=4) ndime,ndofe,ng,nstre,nn,il,istop,ielem
 REAL (kind=8) rprop(:),shape(nn,ng),deriv(nn,ng),weigh(ng),             &
               jac(ng),stra0(nstre,ng),stran(nstre,ng),stres(nstre,ng),  &
               epdef(ng),sedef(ng),bmsig(ndofe,nn),x(3,nn),              &
               v(6,nn),lbd(3,3,nn),llbd(3,3,nn)
 LOGICAL :: isotro,visco,plast
 !***  local variables
 INTEGER (kind=4) n,g,k,nv
 REAL (kind=8) iro(3),vg(6),t(6),st(6),fu,dn,dstra(nstre,ng),       &
               bts(ndofe,nn),lb(3,3,nn),dx(3,ng),vb(6,nn),          &
               bm(6,nn,nstre),vst(nstre)

 REAL (kind=8),SAVE :: h(3) = (/0d0,0d0,0d0/)

 !     compute the local cartesian system
 IF(il == 0) THEN  !same systems
   lb = lbd             !local system = global system
   vb = v               !local velocities = global velocities
 ELSE              !relative systems
   DO n=1,nn       !for each node
     k = MIN(il,n)
     lb(1:3,1:3,n) = MATMUL(lbd(1:3,1:3,n),llbd(1:3,1:3,k))     !element nodal local systems
     vb(1:3,n) = v(1:3,n)                                       !translational velocities
     vb(4:6,n) = MATMUL(TRANSPOSE(llbd(1:3,1:3,k)),v(4:6,n))    !rotation in element local system
   END DO
 END IF
 !     compute actual angles and curvatures
 IF( plast )dstra = stran
 CALL stran8(nn,ng,ndime,nstre,x,dx,lb,stran,shape,deriv,jac,lbd)
 !     evaluates incremental/total strain measures
 IF( plast )THEN !incremental strains
   dstra = stran - dstra
 ELSE            !total strains
   dstra = stran - stra0
  END IF
 ! rotational inertia
 iro(1) =  rprop(4)-rprop(3)   !I3 - I2
 iro(2) =  rprop(2)-rprop(4)   !I1 - I3 =  I2
 iro(3) =  rprop(3)-rprop(2)   !I2 - I1 = -I3
 !     evaluates incremental and total stresses at gauss points
 IF( isotro )THEN
   CALL istg18(ng,nstre,stres,dstra,epdef,sedef,rprop,istop,ielem,plast)
   IF(istop == 1) RETURN
   IF( visco )nv = 10
 ELSE
   CALL istg28(ng,nstre,stres,dstra,sedef,rprop)
   IF( visco )nv = 25
 END IF
 bts = 0d0
 !     sums contribution to residual forces
 DO g=1,ng
   vg = MATMUL(vb,shape(1:nn,g))
   h(1) = iro(1)*vg(5)*vg(6)       ! w x h = w x (Iw)
   h(2) = iro(2)*vg(6)*vg(4)
   h(3) = iro(3)*vg(4)*vg(5)
   t  = stran(1:6,g)
   IF( visco )THEN
     CALL bmatx8(nstre,nn,shape(:,g),deriv(:,g),jac(g),t,lbd(:,:,g),bm)
     CALL proma1(vst,vb,bm,1,nstre,ndofe*nn)
     !IF( isotro ) THEN
       DO k=1,nstre
         stres(k,g) = stres(k,g) + rprop(nv+k)*vst(k)
       END DO
     !ELSE   ! full version
     !  vst = vst*rprop(32)
     !  stres(1,g) = stres(1,g) +rprop( 5)*vst(1) +rprop(11)*vst(2) +rprop(12)*vst(3) &
     !                          +rprop(13)*vst(4) +rprop(14)*vst(5) +rprop(15)*vst(6)
     !  stres(2,g) = stres(2,g) +rprop(11)*vst(1) +rprop( 6)*vst(2) +rprop(16)*vst(3) &
     !                          +rprop(17)*vst(4) +rprop(18)*vst(5) +rprop(19)*vst(6)
     !  stres(3,g) = stres(3,g) +rprop(12)*vst(1) +rprop(16)*vst(2) +rprop( 7)*vst(3) &
     !                          +rprop(20)*vst(4) +rprop(21)*vst(5) +rprop(22)*vst(6)
     !  stres(4,g) = stres(4,g) +rprop(13)*vst(1) +rprop(17)*vst(2) +rprop(20)*vst(3) &
     !                          +rprop( 8)*vst(4) +rprop(23)*vst(5) +rprop(24)*vst(6)
     !  stres(5,g) = stres(5,g) +rprop(14)*vst(1) +rprop(18)*vst(2) +rprop(21)*vst(3) &
     !                          +rprop(23)*vst(4) +rprop( 9)*vst(5) +rprop(25)*vst(6)
     !  stres(6,g) = stres(6,g) +rprop(15)*vst(1) +rprop(19)*vst(2) +rprop(22)*vst(3) &
     !                          +rprop(24)*vst(4) +rprop(25)*vst(5) +rprop(10)*vst(6)
     !END IF
   END IF
   st = jac(g)*weigh(g)*stres(1:6,g)
   h  = jac(g)*weigh(g)*h

   ! CALL proma3(bts(1,1),bm(1,1,1),st(1),ndofe*nn,nstre,.FALSE.)  !nodal contributions
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

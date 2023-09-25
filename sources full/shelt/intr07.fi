 SUBROUTINE intr07(ndofe,nstre,ngaus,ngamm,ap1t,plast,matty,nlayr,dmatx,props,sf,x,  &
                   t,lb,tgp,dvolu,cartd,jacin,cd,mems,stra0,gamm0,strsg,q,psi,nangl, &
                   jas,ambda,ehist,strap,stres,bmsig,                                &
                   energ,thick,istop,ielem,stabq,quad,ds)
 !***********************************************************************
 !
 !****this routine compute the internal nodal forces for element no. 3
 !    6-node linear shell element TLLL
 !***********************************************************************
 IMPLICIT NONE
 !***routine parameters
 INTEGER (kind=4), INTENT(IN) :: ndofe,nstre,ngaus,ngamm,matty,nlayr,ielem
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: props(:),dmatx(:),x(:,:),t(:,:), &
                lb(:,:,:),dvolu,cartd(:,:),ap1t(:,:),psi(:,:), &
                jacin(:,:),stra0(:),gamm0(:),stabq,ds

 REAL (kind=8), INTENT(IN) :: sf(3)
 REAL (kind=8), INTENT(IN OUT) :: thick,strsg(:),ambda(:),energ(:),    &
                bmsig(:,:),q(:),tgp(:)
 REAL(kind=8), POINTER :: ehist(:,:),strap(:,:),stres(:,:,:),cd(:,:,:),nangl(:,:),jas(:,:)
 LOGICAL, INTENT(IN) :: plast,quad,mems(:)
 !***local variables
 INTEGER (kind=4) i

 REAL (kind=8) s(2),def,norml(3),bsmtx(5,3,ngamm),bbe(5,3,3),str(nstre,ngaus),sta(nstre,ngaus), &
               auxil(nstre),gamma(ngamm),dx(3,2),dt(3,2),tm(3,2,3),daux(6),f2,     &
               aa,bb,cc,dm(3),a6,bm(3,6,3),bfs(2,3,ngamm),bfb(2,3,4),psib(ngamm),qf(ngamm)

 INTERFACE
   INCLUDE 'vecuni.h'
   INCLUDE 'stgp06.h'
 END INTERFACE

 ! director at center
 tgp = (t(:,1)+t(:,2)+t(:,3))/3d0
 CALL vecuni(3,tgp,def)
 !     evaluates first and second fundamental forms and distorsions
 CALL stra07(cartd,x,t,auxil,dx,dt,gamma,quad,mems,tm,cd,psi,psib,jas) !=>auxil & gamma

 !     computes shear strains substracting the initial distorsion -->dstra
 gamma = gamma - gamm0     !natural shear strains at selected points
 CALL vecpro(dx(:,1),dx(:,2),norml)
 ambda(2)= ambda(1)/DOT_PRODUCT(norml,tgp)
 !ambda(2)= ambda(1)/(auxil(1)*auxil(2)*4d0-auxil(3)**2)
 ! computes strains substracting the initial configuration -->auxil
 IF(nlayr == 1 )THEN
   auxil(1:6) = auxil(1:6) - stra0(:)
   auxil(1:2) = auxil(1:2)/2d0
 ELSE ! modify curvatures only for large-strain multilayer
   auxil(4:5) = (auxil(4:5) - stra0(4:5))*2d0  !twice the curvatures
   auxil(6) = auxil(6) - stra0(6)
 END IF
 ! compute shear strains
 s = MATMUL(ap1t,gamma)               !natural shear strains at element center
 auxil(7:8) = MATMUL(jacin,s)    !cartesyan shear strains at element center
 IF( ndofe == 7 ) THEN    ! psi gradient
   auxil( 9) = DOT_PRODUCT(psi(1,:),cartd(:,1))
   auxil(10) = DOT_PRODUCT(psi(2,:),cartd(:,2))
   auxil(11) = DOT_PRODUCT(psi(1,:),cartd(:,2))
   auxil(12) = DOT_PRODUCT(psi(2,:),cartd(:,1))
   s = MATMUL(ap1t,psib)
   auxil(13:14) = MATMUL(jacin,s)
 END IF

 ! evaluates total stresses at gauss points --> strsg
 !WRITE(55,"('def',3e15.6,/,5E15.6)")auxil
 IF( plast )THEN
   strap(:,1) =  strap(:,1)     !initializes current values
   ehist(1,1) = ehist(3,1)          !equivalent plastic strain
   ehist(2,1) = 0d0               !equivalent stress
   ehist(4:5,1) = 0d0             !consistent parameters
   IF(nlayr == 1)THEN
     ! computes elastic strains --> auxil
     auxil = auxil - strap(:,1)
     CALL radre6(matty,dx(1,1),dx(1,2),ambda(1),           &
                 tgp(1),auxil(1),strsg(1),strap(1,1),      &
                 ehist(1,1),dmatx(1),props(1),istop)
   ELSE
     CALL slayer(1,nstre,strsg(1),auxil(1),dmatx(1),            &
                 props(1),strap(1,1),stres(1,1,1),ambda(1),nlayr, &
                 thick,istop,ehist(1,1),ielem,plast)
   END IF

 ELSE     !elastic
   IF(nlayr == 1)THEN
     IF( ndofe == 5 )THEN
       CALL istgp6(1,nstre,strsg(1),auxil(1),dmatx(1), &
                   ambda(1),matty,props(1),sf(1))
     ELSE      ! zigzag
       sta(:,1) = auxil
       CALL stgp06(ngaus,str,sta,dmatx)
       strsg = str(:,1)
     END IF
   ELSE
     CALL slayer(1,nstre,strsg(1),auxil(1),dmatx(1),       &
                 props(1),s(1),s(2),ambda(1),nlayr,        &
                 thick,istop,s(2),ielem,plast)
   END IF
   !DO n=1,nstre
   !  energ(n) = energ(n) + dvolu*strsg(n)*auxil(n)
   !END DO
 END IF

 !       scales stresses for numerical integration
 auxil(:) = dvolu*strsg(:)
 !       computes b matrix , bending & membrane
 CALL bmat07(cartd,lb,bbe,bm,dx,dt,quad,mems,tm,cd,ndofe,nangl,bfb)
 !                computes bsshear matrix
 CALL bshe07(x,t,lb,bsmtx,bfs,jas,nangl)
 IF( ndofe == 5 )THEN
   !       equivalent q for transverse shear (one-point-quadrature)
   s = MATMUL(TRANSPOSE(jacin),auxil(7:8))
   q = MATMUL(TRANSPOSE(ap1t),s)           !A(C)*(J^(-T)*Q*area)
   ! + stabilization
   dm(1:2) = ds; dm(3) = 0d0     !matty = 1-2
   aa = dm(1)*jacin(1,1)**2+dm(2)*jacin(2,1)**2+dm(3)*jacin(1,1)*jacin(2,1)*2d0
   bb = dm(1)*jacin(1,2)**2+dm(2)*jacin(2,2)**2+dm(3)*jacin(2,2)*jacin(1,2)*2d0
   cc = dm(1)*jacin(1,1)*jacin(1,2)+dm(2)*jacin(2,2)*jacin(2,1)+dm(3)*jacin(1,2)*jacin(2,1)*2d0
   a6 = dvolu*stabq/18d0*(+aa+bb+cc)*SUM(gamma) !(gamma(1)+gamma(2)+gamma(3))
   q = q + a6
 ELSE
   a6 = dvolu/6d0
   aa = dmatx(13)*jacin(1,1)**2+dmatx(15)*jacin(2,1)**2+dmatx(14)*jacin(1,1)*jacin(2,1)*2d0
   bb = dmatx(13)*jacin(1,2)**2+dmatx(15)*jacin(2,2)**2+dmatx(14)*jacin(2,2)*jacin(1,2)*2d0
   cc = dmatx(13)*jacin(1,1)*jacin(1,2)+dmatx(15)*jacin(2,2)*jacin(2,1)+dmatx(14)*jacin(1,2)*jacin(2,1)*2d0

   daux( 1) = (+    aa+    bb-    cc)*a6  !ds(1,1)
   daux( 2) = (+    aa-    bb+    cc)*a6  !ds(2,1)
   daux( 3) = (-    aa+    bb+    cc)*a6  !ds(3,1)
   daux( 4) = (+    aa+3d0*bb+3d0*cc)*a6  !ds(2,2)
   daux( 5) = (-    aa-    bb-3d0*cc)*a6  !ds(3,2)
   daux( 6) = (+3d0*aa+    bb+3d0*cc)*a6  !ds(3,3)
   q(1) = daux(1)*gamma(1)+daux(2)*gamma(2)+daux(3)*gamma(3)
   q(2) = daux(2)*gamma(1)+daux(4)*gamma(2)+daux(5)*gamma(3)
   q(3) = daux(3)*gamma(1)+daux(5)*gamma(2)+daux(6)*gamma(3)

   f2 = dmatx(63)+dmatx(64)         !this is a simplification? check
   aa = dmatx(62)*jacin(1,1)**2+dmatx(65)*jacin(2,1)**2+f2*jacin(1,1)*jacin(2,1)
   bb = dmatx(62)*jacin(1,2)**2+dmatx(65)*jacin(2,2)**2+f2*jacin(2,2)*jacin(1,2)
   cc = dmatx(62)*jacin(1,1)*jacin(1,2)+dmatx(65)*jacin(2,2)*jacin(2,1)+f2*jacin(1,2)*jacin(2,1)

   daux( 1) = (+    aa+    bb-    cc)*a6  !ds(1,1)
   daux( 2) = (+    aa-    bb+    cc)*a6  !ds(2,1)
   daux( 3) = (-    aa+    bb+    cc)*a6  !ds(3,1)
   daux( 4) = (+    aa+3d0*bb+3d0*cc)*a6  !ds(2,2)
   daux( 5) = (-    aa-    bb-3d0*cc)*a6  !ds(3,2)
   daux( 6) = (+3d0*aa+    bb+3d0*cc)*a6  !ds(3,3)

   q(1) = q(1) +daux(1)*psib(1)+daux(2)*psib(2)+daux(3)*psib(3)
   q(2) = q(2) +daux(2)*psib(1)+daux(4)*psib(2)+daux(5)*psib(3)
   q(3) = q(3) +daux(3)*psib(1)+daux(5)*psib(2)+daux(6)*psib(3)

   qf(1) = +daux(1)*gamma(1)+daux(2)*gamma(2)+daux(3)*gamma(3)
   qf(2) = +daux(2)*gamma(1)+daux(4)*gamma(2)+daux(5)*gamma(3)
   qf(3) = +daux(3)*gamma(1)+daux(5)*gamma(2)+daux(6)*gamma(3)

   aa = dmatx(59)*jacin(1,1)**2+dmatx(61)*jacin(2,1)**2+dmatx(60)*jacin(1,1)*jacin(2,1)*2d0
   bb = dmatx(59)*jacin(1,2)**2+dmatx(61)*jacin(2,2)**2+dmatx(60)*jacin(2,2)*jacin(1,2)*2d0
   cc = dmatx(59)*jacin(1,1)*jacin(1,2)+dmatx(61)*jacin(2,2)*jacin(2,1)+dmatx(60)*jacin(1,2)*jacin(2,1)*2d0

   daux( 1) = (+    aa+    bb-    cc)*a6  !ds(1,1)
   daux( 2) = (+    aa-    bb+    cc)*a6  !ds(2,1)
   daux( 3) = (-    aa+    bb+    cc)*a6  !ds(3,1)
   daux( 4) = (+    aa+3d0*bb+3d0*cc)*a6  !ds(2,2)
   daux( 5) = (-    aa-    bb-3d0*cc)*a6  !ds(3,2)
   daux( 6) = (+3d0*aa+    bb+3d0*cc)*a6  !ds(3,3)

   qf(1) = qf(1) + +daux(1)*psib(1)+daux(2)*psib(2)+daux(3)*psib(3)
   qf(2) = qf(2) + +daux(2)*psib(1)+daux(4)*psib(2)+daux(5)*psib(3)
   qf(3) = qf(3) + +daux(3)*psib(1)+daux(5)*psib(2)+daux(6)*psib(3)

 END IF

  DO i=1,3
   bmsig(1:3,i) =   bmsig(1:3,i) + MATMUL(bm(:,i,:),auxil(1:3))   !membrane part
   IF( quad ) bmsig(1:3,i+3) =   bmsig(1:3,i+3) + MATMUL(bm(:,i+3,:),auxil(1:3))
   bmsig(1:5,i) =   bmsig(1:5,i) + MATMUL(bbe(:,i,:),auxil(4:6))  !bending part
   bmsig(1:5,i) =   bmsig(1:5,i) + MATMUL(bsmtx(:,i,:),q)         !shear part
   IF(ndofe > 5 ) bmsig(6:7,i) = bmsig(6:7,i) + MATMUL(bfs(:,i,:),qf) + MATMUL(bfb(:,i,:),auxil(9:12))
 END DO

 thick = thick*ambda(2)  !present thickness
 RETURN
 END SUBROUTINE intr07

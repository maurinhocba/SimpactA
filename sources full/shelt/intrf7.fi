 SUBROUTINE intrf7(nvarl,plast,ndofe,nstre,        &
                   shape,posgp,ngaus,ngamm,nnode,  &
                   stype,matty,nlayr,dmatx,props,sf,x,t,lb,psi,nangl,jas,&
                   dvolu,cartd,jacin,stra0,gamm0,strsg,q,beta, &
                   ambda,ehist,strap,stres,bmsig,              &
                   dn,ap1,energ,thick,istop,ielem)
 !***********************************************************************
 !
 !****this routine compute the internal nodal forces for element no. 7
 !    6-node triangular shell element
 !***********************************************************************
 IMPLICIT NONE
 !***routine parameters
 INTEGER (kind=4), INTENT(IN) :: nvarl,stype,matty,nlayr,ielem,ngaus,ngamm,nnode,ndofe,nstre
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: props(:),dmatx(:),x(:,:),t(:,:),             &
                lb(:,:,:),dvolu(:),cartd(:,:,:),psi(:,:),                  &
                jacin(:,:,:),stra0(:,:),gamm0(:),                          &
                dn(:,:,:),ap1(:,:,:),shape(:,:),posgp(:,:)
 REAL (kind=8), INTENT(IN) :: sf(3)

 REAL (kind=8), INTENT(IN OUT) :: thick,strsg(:,:),ambda(:,:),energ(:),    &
                bmsig(:,:),q(:)
 REAL(kind=8), POINTER :: ehist(:,:),strap(:,:),stres(:,:,:),beta(:),nangl(:,:),jas(:,:)
 LOGICAL, INTENT(IN) :: plast
 !***local variables
 INTEGER (kind=4) n,g,i,j,k
 REAL (kind=8)    bmatx(ndofe,nnode,nstre),             &
                  bsmtx(5,nnode,ngamm),gmatx(3,nnode,nbeta)
 REAL (kind=8) s(3),xi(3),def,norml(3),tgp(3,ngaus),qf(ngamm),            &
               auxil(nstre,ngaus),qaux(ngamm),gamma(ngamm),psib(ngamm),   &
               dx(3,2,ngaus),dt(3,2,ngaus),baux(nbeta),e(3,ngaus),bfmtx(2,nnode,ngamm)

 INTERFACE
   INCLUDE 'vecuni.h'
   INCLUDE 'stgp06.h'
 END INTERFACE

 bmatx = 0d0
 bsmtx = 0d0
 bmsig = 0d0
 IF( stype > 1 )gmatx = 0d0
 IF( ndofe > 5 )bfmtx = 0d0
 !IF(ndofe == 7 )THEN
 !  bsmtx(6:ndofe,:,:) = 0d0
 !  bmatx(1:5,:,9:14) = 0d0
 !END IF
 tgp = MATMUL(t,shape)  !director at Gauss points (not unit vector yet)
 !     evaluates first and second fundamental forms and distorsions
 CALL stran7(stype,cartd,dn,x,t,auxil,dx,dt,gamma,e,shape,ngaus,nnode,psi,psib,jas,ngamm) !=>auxil & gamma
 !     computes shear strains substracting the initial distorsion -->dstra
 !                computes bsshear matrix
 CALL bshem7(x,t,lb,bsmtx,bfmtx,jas,nangl,ngamm)
 gamma = gamma - gamm0     !sytpe == 1 not possible
 DO g=1,ngaus
   CALL vecuni(3,tgp(:,g),def)
   CALL vecpro(dx(:,1,g),dx(:,2,g),norml)
   ambda(2,g)= ambda(1,g)/DOT_PRODUCT(norml,tgp(1:3,g))
   ! computes strains substracting the initial configuration -->auxil
   IF(nlayr == 1 )THEN
     auxil(1:6,g) = auxil(1:6,g) - stra0(:,g)
     auxil(1:2,g) = auxil(1:2,g)/2d0
     IF( ndofe == 7 ) THEN    ! psi gradient
       auxil( 9,g) = DOT_PRODUCT(psi(1,:),cartd(:,1,g))
       auxil(10,g) = DOT_PRODUCT(psi(2,:),cartd(:,2,g))
       auxil(11,g) = DOT_PRODUCT(psi(1,:),cartd(:,2,g))
       auxil(12,g) = DOT_PRODUCT(psi(2,:),cartd(:,1,g))
       s(1:2) = MATMUL(ap1(:,:,g),psib)
       auxil(13:14,g) = MATMUL(jacin(:,:,g),s(1:2))
     END IF
   ELSE ! modify curvatures only for large-strain multilayer
     auxil(4:5,g) = (auxil(4:5,g) - stra0(4:5,g))*2d0  !twice the curvatures
     auxil(6,g) = auxil(6,g) - stra0(6,g)
   END IF
   ! evaluates shear strains at gauss integration points --> auxil
   s(1:2) = MATMUL(ap1(:,:,g),gamma)
   auxil(7:8,g) = MATMUL(jacin(:,:,g),s(1:2))
 END DO

 !   assumed strain for membrane
 IF (stype == 3) THEN
   DO g=1,3
     auxil(3,g) = auxil(3,g)/2.0d0
   END DO
   baux(1:nbeta) = 0d0
 END IF
 !     assumed strains for membrane
 IF(stype == 2 .OR. stype == 3) THEN
   !       interpolates the nodal strains
   DO g=1,ngaus
     CALL intem7(posgp(1,g),posgp(2,g),stype,e(:,g),auxil)
   END DO
   !                     transform from convected to cartesian
   DO g=1,ngaus
     CALL concar(e(:,g),auxil(:,g),jacin(:,:,g),1)
   END DO
   IF(stype == 2) beta = 0d0
 END IF

 ! evaluates total stresses at gauss points --> strsg
 IF( plast )THEN
   !strap(:,:,1) =  strap(:,:,2)   !initializes current values
   !ehist(1,:) = ehist(3,:)        !equivalent plastic strain
   ehist(2,:) = 0d0               !equivalent stress
   ehist(4:5,:) = 0d0             !consistent parameters
   IF(nlayr == 1)THEN
     ! computes elastic strains --> auxil
     auxil = auxil - strap
     DO g=1,ngaus
       CALL radre6(matty,dx(1,1,g),dx(1,2,g),ambda(1,g),             &
                   tgp(1,g),auxil(1,g),strsg(1,g),strap(1,g),        &
                   ehist(1,g),dmatx(1),props(1),istop)
     END DO
   ELSE
     CALL slayer(ngaus,nstre,strsg(1,1),auxil(1,1),dmatx(1),          &
                 props(1),strap(1,1),stres(1,1,1),ambda(1,1),nlayr,   &
                 thick,istop,ehist(1,1),ielem,plast)
   END IF

 ELSE     !elastic
   IF(nlayr == 1)THEN
     IF( ndofe == 5 )THEN
       CALL istgp6(ngaus,nstre,strsg(1,1),auxil(1,1),dmatx(1), &
                   ambda(1,1),matty,props(1),sf(1))
     ELSE      ! zigzag
       CALL stgp06(ngaus,strsg,auxil,dmatx)
     END IF
   ELSE
     CALL slayer(ngaus,nstre,strsg(1,1),auxil(1,1),dmatx(1),       &
                 props(1),s(1),s(2),ambda(1,1),nlayr,              &
                 thick,istop,s(2),ielem,plast)

   END IF
   DO g = 1,ngaus
     DO n=1,8
       energ(n) = energ(n) + dvolu(g)*strsg(n,g)*auxil(n,g)
     END DO
   END DO
 END IF
 !     initializes qaux
 !IF(stype /= 1) qaux = 0d0 !not possible
 qaux = 0d0
 qf   = 0d0
 !     sums bending contribution & membrane contribution
 !WRITE(55,"('sig',3e15.6,/,5E15.6)")strsg
 DO g=1,ngaus
   !       scales stresses for numerical integration
   auxil(:,g) = dvolu(g)*strsg(:,g)
   !       computes b matrix , bending & membrane
   CALL bmatx7(stype,cartd(:,:,g),lb,bmatx,dx(:,:,g),dt(:,:,g),ndofe,nangl)
   !       equivalent q for transverse shear
   s(1)= jacin(1,1,g)*auxil(7,g)+jacin(2,1,g)*auxil(8,g)
   s(2)= jacin(1,2,g)*auxil(7,g)+jacin(2,2,g)*auxil(8,g)
   xi(1)=1d0
   xi(2)=posgp(1,g)
   xi(3)=posgp(2,g)
   k=0
   DO i=1,2
     DO j=1,3
       k=k+1
       qaux(k)=qaux(k)+s(i)*xi(j)
     END DO
   END DO
   IF( ndofe == 7 )THEN
     CALL proma3(bmsig(1,1),bmatx(1,1,9),auxil(9,g),nvarl,4,.FALSE.)    !Psi moments
     s(1)= jacin(1,1,g)*auxil(13,g)+jacin(2,1,g)*auxil(14,g)
     s(2)= jacin(1,2,g)*auxil(13,g)+jacin(2,2,g)*auxil(14,g)
     xi(1)=1d0
     xi(2)=posgp(1,g)
     xi(3)=posgp(2,g)
     k=0
     DO i=1,2
       DO j=1,3
         k=k+1
         qf(k)=qf(k)+s(i)*xi(j)
       END DO
     END DO
   END IF
   !       equivalent beta for membrane stresses
   IF(stype == 2 .OR. stype == 3) THEN
     CALL proma3(bmsig(1,1),bmatx(1,1,4),auxil(4,g),nvarl,3,.FALSE.)
     CALL concar(s,auxil(:,g),jacin(:,:,g),4)
     k=0
     IF(stype == 2) THEN
       xi(1) = xi(1) - xi(2) - xi(3)
       DO i=1,3
         DO j=1,3
           k=k+1
           beta(k) = beta(k) + s(i)*xi(j)
         END DO
       END DO
     ELSE
       DO i=1,3
         DO j=1,3
           k=k+1
           baux(k) = baux(k) + s(i)*xi(j)
         END DO
       END DO
     END IF
   ELSE
     CALL proma3(bmsig(1,1),bmatx(1,1,1),auxil(1,g),nvarl,6,.FALSE.)
   END IF
 END DO  !Gauss point loop
 IF(stype == 2 .OR. stype == 3) THEN
   !                computes b matrix for assumed strain elements
   IF(stype == 2) THEN
     CALL bmmt27(nnode,x,gmatx,dn)
   ELSE
     !           multiplies by  the transpose of p**(-1)
     beta(1)= r3p1*baux(1) -r3*baux(2) -r3p3*baux(3) +r3p1*baux(7) &
             -r3*baux(8)-r3p3*baux(9)
     beta(2)=-r3m1*baux(1) +r3*baux(2) +r3m3*baux(3) -r3m1*baux(7) &
             +r3*baux(8)+r3m3*baux(9)
     beta(3)= 2d0*baux(7)-r3p3*2d0*baux(8)+r3m3*2d0*baux(9)
     beta(4)= 2d0*baux(7)+r3m3*2d0*baux(8)-r3p3*2d0*baux(9)
     beta(5)=-r3m1*baux(4)+r3m3*baux(5)+r3*baux(6)-r3m1*baux(7)    &
             +r3m3*baux(8)+r3*baux(9)
     beta(6)= r3p1*baux(4)-r3p3*baux(5)-r3*baux(6)+r3p1*baux(7)    &
             -r3p3*baux(8)-r3*baux(9)
     beta(7)= 3d0*baux(3)-3d0*baux(7)+3d0*baux(8)+6d0*baux(9)
     beta(8)= 3d0*baux(5)-3d0*baux(7)+6d0*baux(8)+3d0*baux(9)
     beta(9)= 3d0*baux(7)-3d0*baux(8)-3d0*baux(9)

     CALL bmmt37(nnode,x,gmatx,nfdas,dn(:,:,1))
   END IF

   DO n=1,nnode
     bmsig(1:3,n) = bmsig(1:3,n) + MATMUL(gmatx(:,n,:),beta)
   END DO

 END IF
 !       shear contribution
 !       evaluates q =  qaux* p**(-T)
 q(1) =  r3p1*qaux(1) - r3*qaux(2) - r3p1*qaux(3) - r3m1*qaux(5)
 q(2) = -r3m1*qaux(1) + r3*qaux(2) + r3m1*qaux(3) + r3p1*qaux(5)
 q(3) =  r3m1*qaux(3) + r3p1*qaux(5)
 q(4) = -r3p1*qaux(3) - r3m1*qaux(5)
 q(5) = -r3p1*qaux(3) + r3m1*qaux(4) - r3m1*qaux(5) - r3*qaux(6)
 q(6) =  r3m1*qaux(3) - r3p1*qaux(4) + r3p1*qaux(5) + r3*qaux(6)

 DO n=1,nnode
   bmsig(1:5,n) = bmsig(1:5,n) + MATMUL(bsmtx(1:5,n,:),q)
 END DO
 IF( ndofe == 7 )THEN
  !       evaluates q =  qf* p**(-1) * t
   q(1) =  r3p1*qf(1) - r3*qf(2) - r3p1*qf(3) - r3m1*qf(5)
   q(2) = -r3m1*qf(1) + r3*qf(2) + r3m1*qf(3) + r3p1*qf(5)
   q(3) =  r3m1*qf(3) + r3p1*qf(5)
   q(4) = -r3p1*qf(3) - r3m1*qf(5)
   q(5) = -r3p1*qf(3) + r3m1*qf(4) - r3m1*qf(5) - r3*qf(6)
   q(6) = +r3m1*qf(3) - r3p1*qf(4) + r3p1*qf(5) + r3*qf(6)

   DO n=1,nnode
     bmsig(6:7,n) = bmsig(6:7,n) + MATMUL(bfmtx(1:2,n,:),q)
   END DO
  END IF

 thick = thick*SUM(ambda(2,1:ngaus))/ngaus  !present average thickness
 RETURN
 END SUBROUTINE intrf7

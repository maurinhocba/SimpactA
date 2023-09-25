 SUBROUTINE intrf6(nvarl,plast,ndofe,nstre,                         &
                   matty,nlayr,dmatx,props,sf,x,t,lb,psi,nangl,jas, &
                   dvolu,cartd,jacin,stra0,gamm0,strsg,ambda,       &
                   ehist,strap,stres,bmsig,                         &
                   energ,thick,istop,ielem)
 !***********************************************************************
 !
 !****this routine compute the internal nodal forces for element no. 6
 !    4-node cuadrilateral shell element
 !***********************************************************************
 IMPLICIT NONE
 !***routine parameters
 INTEGER (kind=4), INTENT(IN) :: nvarl,matty,nlayr,ielem,ndofe,nstre
 INTEGER (kind=4), INTENT(OUT) :: istop
 REAL (kind=8), INTENT(IN) :: props(:),dmatx(:),x(:,:),t(:,:),   &
                lb(:,:,:),dvolu(:),cartd(:,:,:),psi(:,:),        &
                jacin(:,:,:),stra0(:,:),gamm0(:)
 REAL (kind=8), INTENT(IN) :: sf(3)

 REAL (kind=8), INTENT(IN OUT) :: thick,strsg(:,:),ambda(:,:),energ(:),    &
                bmsig(:,:)
 REAL(kind=8), POINTER :: ehist(:,:),strap(:,:),stres(:,:,:),nangl(:,:),jas(:,:)
 LOGICAL, INTENT(IN) :: plast
 !***local variables
 INTEGER (kind=4) n,g
 REAL    (kind=8)  bmatx(ndofe,nnode,nstre),bsmtx(ndofe,nnode,ngamm)
 REAL (kind=8) s(2),def,norml(3),tgp(3,ngaus),dstra(ngamm),                &
               auxil(nstre,ngaus),qaux(ngamm),psib(ngamm),                 &
               dx(3,2,ngaus),dt(3,2,ngaus)

 INTERFACE
   INCLUDE 'vecuni.h'
   INCLUDE 'stgp06.h'
 END INTERFACE

 bmatx = 0d0
 bsmtx = 0d0
 bmsig = 0d0
 tgp = MATMUL(t,shape)  !Director at Gauss points (not a unit vector)
 !     evaluates first and second fundamental forms and distorsions
 CALL stran6(ngaus,cartd,x,t,auxil,dx,dt) !==>auxil
 !                  evaluates distorsion at the assumed strain points
 CALL asstr6(ngamm,x,t,dstra)
 !     computes shear strains substracting the initial distorsion -->dstra
 dstra = dstra - gamm0
 IF( ndofe == 7 ) CALL psib06(psi,psib,jas)
 CALL bshem6(x,t,lb,bsmtx,ndofe,jas,nangl)
 DO g=1,ngaus  ! for each integration point
   CALL vecuni(3,tgp(:,g),def)             !unit vector
   CALL vecpro(dx(:,1,g),dx(:,2,g),norml)  !mid-surface normal
   ambda(2,g)= ambda(1,g)/DOT_PRODUCT(norml,tgp(1:3,g))  !thickness ratio
   ! computes strains substracting the initial configuration -->auxil
   IF(nlayr == 1 )THEN    !no TTT integration
     auxil(1:6,g) = auxil(1:6,g) - stra0(:,g)
     auxil(1:2,g) = auxil(1:2,g)/2d0           !Green Lagrange strains
     !----------------
     IF( ndofe == 7 ) THEN    ! psi gradient
       auxil( 9,g) = DOT_PRODUCT(psi(1,:),cartd(:,1,g))
       auxil(10,g) = DOT_PRODUCT(psi(2,:),cartd(:,2,g))
       auxil(11,g) = DOT_PRODUCT(psi(1,:),cartd(:,2,g))
       auxil(12,g) = DOT_PRODUCT(psi(2,:),cartd(:,1,g))
       s(1) = ((1d0-posgp(2,g))*psib(2)+(1d0+posgp(2,g))*psib(4))/2d0
       s(2) = ((1d0-posgp(1,g))*psib(1)+(1d0+posgp(1,g))*psib(3))/2d0
       auxil(13,g) = jacin(1,1,g)*s(1)+jacin(1,2,g)*s(2)
       auxil(14,g) = jacin(2,1,g)*s(1)+jacin(2,2,g)*s(2)
     END IF
   ELSE ! modify curvatures only for large-strain multilayer
     auxil(4:5,g) = (auxil(4:5,g) - stra0(4:5,g))*2d0  !twice the curvatures
     auxil(6,g) = auxil(6,g) - stra0(6,g)
   END IF
   ! evaluates shear strains at gauss integration points --> auxil
   s(1) = ((1d0-posgp(2,g))*dstra(2)+(1d0+posgp(2,g))*dstra(4))/2d0
   s(2) = ((1d0-posgp(1,g))*dstra(1)+(1d0+posgp(1,g))*dstra(3))/2d0
   auxil(7,g) = jacin(1,1,g)*s(1)+jacin(1,2,g)*s(2)
   auxil(8,g) = jacin(2,1,g)*s(1)+jacin(2,2,g)*s(2)
 END DO
 ! evaluates total stresses at gauss points --> strsg
 IF( plast )THEN
   ehist(2,:) = 0d0               !equivalent stress
   ehist(4:5,:) = 0d0             !consistent parameters
   IF(nlayr == 1)THEN
     ! computes elastic strains --> auxil
     auxil = auxil - strap
     DO g=1,ngaus
       CALL radre6(matty,dx(1,1,g),dx(1,2,g),ambda(1,g),             &
                   tgp(1,g),auxil(1,g),strsg(1,g),strap(1,g),      &
                   ehist(1,g),dmatx(1),props(1),istop)
     END DO
   ELSE
     CALL slayer(ngaus,nstre,strsg(1,1),auxil(1,1),dmatx(1),          &
                 props(1),strap(1,1),stres(1,1,1),ambda(1,1),nlayr, &
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
 qaux = 0d0
 IF( ndofe == 7 )psib = 0d0
 !     sums bending contribution & membrane contribution
 !WRITE(55,"('sig',3e15.6,/,5E15.6)",ERR=9999)strsg
 DO g=1,ngaus
   !       scales stresses for numerical integration
   auxil(:,g) = dvolu(g)*strsg(:,g)
   !       computes b matrix , bending & membrane
   CALL bmatx6(nnode,cartd(:,:,g),lb(:,:,:),bmatx(:,:,:),dx(:,:,g),dt(:,:,g),ndofe,nangl)
   CALL proma3(bmsig(1,1),bmatx(1,1,1),auxil(1,g),nvarl,6,.FALSE.)
   !       equivalent q for transverse shear
   s(1) = jacin(1,1,g)*auxil(7,g)+jacin(2,1,g)*auxil(8,g)
   s(2) = jacin(1,2,g)*auxil(7,g)+jacin(2,2,g)*auxil(8,g)
   qaux(1) = qaux(1) + (1d0-posgp(1,g))*s(2)/2d0
   qaux(2) = qaux(2) + (1d0-posgp(2,g))*s(1)/2d0
   qaux(3) = qaux(3) + (1d0+posgp(1,g))*s(2)/2d0
   qaux(4) = qaux(4) + (1d0+posgp(2,g))*s(1)/2d0
   IF( ndofe == 7 )THEN
     CALL proma3(bmsig(1,1),bmatx(1,1,9),auxil(9,g),nvarl,4,.FALSE.)
     s(1) = jacin(1,1,g)*auxil(13,g)+jacin(2,1,g)*auxil(14,g)
     s(2) = jacin(1,2,g)*auxil(13,g)+jacin(2,2,g)*auxil(14,g)
     psib(1) = psib(1) + (1d0-posgp(1,g))*s(2)/2d0
     psib(2) = psib(2) + (1d0-posgp(2,g))*s(1)/2d0
     psib(3) = psib(3) + (1d0+posgp(1,g))*s(2)/2d0
     psib(4) = psib(4) + (1d0+posgp(2,g))*s(1)/2d0
   END IF
 END DO
 !      add shear contributions
 DO n=1,nnode
   bmsig(1:5,n) = bmsig(1:5,n) + MATMUL(bsmtx(1:5,n,:),qaux)
   IF( ndofe == 7) bmsig(6:7,n) = bmsig(6:7,n) + MATMUL(bsmtx(6:7,n,:),psib)
 END DO

 thick = thick*SUM(ambda(2,1:ngaus))/ngaus  !present average thickness
 RETURN
 END SUBROUTINE intrf6

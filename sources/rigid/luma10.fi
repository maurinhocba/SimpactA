 SUBROUTINE luma10 (nelem,ntype,ndime,nnode,lnods,matno,coord,   &
                    emass,iwrit,sumat,tmass,lumped,mass,ifpre)
 !
 ! computes nodal mass for a set of rigid elements
 !
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem, & !number of elements
                                 ntype, & !type of rigid body
                                 ndime, & !problem dimension
                                 nnode, & !number of nodes per element
                                 iwrit, & !flag for echo
                                 lnods(:,:)  !connectivities
 INTEGER (kind=4), POINTER :: matno(:)      !properties
 REAL (kind=8), INTENT (IN) :: coord(:,:)      !nodal coordinates
 REAL (kind=8), INTENT (IN OUT) :: emass(:,:), & !nodal mass
                                   sumat,      & !total mass
                                   tmass         !total body mass
 LOGICAL :: lumped
 REAL (kind=8), INTENT (IN OUT) ::   mass(:)     !mass matrix
 INTEGER (kind=4), POINTER :: ifpre(:,:)         !Associated DOFs

 INTEGER (kind=4) :: ng,ngaus,ndiml,i,j,k,g
 REAL (kind=8) weigp(3),posgp(3),a,b
 REAL (kind=8), ALLOCATABLE :: shape(:,:),deriv(:,:,:),we(:),pos(:,:),grad(:,:)
 LOGICAL :: quad,flag

 REAL (kind=8), PARAMETER :: twopi=6.283185307179586

 TYPE(section), POINTER :: sec
 INTEGER (kind=4) isec,iel,osec
 REAL (kind=8) :: dens,thick,jac,sumas,r,sume1
 REAL (kind=8) :: diagm(nnode),x(ndime,nnode),t(ndime)
 LOGICAL :: newm

 INTERFACE
   INCLUDE 'shape9.h'
 END INTERFACE

 flag = ASSOCIATED( matno )     !Associated material
 quad = .FALSE.
 ! select variables
 SELECT CASE (ndime)
 CASE (2)                 !2 dimensional problems
   SELECT CASE (ntype)
   CASE (1,2,3)           !Solid elements
     ndiml = 2
     SELECT CASE (nnode)
     CASE (3)             !linear triangles
       ng    =-1
       ngaus = 1
     CASE (4)             !bilinear quads
       ng    = 2
       ngaus = 4
     CASE (6)             !cuadratic triangles
       ng    =-3
       ngaus = 3
       quad  = .TRUE.
     CASE (8)             !cuadratic quads
       ng    = 3
       ngaus = 9
       quad  = .TRUE.
     END SELECT
   CASE (4,5,6)           !Shell elements
     ndiml = 1
     SELECT CASE (nnode)
     CASE (2)             !linear
       ng    = 2
       ngaus = 2
     CASE (3)             !cuadratic
       ng    = 3
       ngaus = 3
       quad  = .TRUE.
     END SELECT
   END SELECT
 CASE (3)                 !3 dimensional problems
   SELECT CASE (ntype)
   CASE (1,2,3)           !Solid elements
     ndiml = 3
     SELECT CASE (nnode)
     CASE (4)             !Tetrahedra
       ng    =-1
       ngaus = 1
     CASE (6)             !Prisma
       ng    =-2
       ngaus = 2
     CASE (8)             !Hexahedra
       ng    = 2
       ngaus = 8
     END SELECT
   CASE (4,5,6)           !Shell elements
     ndiml = 2
     SELECT CASE (nnode)
     CASE (3)             !linear triangles
       ng    =-1
       ngaus = 1
     CASE (4)             !bilinear quads
       ng    = 2
       ngaus = 4
     CASE (6)             !cuadratic triangles
       ng    =-3
       ngaus = 3
       quad  = .TRUE.
     END SELECT
   END SELECT
 END SELECT

 ALLOCATE (shape(nnode,ngaus), deriv(nnode,ndiml,ngaus),we(ngaus))
 ALLOCATE( pos(ngaus,ndiml),grad(ndime,ndiml))
 IF( ng > 0 )THEN !Lagrangean Quads
   CALL gaussq (ng,posgp,weigp)
   SELECT CASE (ndiml)
   CASE (1)
     CALL shape9 (we,pos(:,1),shape,deriv(:,1,:),nnode,ngaus)
   CASE (2) !2-D
     g = 0
     DO i=1,ng
       DO j=1,ng
         g = g+1
         CALL shape3 (deriv(1,1,g) ,shape(1,g) ,posgp(i),posgp(j),nnode)
         we(g) = weigp(i)*weigp(j)
       END DO
     END DO
   CASE (3) !3-D
     g = 0
     DO i=1,ng
       DO j=1,ng
         DO k=1,ng
           g = g+1
           CALL shape5(deriv(1,1,g),shape(1,g),posgp(i),posgp(j),posgp(k),nnode)
           we(g) = weigp(i)*weigp(j)*weigp(k)
         END DO
       END DO
     END DO
   END SELECT
 ELSE             !triangles
   IF( ng == -1 )THEN        !linear triangles
     shape(:,:) = 1d0/REAL(nnode,8)
     IF( nnode == 3 )THEN      !triangles
       deriv = RESHAPE( (/ 1d0,0d0,-1d0,0d0,1d0,-1d0 /),(/ 3,2,1 /))
       we = 0.5d0
     ELSE                      !Tetrahedra
       deriv = RESHAPE(  (/ -1d0, 1d0, 0d0, 0d0,      &
                            -1d0, 0d0, 1d0, 0d0,      &
                            -1d0, 0d0, 0d0, 1d0/) , (/ 4,3,1 /))
       !deriv = RESHAPE( (/ 1d0, 0d0, 0d0,-1d0, 0d0, 1d0, 0d0,-1d0,   &
       !                    0d0, 0d0, 1d0,-1d0 /), (/ 4,3,1 /))
       we = 1d0/6d0
     END IF

   ELSE IF( ng == -2 )THEN  !prisma
     a = 1d0/3d0
     b = 1d0/SQRT(3d0)
     pos(1,:)= (/ a,a,-b /)
     pos(2,:)= (/ a,a,+b /)
     we = 0.5d0

     ! gauss points shape and derivatives of nodal functions

     DO g=1,ngaus
       CALL shape4(6,shape(1,g),deriv(1,1,g),pos(g,1),pos(g,2),pos(g,3))
     END DO

   ELSE IF( ng == -3 )THEN  !cuadratic triangles
     a = 1d0/6d0            !a = 1d0/6d0
     b = 2d0/3d0            !b = 1d0/2d0
     pos(1,:)= (/ a,a /)    !pos(1,:) = (/ b  ,0d0 /)
     pos(2,:)= (/ b,a /)    !pos(2,:) = (/ b  ,b   /)
     pos(3,:)= (/ a,b /)    !pos(3,:) = (/ 0d0,b   /)
     we = a

     ! gauss points shape and derivatives of nodal functions

     DO g=1,ngaus
       CALL shape7(pos(g,1),pos(g,1),shape(1,g),deriv(1,1,g))
     END DO
   END IF

 END IF

 ng = ABS(ng)
 sume1 = 0d0 !initializes set mass
 osec = -1
 DO iel=1,nelem
   IF( flag )THEN
     isec = matno(iel)                    !associated material
     newm = isec /= osec
     IF ( newm )THEN
       sec => psecs(isec)%p
       dens = sec%mtbas%prope(5)
       IF(ntype > 3 .OR. ntype == 1 )THEN !surfaces or plane-stress type
         thick  = sec%rprop(1)
         dens   = thick*dens              !surface density
       END IF
       osec = isec
     END IF
   ELSE
     dens = -1d0      !negative mass is a flag for scaling
   END IF
   ! ---computes jacobian at Gauss point
   x = coord(:,lnods(:,iel))
   sumas = 0d0 !element mass initialization
   diagm = 0d0 !nodal mass initialization
   DO g=1,ngaus
     grad = MATMUL(x,deriv(:,:,g))
     SELECT CASE (ndiml)
     CASE (1)
       jac = SQRT(DOT_PRODUCT(grad(:,1),grad(:,1)))
     CASE (2)
       IF( ndime == 2 )THEN
         jac = grad(1,1)*grad(2,2) - grad(1,2)*grad(2,1)
       ELSE
         CALL vecpro(grad(1,1),grad(1,2),t(1))
         jac = SQRT(DOT_PRODUCT(t,t))
       END IF
     CASE (3)
      jac =  grad(1,1)*grad(2,2)*grad(3,3)  &
            +grad(1,3)*grad(2,1)*grad(3,2)  &
            +grad(3,1)*grad(1,2)*grad(2,3)  &
            -grad(3,1)*grad(2,2)*grad(1,3)  &
            -grad(3,3)*grad(1,2)*grad(2,1)  &
            -grad(1,1)*grad(2,3)*grad(3,2)
     END SELECT
     IF( MOD(ntype,3)  == 0 ) THEN !axilsymmetric problems
       r = DOT_PRODUCT(x(1,:),shape(:,g))
       jac = jac*r*twopi
     END IF
     IF( jac <= 0d0 )THEN
       IF( iwrit > 0 )THEN
         WRITE(55,"()",ERR=9999)
       END IF
       CALL runen3('LUMA10: error in rigid body definition')
     END IF
     sumas = sumas + jac     !element mass
     ! ---generates lumped mass matrix proportional to diagonal
     DO j=1,nnode !for each node
       r = shape(j,g)
       IF( quad ) r = r*r             !quadratic elements
       diagm(j) = diagm(j) + jac*r
     END DO
     IF( quad )THEN                   !quadratic elements
       r = sumas/SUM(diagm)
       diagm = diagm*r
     END IF
   END DO
   ! ---add to global mass
   diagm = diagm*dens
   DO i=1,nnode                        !for each element node
     j = lnods(i,iel)             !node number
     emass(1:ndime,j) = emass(1:ndime,j) + diagm(i)
     !for consistent mass matrix use lumped matrix for rigid bodies
     IF( .NOT.lumped )THEN
       DO k=1,ndime
         CALL ensmat(1,ifpre(k,j),diagm,mass(1))
       END DO
     END IF
   END DO
   sume1 = sume1 + sumas*dens            !add to total mass

 END DO

 ! When Rigid body mass is given
 IF( tmass > 0d0 )THEN
   dens  = tmass/sume1
   DO iel=1,nelem
     DO i=1,nnode                 !for each element node
       j = lnods(i,iel)             !node number
       IF( emass(1,j) < 0d0 ) emass(1:ndime,j) = emass(1:ndime,j)*dens
     END DO
   END DO
   sume1 = tmass
 END IF
 sumat = sumat + sume1
 IF( iwrit > 0 )  &
 WRITE(lures,"(//'  total mass for rigid body :', e15.7)",ERR=9999) sume1

 DEALLOCATE (shape,deriv,we)
 DEALLOCATE (pos,grad)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE luma10

 SUBROUTINE luma10p (ndime,nnode,lnods,emass,sumat,tmass)
 !
 ! computes nodal mass for a set of rigid node
 !
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: ndime, & !problem dimension
                                 nnode, & !number of nodes per element
                                 lnods(:,:)    !connectivities
 REAL (kind=8), INTENT (IN) :: tmass           !total mass
 REAL (kind=8), INTENT (IN OUT) :: emass(:,:), & !nodal mass
                                   sumat         !total mass
 ! local variables
 REAL (kind=4) :: nmass
 INTEGER (kind=4) :: i,n


 sumat = sumat + tmass

 nmass = tmass/nnode

 DO i=1,nmass
   n = lnods(i,1)
   emass(1:ndime,n ) = emass(1:ndime, n) + nmass
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE luma10p

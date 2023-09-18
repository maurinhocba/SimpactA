 SUBROUTINE gaus10(nelem, ntype, ngaus, nnode, lnods, coord,         &
                   cartd, dvolu, shape, istop)
!********************************************************************
!
!*** evaluates gauss points constants for RIGID-HEAT elements
!
!********************************************************************
 USE ctrl_db, ONLY : ndime
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: nelem, & !number of elements
                                 ntype, & !type of rigid body
                                 nnode, & !number of nodes per element
                                 lnods(:,:)  !connectivities
 REAL (kind=8), INTENT (IN) :: coord(:,:)      !nodal coordinates
 INTEGER (kind=4), INTENT(OUT) :: ngaus,istop
 REAL    (kind=8), POINTER :: cartd(:,:,:,:),dvolu(:,:),shape(:,:)

 ! Local variables
 INTEGER (kind=4) :: ng,i,j,k,g,iel
 REAL (kind=8) weigp(3),posgp(3),a,b
 REAL (kind=8), ALLOCATABLE :: deriv(:,:,:),we(:),pos(:,:)

 REAL (kind=8), PARAMETER :: twopi=6.283185307179586

 REAL (kind=8) :: jac,r
 REAL (kind=8) :: x(ndime,nnode)

 ! select variables
 SELECT CASE (ndime)
 CASE (2)                 !2 dimensional problems
   SELECT CASE (nnode)
   CASE (3)             !linear triangles
     ng    =-1
     ngaus = 1
   CASE (4)             !bilinear quads
     ng    = 2
     ngaus = 4
   END SELECT
 CASE (3)                 !3 dimensional problems
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
 END SELECT

 ALLOCATE ( deriv(nnode,ndime,ngaus),we(ngaus), pos(ngaus,ndime) )
 ALLOCATE( shape(nnode,ngaus),cartd(nnode,ndime,ngaus,nelem),dvolu(ngaus,nelem))

 IF( ng > 0 )THEN !Lagrangean Quads
   CALL gaussq (ng,posgp,weigp)
   SELECT CASE (ndime)
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
       deriv = RESHAPE( (/ 1d0, 0d0, 0d0,-1d0, 0d0, 1d0, 0d0,-1d0,   &
                           0d0, 0d0, 1d0,-1d0 /), (/ 4,3,1 /))
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

   END IF

 END IF

 ng = ABS(ng)

 DO iel=1,nelem
   ! ---computes jacobian at Gauss point
   x = coord(:,lnods(:,iel))
   DO g=1,ngaus
     SELECT CASE (ndime)
     CASE (2)
       CALL jacob3(cartd(1,1,g,iel),deriv(1,1,g),jac,x,iel,nnode,istop)
       IF( ntype  == 3 ) THEN !axilsymmetric problems
         r = DOT_PRODUCT(x(1,:),shape(:,g))
         jac = jac*r*twopi
       END IF
     CASE (3)
       CALL jacob5(cartd(1,1,g,iel),deriv(1,1,g),jac,x,nnode,iel,istop)
     END SELECT
     IF( jac <= 0d0 )THEN
       CALL runen3('GAUS10: error in rigid body definition')
     END IF
     dvolu(g,iel) = jac*we(g) !does not include thickness
   END DO
 END DO

 DEALLOCATE (deriv,we,pos)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE gaus10

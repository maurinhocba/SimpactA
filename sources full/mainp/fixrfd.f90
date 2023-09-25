 SUBROUTINE fixrfd(xa)
 !***********************************************************************
 !
 !     updates configuration and factors of rotation-free dependant displacements
 !
 !***********************************************************************
 USE ctrl_db, ONLY : ndime
 USE kinc_db, ONLY : npsdf,nndpd,ftsdf
 USE npo_db,  ONLY : ifpre
 USE rfdp_db
 IMPLICIT NONE
 REAL (kind=8),INTENT(IN OUT) :: xa(:,:)

 INTEGER (kind=4) :: i,j,jj,k,l,m,sl
 REAL    (kind=8) :: ll,s,f1,f2,ld,t(ndime),n(ndime),x(ndime,6),xc(ndime),ls(3,2)
 REAL    (kind=8) :: g(ndime,3),tn(2,2),p(3,3),pa(3,3,2),nfd(6,3)
 TYPE (rfdp_nod), POINTER :: rfd
 INTEGER (kind=4), PARAMETER ::  &  !side element connectivities for triangles
           hh(2,3) = RESHAPE((/ 3,2, 1,3, 2,1 /), (/2,3/) )


 rfd => rfd_head       !point to first slave-master pair
 DO m=1,nrfdp          !loop for each pair
   sl = rfd%slave         !slave node

   IF( ndime == 2 )THEN   ! T W O  D I M E N S I O N A L   P R O B L E M S

     s  = rfd%rfpar(1)      !local coordinate
     ld = rfd%rfpar(2)      !original distance

     IF( planar ) THEN      ! 2 - NODE PLANAR APPROACH

       t = (xa(:,rfd%lnods(2)) - xa(:,rfd%lnods(1)))   !tangent vector between nodes (master)
       CALL vecuni(2,t,ll)                 ! ll : present length
       n = (/ -t(2),t(1) /)                !segment normal
       f1 = s                              !factor associated to first node
       f2 = 1d0-s                          !factor associated to second node
       xa(:,sl) = xa(:,rfd%lnods(1))*f1 + xa(:,rfd%lnods(2))*f2 + ld*n !updated slave coordinates
       ld = -ld/ll                         !factor
       tn = RESHAPE( (/ t(1)* n(1),  t(2)* n(1),  t(1)* n(2), t(2)* n(2)/),(/2,2/))*ld
       DO i=1,2  !for each slave DOF
         j  = -ifpre(i,sl)                 !position in the pointer array of first DOF
         k  = npsdf(j)                     !pointer to first DOF in ftsdf
         !new factors
         ftsdf(k)     =  - tn(i,1)
         ftsdf(k+1)   =  - tn(i,2)
         ftsdf(k+2)   =  + tn(i,1)
         ftsdf(k+3)   =  + tn(i,2)
         j = k+i-1
         ftsdf(j)   = ftsdf(j)   + f1
         ftsdf(j+2) = ftsdf(j+2) + f2
       END DO

     ELSE                 ! 3 - NODE CURVED APPROACH

       ls(1,1) = (s**2-s)/2d0                 !first local function
       ls(2,1) = (1d0-s**2)                   !second local function
       ls(3,1) = (s**2+s)/2d0                 !third local function
       ls(1,2) = s-5d-1                       !first local function derivative
       ls(2,2) = -2d0*s                       !second local function derivative
       ls(3,2) = s+5d-1                       !third local function derivative
       x(:,1:3) = xa(:,rfd%lnods)             !gather actual coordinates
       t = MATMUL(x(:,1:3),ls(:,2))           !tangent vector
       CALL vecuni(2,t,s)                     !unit tangent vector
       n = (/ -t(2),t(1) /)                   !nit normal vector
       xa(:,sl) = MATMUL(x(:,1:3),ls(:,1)) + ld*n  !updated slave coordinates
       ! compute proyection matrix
       tn = RESHAPE( (/ t(1)* n(1),  t(2)* n(1),  t(1)* n(2), t(2)* n(2)/),(/2,2/))*ld/s
       !COMPUTE new factors
       DO i=1,2  ! for each DOF
         k  = npsdf(-ifpre(i,sl))            !pointer to first DOF in ftsdf
         DO j=1,3                            !for each master node
           DO l=1,2                            !for each master DOF
             ftsdf(k) = - tn(i,l)*ls(j,2)      ! -t x n  N,xi
             IF( l == i ) ftsdf(k) = ftsdf(k) + ls(j,1) !include shape function
             k = k+1
           END DO
         END DO
       END DO
     END IF

   ELSE !  T H R E E   D I M E N S I O N A L   P R O B L E M S

     x(:,1:3) = xa(:,rfd%lnods(1:3))                 !3-node triangle coordinates
     ld = rfd%rfpar(3)                               !original distance

     IF( planar ) THEN  !  3-NODE PLANAR APPROACH

       xc(1:2) = rfd%rfpar(1:2)                      !local triangular coordinates
       xc(3) = 1d0 - xc(1) - xc(2)                   !third triangular coordinate
       DO i=1,3  !for each side
         g(1:3,i) = x(:,hh(1,i)) - x(:,hh(2,i))      ! lI
       END DO
       CALL vecpro(g(1,1),g(1,2),t(1))               ! t 2a
       CALL vecuni(3,t(1),ll)                        ! t  & ll=2A
       xa(:,sl) = MATMUL(x(:,1:3),xc) + ld*t         !update slave position
       ! compute side derivatives
       DO i=1,3
         CALL vecpro(g(1,i),t(1),p(1,i))             !normal direction times side length
       END DO
       ! compute slave factors for each node
       t = ld*t/ll                          ! d t /2A
       DO i=1,3                             !for each slave DOF
         jj = -ifpre(i,sl)                    !slave DOF order
         k  = npsdf(jj)                       !pointer to first slave DOFs
         DO j=1,3                             !for each master node
           DO l=1,3                             !for each DOF in the master node
             ftsdf(k) = p(i,j)*t(l)                   ! n,i t d /h
             IF( l == i ) ftsdf(k) = ftsdf(k) + xc(j) !include local coordinate
             k = k + 1                                !update master DOF
           END DO
         END DO
       END DO

     ELSE           ! 6-NODE CURVED APPROACH


       f1 = rfd%rfpar(1)       !xita
       f2 = rfd%rfpar(2)       !eta
       s  = 1d0 - f1 - f2      !zeta
       ! nodal functions and functions derivatives
       nfd(:,1) = (/s+f1*f2,f1+f2*s,f2+s*f1,s*(s-1d0)/2d0,f1*(f1-1d0)/2d0,f2*(f2-1d0)/2d0 /)
       nfd(:,2) = (/-1d0+f2, 1d0-f2,   s-f1,      -s+5d-1,        f1-5d-1,      0d0 /)
       nfd(:,3) = (/-1d0+f1,   s-f2, 1d0-f1,      -s+5d-1,            0d0,    f2-5d-1 /)
       ! get additional coordinates
       DO i=4,6
         j = i-3                                 !side
         k = rfd%lnods(i)                        !node
         IF( k > 0 )THEN                         !if side exist
           x(:,i) = xa(:,k)                           !additional node coordinates
         ELSE
           x(:,i) = 0d0                               !null coordinates
           nfd(j,:) = nfd(j,:) - nfd(i,:)             !include auxiliar node factors
           nfd(hh(1,j),:) = nfd(hh(1,j),:) + nfd(i,:) !in shape functions
           nfd(hh(2,j),:) = nfd(hh(2,j),:) + nfd(i,:) !and shape functions derivatives
           nfd(i,:) = 0d0                             !null node contribution
         END IF
       END DO

       g = MATMUL(x,nfd)                         !proyection and derivatives
       CALL vecpro(g(1,2),g(1,3),t(1))           !normal vector
       CALL vecuni(3,t(1),ll)                    !unit vector and twice the area
       xa(:,sl) = g(:,1) + ld*t                  !update slave node position

       ! compute proyection matrix 1 - t x t
       DO i=1,3
         DO j=1,3
           p(i,j) = -t(i)*t(j)
         END DO
         p(i,i) = p(i,i) + 1d0
       END DO
       p = p * ld/ll     ! d/2A (1- t x t)
       DO l=1,2          !compute curl matrices
         j = l+1         !correct pointer to derivatives
         pa(:,:,l) = RESHAPE((/ 0d0,g(3,j),-g(2,j),-g(3,j),0d0,g(1,j),g(2,j),-g(1,j),0d0 /),(/3,3/))
         pa(:,:,l) = MATMUL(p,pa(:,:,l))  !d/h (1 - t x t)*( X,xi x )
       END DO
       ! compute slave factors for each node
       DO i=1,3                                !for each slave DOF
         k  = npsdf(-ifpre(i,sl))              !pointer to first slave DOFs
         DO j=1,6                                !for each master node
           IF( rfd%lnods(j) <= 0 )CYCLE          !skip auxiliar nodes
           DO l=1,3                              !for each DOF in the master node
             ftsdf(k) = pa(i,l,1)*nfd(j,3) - pa(i,l,2)*nfd(j,2)  !normal variation
             IF( l == i ) ftsdf(k) = ftsdf(k) + nfd(j,1)         !include shape function
             k = k + 1
           END DO
         END DO
       END DO
     END IF
   END IF
   rfd => rfd%next                        !point to next pair
 END DO
 RETURN

 END SUBROUTINE fixrfd

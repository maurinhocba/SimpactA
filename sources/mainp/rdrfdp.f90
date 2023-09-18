 SUBROUTINE rdrfdp(npsdf,nesdf,ftsdf,jj,np)
 !**********************************************************************
 !
 !     Generates Data for Rotation-free Dependant Nodes
 !     Compute local coordinates
 !     Dependencies
 !     Initial factors
 !
 !**********************************************************************
 USE lispa0, ONLY : lures
 USE rfdp_db
 USE ctrl_db, ONLY : ndime,ndofn,npoin
 USE npo_db,  ONLY : ifpre,coord,label
 USE outp_db, ONLY: iwrit

 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4),INTENT(IN) :: np   !maximum number of dependencies (to check)
 INTEGER (kind=4),INTENT(IN OUT) :: jj, &        !present number of slave nodes
                                    npsdf(np), & !pointer to NESDF & FTSDF
                                    nesdf(2*np)  !master DOFs
 REAL (kind=8), INTENT(IN OUT) :: ftsdf(2*np)    !factors

 ! local variables
 INTEGER (kind=4) :: i,j,k,l,n,sl,ln(6)
 REAL (kind=8) :: l0,r,s,f,a2,dd,z
 REAL (kind=8) :: t(ndime),d(ndime),x(ndime,6),lc(ndime),xa(ndime),ls(3,3),lr(ndime,2), &
                  tn(ndime,ndime),p(3,3),pa(3,3,2),nfd(6,3)
 TYPE (rfdp_nod), POINTER :: rfd
 LOGICAL :: sides(3)
 INTEGER (kind=4), PARAMETER ::  &  !side element connectivities for triangles
           hh(2,3) = RESHAPE((/ 3,2, 1,3, 2,1 /), (/2,3/) )
 REAL (kind=8), PARAMETER :: tol = 0.0001


 ! print heading
 IF (iwrit == 1) WRITE(lures,"(//,3x,'Slave/Master nodes Data',//&
                 & '   Slave Node  Master Nodes')",ERR=9999)

 rfd => rfd_head       !point to first slave-master pair
 DO n=1,nrfdp          !loop for each dependant node

   sl = rfd%slave               !slave internal node
   DO i = 1,ndime               !for each (slave) DOF
     IF(ifpre(i,sl) < 0 ) THEN  !check that the DOF is not already slave
       WRITE(lures,"(' DOF already declared Slave:  ',&
                  &  'Node',i6,' DOF',i2)",ERR=9999) label(sl) ,i
       CALL runend('RDRFPD: Inconsistent Input Data    ')
     END IF
   END DO

   IF( ndime == 2 )THEN       ! T W O  D I M E N S I O N A L   P R O B L E M S

     IF( planar ) THEN ! 2 -NODE PLANAR APPROXIMATION

       ln(1:2) = rfd%lnods(1:2)               !nodes of the master element
       t = coord(:,ln(2)) - coord(:,ln(1))    !distance vector between nodes (master)
       CALL vecuni(ndime,t,l0)                !unit tangent vector
       d  = coord(:,sl) - coord(:,ln(1))      !local coordinates of slave node
       s = 2d0*DOT_PRODUCT(t,d)/l0-1d0        !proyected local component [-1:1]
       IF( ABS(s) > 1d0 ) WRITE(55,"(' Warning slave node ',i5, &
                                   & ' proyects out of the segment',f7.4)")label(sl),s
       rfd%rfpar(2) = -t(2)*d(1)+t(1)*d(2)    !distance from slave node to segment
       lc(1) = (1d0-s)/2d0                    !first local function
       lc(2) = 1d0-lc(1)                      !second local function
       rfd%rfpar(1) = lc(1)                   !keep first nodal function [0:1]
       f  = rfd%rfpar(2)/l0                   !d/l0
       d  = (/ -t(2),t(1) /)*f                !normal vector (scaled)
       tn = -RESHAPE((/ t(1)*d(1), t(2)*d(1),t(1)*d(1),t(2)*d(2) /), (/2,2/))
       DO i=1,2       !for each slave DOF
         k = npsdf(jj)                        !pointer to first slave DOFs
         ifpre(i,sl) = -jj                    !slave DOF order
         jj = jj+1                            !increase counter of slave DOFs
         npsdf(jj)  = k + 4                   !compute next pointer

         f = -1d0                             !factor
         DO j=1,2
           DO l=1,2
             nesdf(k)   = l+10*ln(j)          !store master node and dof until later
             ftsdf(k)   = tn(i,l) * f
             IF(i == l) ftsdf(k) = ftsdf(k) + lc(j)
             k = k+1
           END DO
           f = -f
         END DO

       END DO

     ELSE  ! 3 NODE CURVED APPROXIMATION

       s = 0  !initializes proyection on central node
       IF( rfd%lnods(1) == -1 )THEN    !check flag to get actual connectivities
         ln(1:3) = rfd%lnods(1:3)      !flag, nset, iel
         CALL point_ele11e(ln)         !get 4-node patch
         t = coord(:,ln(3)) - coord(:,ln(2))      !distance vector between nodes (master)
         CALL vecuni(2,t,l0)                      !unit tangent vector
         d  = coord(:,sl) - coord(:,ln(2))        !local coordinates of slave node
         s = 2d0*DOT_PRODUCT(t,d)/l0-1d0          !proyected local component [-1:1]
         IF( s < 0d0 )THEN        ! left node is nearer
           IF(ln(1) > 0 )THEN     ! if left node exists
             rfd%lnods = ln(1:3)      !use first three nodes
             s = s + 1d0              !shift local coordinate
           ELSE
             rfd%lnods(1:3) = ln(2:4) !shift to the right
             s = s - 1d0              !shift coordinate
             WRITE(55,"(' Warning slave node ',i5,' near a border of the patch',f7.4)")label(sl),s
           END IF
         ELSE                     ! right node is nearer
           IF(ln(4) > 0 )THEN     ! if right node exist
             rfd%lnods = ln(2:4)      !use the last three nodes
             s = s - 1d0              !shift coordinate
           ELSE
             rfd%lnods(1:3) = ln(1:3) !shift to the left
             s = s + 1d0              !shift coordinate
             WRITE(55,"(' Warning slave node ',i5,' near a border of the patch',f7.4)")label(sl),s
           END IF
         END IF
       END IF
       x(:,1:3) = coord(:,rfd%lnods(1:3))       !recover first two node coordinates
       lr(:,1) = (x(:,3) - x(:,1))/2d0          !average tangent vector
       lr(:,2) = x(:,1) - 2d0*x(:,2) + x(:,3)   !cuvature factor
       DO ! iterate until convergence
         xa = x(:,2) + s*lr(:,1) + s**2/2d0*lr(:,2)  !surface point
         d  = coord(:,sl) - xa             !distance
         t = lr(:,1)+s*lr(:,2)             !tangent vector
         l0 = DOT_PRODUCT(t,t)             !square length
         f = + DOT_PRODUCT(d,t)/l0         !correction
         s = s + f                         !update local coordinate
         IF( ABS(f) < tol )EXIT         !check convergence
       END DO
       l0 = SQRT(l0)                       !length
       rfd%rfpar(2) = (-t(2)*d(1)+t(1)*d(2))/l0     !distance from slave node to segment
       ls(1,1) = (s**2-s)/2d0                 !first local function
       ls(2,1) = (1d0-s**2)                   !second local function
       ls(3,1) = (s**2+s)/2d0                 !third local function
       ls(1,2) = s-5d-1                       !first local function derivative
       ls(2,2) = -2d0*s                       !second local function derivative
       ls(3,2) = s+5d-1                       !third local function derivative
       rfd%rfpar(1) = s                       !keep local coordinate
       d  = (/ -t(2),t(1) /)*rfd%rfpar(2)/l0**3 !normal vector (scaled)
       tn = -RESHAPE((/ t(1)*d(1), t(2)*d(1),t(1)*d(2),t(2)*d(2) /), (/2,2/))
       DO i=1,2       !for each slave DOF
         !WRITE(55,"('CONNOD',i5,' CONDOF ', i3)")label(sl),i
         k = npsdf(jj)                        !pointer to first slave DOFs
         ifpre(i,sl) = -jj                    !slave DOF order
         jj = jj+1                            !increase counter of slave DOFs
         npsdf(jj)  = k + 6                   !compute next pointer

         DO j=1,3                          !for each master node
           DO l=1,2                          !for each DOF
             nesdf(k)   = l+10*rfd%lnods(j)    !store master node and dof until later
             ftsdf(k)   = tn(i,l) * ls(j,2)            !factor (n x t)d N,xi
             IF(i == l) ftsdf(k) = ftsdf(k) + ls(j,1)  !add nodal function
             !IF( ABS(ftsdf(k)) > 1e-6 )WRITE(55,"('INDNOD ',i5,' INDDOF', i3,' FACTOR ',e15.5)")&
             !                                   label(rfd%lnods(j)),l,ftsdf(k)
             k = k+1
           END DO
         END DO

       END DO

     END IF

   ELSE         ! T H R E E   D I M E N S I O N A L   P R O B L E M S

     x(:,1:3)  = coord(:,rfd%lnods(1:3))                    !triangle coordinates

     IF( planar )THEN   ! 3-NODE PLANAR APPROXIMATION

       DO i=1,3         ! compute sides
         ls(1:3,i) = x(:,hh(1,i)) - x(:,hh(2,i))       !side lI
       END DO
       CALL vecpro(ls(1,1),ls(1,2),t(1))               ! normal t times 2A
       CALL VECUNI(3,t(1),a2)                          ! unit vector t  2A
       d = coord(:,sl) - x(:,1)                        ! distance vector of slave node to element center
       dd = DOT_PRODUCT(d,t)                           ! distance to the plane
       rfd%rfpar(3)   = dd                             ! keep local normal coordinate d
       d = d - dd*t                                    ! in plane vector
       CALL vecpro(ls(1,3),d(1),xa(1))                 ! area3
       lc(3) = DOT_PRODUCT(xa,t)/a2                    ! area coordinate 3
       CALL vecpro(ls(1,2),d(1),xa(1))                 ! area2
       lc(2) = DOT_PRODUCT(xa,t)/a2                    ! area coordinate 2
       lc(1) = 1d0 - lc(2) - lc(3)                     ! area coordinate 1
       rfd%rfpar(1:2) = lc(1:2)                        ! keep local coordinates (xita,eta)
       IF( ANY (lc < 0d0 ) .OR. ANY (lc > 1d0 ))  &
             WRITE(55,"(' Warning slave node ',i5,' proyects out of the segment',3f7.4)")label(sl),lc

       ! compute side derivatives
       DO i=1,3
         CALL vecpro(ls(1,i),t(1),tn(1,i))    !normal direction x side lenght (outward)
       END DO                                 ! when divided by 2A = side derivative 1/h
       d = dd/a2 * t                          ! normal vector t   x d / twice the area
       DO i=1,3                               !for each slave DOF
         k  = npsdf(jj)                       !pointer to first slave DOFs
         ifpre(i,sl) = -jj                    !slave DOF order
         !WRITE(55,"('CONNOD',i5,' CONDOF ', i3)")label(sl),i
         DO j=1,3                             !for each master node
           DO l=1,3                           !for each DOF in the master node
             nesdf(k) = l+10*rfd%lnods(j)     !store master node and dof until later
             ftsdf(k) = tn(i,j)*d(l)          !1/h  n,i t
             IF( l == i ) ftsdf(k) = ftsdf(k) + lc(j) !include local coordinate
             !IF( ABS(ftsdf(k)) > 1e-6 )WRITE(55,"('INDNOD ',i5,' INDDOF', i3, ' FACTOR',e15.5)")&
             !                                     label(rfd%lnods(j)),l,ftsdf(k)
             k = k + 1
           END DO
         END DO
         jj = jj+1                            !increase counter of slave DOFs
         npsdf(jj)  = k                       !compute next pointer to NESDF & FTSDF
       END DO

     ELSE  ! 6-NODE CURVED APPROXIMATION

       IF(rfd%lnods(4) == -1 )CALL point_ele13e(rfd%lnods) !get neighbour elements

       DO i=4,6  !for each extra node

         l = rfd%lnods(i)         !node number
         j = i-3                  !side
         IF( l > 0 )THEN          !if node exists
           x(:,i) = coord(:,l)        !node coordinates
           sides(j) = .TRUE.          !side exists
         ELSE                     !node does not exist
           x(:,i) = -x(:,i-3)+x(:,hh(1,i-3))+x(:,hh(2,i-3))  !auxiliar node
           sides(j) = .FALSE.         !side does not exist
         END IF
       END DO
       r = 1d0/3d0; s = r; z = r      !initializes coordinates (element center)
       DO                  !loop to compute proyection
         ! shape functions and shape function derivatives
         nfd(:,1) = (/ z+r*s, r+s*z, s+z*r, z*(z-1d0)/2d0, r*(r-1d0)/2d0, s*(s-1d0)/2d0 /)
         nfd(:,2) = (/-1d0+s, 1d0-s,   z-r,       -z+5d-1,        r-5d-1,           0d0 /)
         nfd(:,3) = (/-1d0+r,   z-s, 1d0-r,       -z+5d-1,           0d0,        s-5d-1 /)
         xa = MATMUL(x,nfd(:,1))                     !proyection
         ls(:,1:2) = MATMUL(x,nfd(:,2:3))            !tangent vectors
         CALL vecpro(ls(1,1),ls(1,2),t(1))           ! normal t times 2A
         CALL vecuni(3,t(1),a2)                      ! unit vector t & twice the triangle area
         d = coord(:,sl) - xa                        ! distance vector of slave node to proyection
         CALL vecpro(ls(1,2),t(1),lr(1,1))           ! X,eta x t
         lr(:,1) = lr(:,1)/a2                        ! first contravariant vector
         CALL vecpro(t(1),ls(1,1),lr(1,2))           ! t x X,xi
         lr(:,2) = lr(:,2)/a2                        ! second contravariant vector
         f = DOT_PRODUCT(d,lr(:,1))                  ! increment in xita
         z = DOT_PRODUCT(d,lr(:,2))                  ! increment in eta
         r = r + f                                   ! update xi
         s = s + z                                   ! update eta
         IF( ABS(f) + ABS(z) < tol )EXIT          ! check convergence
         z = 1d0 - r - s                             ! update z
       END DO
       rfd%rfpar(1) = r                              ! keep local coordinates (xita)
       rfd%rfpar(2) = s                              ! keep local coordinates (eta)
       z = 1d0 - r - s                               ! update z to check
       IF( r < 0d0 .OR. s < 0d0 .OR. z < 0d0 .OR. r > 1d0 .OR. s > 1d0 .OR. z > 1d0 ) &
             WRITE(55,"(' Warning slave node ',i5,' proyect out of the segment',3f7.4)")label(sl),r,s,z
       dd = DOT_PRODUCT(d,t)                         ! distance to the plane
       rfd%rfpar(3) = dd                             ! keep local normal coordinate d
       ! compute proyection matrix      1 - t x t
       DO i=1,3
         DO j=1,3
           p(i,j) = -t(i)*t(j)
         END DO
         p(i,i) = p(i,i) + 1d0
         IF (sides(i) ) CYCLE        !node exists
         nfd(i,:) = nfd(i,:) - nfd(i+3,:)             !include auxiliar node factors
         nfd(hh(1,i),:) = nfd(hh(1,i),:) + nfd(i+3,:) !on shape functions
         nfd(hh(2,i),:) = nfd(hh(2,i),:) + nfd(i+3,:) !and shape functions derivatives
       END DO
       p = p * dd/a2
       DO l=1,2      !compute curl matrices
         pa(:,:,l) = RESHAPE((/ 0d0,ls(3,l),-ls(2,l),-ls(3,l),0d0,ls(1,l),ls(2,l),-ls(1,l),0d0 /),(/3,3/))
         pa(:,:,l) = MATMUL(p,pa(:,:,l))
       END DO
       DO i=1,3                                !for each slave DOF
         k  = npsdf(jj)                        !pointer to first slave DOFs
         ifpre(i,sl) = -jj                     !slave DOF order
         !WRITE(55,"('CONNOD',i5,' CONDOF ', i3)")label(sl),i
         DO j=1,6                                !for each master node
           IF( rfd%lnods(j) <= 0 )CYCLE          !skip auxiliar node
           DO l=1,3                              !for each DOF in the master node
             nesdf(k) = l+10*rfd%lnods(j)          !store master node and dof until later
             ftsdf(k) = pa(i,l,1)*nfd(j,3) - pa(i,l,2)*nfd(j,2)  ! normal variation factors
             IF( l == i ) ftsdf(k) = ftsdf(k) + nfd(j,1)         !include shape function
             !IF( ABS(ftsdf(k)) > 1e-6 )WRITE(55,"('INDNOD ',i5,' INDDOF', i3, ' FACTOR',e15.5)")&
             !           label(rfd%lnods(j)),l,ftsdf(k)
             k = k + 1
           END DO
         END DO
         jj = jj+1                             !increase counter of slave DOFs
         npsdf(jj)  = k                        !compute next pointer to NESDF & FTSDF
       END DO
     END IF
   END IF
   IF (iwrit == 1) THEN
     ln = 0
     j = SIZE(rfd%lnods)
     DO i=1,j
       IF( rfd%lnods(i) > 0 )ln(i) = label(rfd%lnods(i))
     END DO
     WRITE(lures, '(i10,6i8)',ERR=9999) label(sl), ln(1:j)   !ECHO
   END IF
   rfd => rfd%next
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rdrfdp

 SUBROUTINE slum16 (nelem,heade,ymass,coord,dtscal,btscal,sms,alp,thl)

 ! compute lumped mass vector for element 16 (3-D solid TLF)
 ! that leads to a unitary time increment
 ! use standard 2 integration rule

 IMPLICIT NONE
 !  dummy arguments

 INTEGER (kind=4), INTENT (IN) :: nelem     !number of elements in the set
 REAL (kind=8), INTENT(IN) :: dtscal, &
                              btscal        !critical time factor
 REAL (kind=8), INTENT(IN OUT) :: alp, thl      !alpha and layer thickness for SMS
 REAL (kind=8) :: ymass(:,:), &    !nodal mass in each direction
                  coord(:,:)       !nodal coordinates
 LOGICAL :: sms                             !selective mass scaling is active for this set
 TYPE (ele16), POINTER  :: heade   !pointer to first element

 ! local variables
 INTEGER, PARAMETER  :: nvert = 9, ngaud = 2
 INTEGER, PARAMETER  :: h(2,nvert) = (/ 1,2, 2,3, 3,1, 4,5, 5,6, 6,4, 1,4, 2,5, 3,6 /)

 REAL (kind=8) :: shape(nnode,ngaud),deriv(nnode,3,ngaud),wei(ngaud), &
                  posgp(ngaud),weigp(ngaud),cd(nnode,3),x(3,nnode)
 INTEGER (kind=4) isec,lnode,inode,i,j,k,osec,err,ielem
 REAL (kind=8) :: diagm(nnode),voli,dens
 REAL (kind=8) :: young,fac,deltm,xlong,poiss,c1,c2,c3,c4,c5,c6,   &
                  xl(nvert),ls(3,nvert),t(3,3),ll,a1,a2,a3
 REAL (kind=8), PARAMETER :: a=0.333333333333333d0
 LOGICAL :: isotro
 TYPE (ele16), POINTER :: elem
 TYPE(section), POINTER :: sec

 ! -----------------------------
 IF( sms )thl = thl**2      !it works with the lengths squared
 fac = (1d0/dtscal/btscal)**2

 !       gauss points shape and derivatives of nodal functions

 CALL gaussq (ngaud ,posgp ,weigp )
 DO i=1,ngaud
   CALL shape4(nnode,shape(1,i),deriv(1,1,i),a,a,posgp(i))
 END DO
 wei = 0.5d0

 osec = -1
 elem => heade    !point to first element
 DO ielem=1,nelem
   x = coord(:,elem%lnods)
   ! ---generates lumped mass matrix proportional to diagonal
   diagm = 0d0
   DO k=1,ngaud
     CALL jacob5(cd(1,1),deriv(1,1,k),voli,x(1,1),nnode,999,err)
     voli = voli*wei(k)               !Gauss point mass
     diagm = diagm + shape(:,k)*voli
   END DO
   ! compute critical distance
   isec = elem%matno                    !associated material
   IF( isec /= osec )THEN               !present and previous material diff
     sec => psecs(isec)%p
     isotro = .FALSE.                   !default
     IF( sec%mabas >= 0) THEN
       isotro = sec%mtbas%mtype == 1 .OR. sec%mtbas%mtype == 6
       IF( isotro )THEN     !for an isotropic material
         young  = sec%mtbas%prope(1)        !young modulus
         poiss  = sec%mtbas%prope(2)        !poisson ratio
         !young  = young / (1d0 - 2d0*poiss) !modify
       ELSE                 !for an orthotropic material
         c1 = sec%mtbas%prope(45)/1.5d0
         c2 = sec%mtbas%prope(52)/1.5d0
         c3 = sec%mtbas%prope(59)/1.5d0
         c4 = 2d0*(sec%mtbas%prope(46)+2d0*sec%mtbas%prope(66)) /1.5d0
         c5 = 2d0*(sec%mtbas%prope(47)+2d0*sec%mtbas%prope(73)) /1.5d0
         c6 = 2d0*(sec%mtbas%prope(53)+2d0*sec%mtbas%prope(80)) /1.5d0
       END IF
     ELSE  !Multi-layered solid
       c1 = sec%rprop( 6)
       c2 = sec%rprop(13)
       c3 = sec%rprop(20)
       c4 = 2d0*(sec%rprop( 7)+2d0*sec%rprop(27))
       c5 = 2d0*(sec%rprop( 8)+2d0*sec%rprop(34))
       c6 = 2d0*(sec%rprop(14)+2d0*sec%rprop(41))
     END IF
     osec  = isec                       !update previous material
   END IF

   DO j = 1,nvert                         !for each side and diagonal
     ls(:,j) = x(:,h(1,j)) - x(:,h(2,j))      ! direction x length
     xl(j) = DOT_PRODUCT(ls(:,j),ls(:,j))     ! length (squared)
     IF( sms )THEN
       IF( xl(j) < thl ) xl(j) = xl(j)/alp
     END IF
   END DO

   IF( isotro )THEN     !isotropic cmaterial
     deltm = MINVAL(xl(1:nvert), MASK= xl(1:nvert) > 0d0)/young    !min length (squared)
   ELSE
     t = MATMUL(x,elem%cartd(:,:,1))        !triad at first Gauss point
     !DO j=1,3
     !  CALL vecuni(3,t(:,i),ll)             !unit vectors
     !END DO
     deltm = 1d30
     DO j = 1,nvert                         !for each side and diagonal
       IF( xl(j) == 0d0 )CYCLE
       a1 = DOT_PRODUCT(t(:,1),ls(:,j))  !proyections
       a2 = DOT_PRODUCT(t(:,2),ls(:,j))
       a3 = DOT_PRODUCT(t(:,3),ls(:,j))
       ll = xl(j)**3                     !length^6
       xlong = ll/ (c1*a1**4        + c2*a2**4         + c3*a3**4            &
                 +  c4*a1**2 *a2**2 + c5* a1**2 *a3**2 + c6*a2**2 *a3**2)
       IF(xlong <= 0d0 ) CYCLE           !check
       IF(xlong < deltm) deltm = xlong   !compare with previous min
     END DO
   END IF
   dens = fac/deltm
   ! ---add to global mass
   DO inode=1,nnode                        !for each element node
     lnode = elem%lnods(inode)             !node number
     ymass(1:3,lnode) = ymass(1:3,lnode) + diagm(inode)*dens
   END DO
   elem => elem%next                     !point to next element
 END DO
 RETURN
 END SUBROUTINE slum16

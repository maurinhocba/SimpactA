 SUBROUTINE delt04(nelem ,heade ,deltc , coord , btscal)

 !  Compute Critical time increment for a set of 3-D solid (TLF) elements

 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem     !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:), & !present mesh coordinates
                              btscal        !critical time factor
 REAL (kind=8), INTENT(OUT) :: deltc        !critical time increment
 TYPE (ele04), POINTER :: heade             !pointer to first element

 INTEGER, PARAMETER  :: nvert = 6
 INTEGER, PARAMETER  :: h(2,nvert) = (/ 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 /)
 INTEGER :: ielem,j,isec,osec
 REAL (kind=8) :: young,densi,fac,deltm,xlong,poiss,c1,c2,c3,c4,c5,c6,   &
                  xl(nvert),x(3,4),ls(3,nvert),t(3,3),ll,a1,a2,a3
 TYPE (ele04), POINTER :: elem
 TYPE (section), POINTER :: sec
 LOGICAL :: isotro

 deltm = (deltc/btscal)**2

 osec = -1                              !previous material number
 elem => heade                          !point to first element
 DO ielem=1,nelem                       !for each element in the set
   isec = elem%matno                    !associated material
   IF( isec /= osec )THEN               !present and previous material diff
     sec => psecs(isec)%p
     isotro = .FALSE.                   !default
     IF( sec%mabas >= 0) THEN
       isotro = sec%mtbas%mtype == 1 .OR. sec%mtbas%mtype == 6
       IF( isotro )THEN     !for an isotropic material
         densi  = sec%mtbas%prope(5)*mscal  !mass density
         young  = sec%mtbas%prope(1)        !young modulus
         poiss  = sec%mtbas%prope(2)        !poisson ratio
         young  = young / (1d0 - 2d0*poiss) !modify
         fac   = densi/young                !inverse of sound velocity (squared)
       ELSE                 !for an orthotropic material
         densi  = 1.5d0*sec%mtbas%prope(5)*mscal  !mass density
         c1 = sec%mtbas%prope(45)/densi
         c2 = sec%mtbas%prope(52)/densi
         c3 = sec%mtbas%prope(59)/densi
         c4 = 2d0*(sec%mtbas%prope(46)+2d0*sec%mtbas%prope(66)) /densi
         c5 = 2d0*(sec%mtbas%prope(47)+2d0*sec%mtbas%prope(73)) /densi
         c6 = 2d0*(sec%mtbas%prope(53)+2d0*sec%mtbas%prope(80)) /densi
       END IF
     ELSE  !Multi-layered solid
       densi = sec%rprop(3)*mscal  !mass density
       c1 = sec%rprop( 6)/densi
       c2 = sec%rprop(13)/densi
       c3 = sec%rprop(20)/densi
       c4 = 2d0*(sec%rprop( 7)+2d0*sec%rprop(27)) /densi
       c5 = 2d0*(sec%rprop( 8)+2d0*sec%rprop(34)) /densi
       c6 = 2d0*(sec%rprop(14)+2d0*sec%rprop(41)) /densi
     END IF
     osec  = isec                       !update previous material
   END IF
   x = coord(:,elem%lnods(1:4))                !nodal coordinates
   DO j = 1,nvert                         !for each side and diagonal
     ls(:,j) = x(:,h(1,j)) - x(:,h(2,j))      ! direction x length
     xl(j) = DOT_PRODUCT(ls(:,j),ls(:,j))     ! length (squared)
   END DO

   IF( isotro )THEN     !isotropic cmaterial
     xlong = MINVAL(xl(1:nvert), MASK= xl(1:nvert) > 0d0)*fac    !min length (squared)
     IF(xlong < deltm) deltm = xlong                     !compare with previous min
   ELSE
     t = MATMUL(x,elem%cd)        !triad at first Gauss point
     !DO j=1,3
     !  CALL vecuni(3,t(:,i),ll)             !unit vectors
     !END DO
     DO j = 1,nvert                         !for each side and diagonal
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
   elem => elem%next                         !point to next element
 END DO

 deltm = SQRT(deltm)*btscal
 IF( deltm < deltc) deltc = deltm

 RETURN
 END SUBROUTINE delt04

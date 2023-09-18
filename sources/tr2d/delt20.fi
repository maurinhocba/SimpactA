 SUBROUTINE delt20(nelem ,heade ,deltc , coord, sms, alp)

 !  Compute Critical time increment for a set of 2-D solid (TLF) elements

 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:) !present mesh coordinates
 REAL (kind=8), INTENT(IN OUT) :: alp    !alpha and layer thickness for SMS
 REAL (kind=8), INTENT(OUT) :: deltc     !critical time increment
 TYPE (ele20), POINTER :: heade          !pointer to first element
 LOGICAL :: sms                             !selective mass scaling is active for this set

 INTEGER :: ielem,n,isec,osec
 REAL (kind=8) :: young,densi,fac,deltm,xl(nnode),xlong,poiss,d(2),x(2,3)
 REAL (kind=8) :: a(2),b(2),area2
 TYPE (ele20), POINTER :: e
 TYPE (section), POINTER :: sec
 INTEGER, PARAMETER  :: h(2,3) = (/ 1,2, 2,3, 3,1 /)


 osec = -1                             !previous material number
 e => heade                         !point to first element
 DO ielem=1,nelem                      !for each element in the set
   isec = e%matno                  !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     young = sec%mtbas%prope(1)       !young modulus
     densi = sec%mtbas%prope(5)*mscal !mass density
     poiss = sec%mtbas%prope(2)       !poisson ratio
     IF( poiss > 0.46 ) young  = young / (1d0 - 2d0*poiss)/12d0 !modify
     fac   = SQRT(densi/young)        !inverse of sound velocity
     osec  = isec                     !update previous material
   END IF
   x = coord(:,e%lnods(1:3))                !nodal coordinates
   DO n=1,nnode                    !for each side
     d(:) = x(:,h(1,n)) - x(:,h(2,n))      ! direction x length
     xl(n) = DOT_PRODUCT(d,d)              ! length squared
   END DO
   IF( sms )THEN
     xl(1) = xl(1)/alp
     xlong = MINVAL(xl(1:nnode))
     xlong = SQRT(xlong)
   ELSE
     xlong = MAX(xl(1),xl(2),xl(3))
     xlong = SQRT(xlong)
     a(1:2) = x(1:2,2)-x(1:2,1)
     b(1:2) = x(1:2,3)-x(1:2,1)
     area2 = ABS (a(1)*b(2) - a(2)*b(1))
     xlong = area2 /xlong
   END IF
   !-----------------
   deltm = xlong*fac
   IF(deltm < deltc) deltc = deltm
   e => e%next                         !point to next element
 END DO

 RETURN
 END SUBROUTINE delt20

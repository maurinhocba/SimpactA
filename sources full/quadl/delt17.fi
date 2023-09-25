 SUBROUTINE delt17(nelem ,heade ,deltc , coord, sms, alp, thl )

 !  Compute Critical time increment for a set of 2-D solid (TLF) elements

 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:) !present mesh coordinates
 REAL (kind=8), INTENT(IN OUT) :: alp, thl      !alpha and layer thickness for SMS
 REAL (kind=8), INTENT(OUT) :: deltc     !critical time increment
 TYPE (ele17), POINTER :: heade          !pointer to first element
 LOGICAL :: sms                             !selective mass scaling is active for this set

 INTEGER :: ielem,n,isec,osec
 REAL (kind=8) :: young,densi,fac,deltm,xl(nnode),xlong,poiss,d(2),x(2,nnode)
 TYPE (ele17), POINTER :: elem
 TYPE (section), POINTER :: sec
 INTEGER, PARAMETER  :: h(2,4) = (/ 1,2, 2,3, 3,4, 1,4 /)


 IF( sms )thl = thl**2      !it works with the lengths squared
 osec = -1                            !previous material number
 elem => heade                        !point to first element
 DO ielem=1,nelem                     !for each element in the set
   isec = elem%matno                  !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     young = sec%mtbas%prope(1)       !young modulus
     poiss  = sec%mtbas%prope(2)        !poisson ratio
     young  = young / (1d0 - 2d0*poiss) !modify
     densi = sec%mtbas%prope(5)*mscal       !mass density
     fac   = SQRT(densi/young)        !inverse of sound velocity
     osec  = isec                     !update previous material
   END IF
   x = coord(:,elem%lnods)                !nodal coordinates
   DO n=1,nnode                    !for each side
     d(:) = x(:,h(1,n)) - x(:,h(2,n))      ! direction x length
     xl(n) = DOT_PRODUCT(d,d)              ! length squared
   END DO
   IF( sms )THEN
     xl(1) = xl(1)/alp
     xl(3) = xl(3)/alp
   END IF
   xlong = MINVAL(xl(1:nnode), mask= xl(1:nnode) > 0d0)
   xlong = SQRT(xlong)
   deltm = xlong*fac
   IF(deltm < deltc) deltc = deltm
   elem => elem%next                         !point to next element
 END DO

 RETURN
 END SUBROUTINE delt17

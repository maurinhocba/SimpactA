 SUBROUTINE delt19 (nelem ,heade ,deltc , coord, sms, alp, thl )

 !  Compute Critical time increment for a set of 2-D solid (TLF) elements

 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE
 INTEGER (kind=4), INTENT (IN) :: nelem  !number of elements in the set
 REAL (kind=8), INTENT(IN) :: coord(:,:) !present mesh coordinates
 REAL (kind=8), INTENT(IN OUT) :: alp, thl      !alpha and layer thickness for SMS
 REAL (kind=8), INTENT(OUT) :: deltc     !critical time increment
 TYPE (ele19), POINTER :: heade          !pointer to first element
 LOGICAL :: sms                             !selective mass scaling is active for this set

 INTEGER :: ielem,i,n,isec,osec,ncurv(2),numct
 REAL (kind=8) :: young,densi,fac,deltm,xl(nnode),xlong,poiss,d(2),x(2,3),val
 REAL (kind=8) :: a(2),b(2),area2
 TYPE (ele19), POINTER :: e
 TYPE (section), POINTER :: sec
 TYPE (curve), POINTER :: crvy,crvn !pointer to curves 
 INTEGER, PARAMETER  :: h(2,3) = (/ 1,2, 2,3, 3,1 /)

 osec = -1                             !previous material number
 e => heade                         !point to first element
 DO ielem=1,nelem                      !for each element in the set
   isec = e%matno                  !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     young = sec%mtbas%prope(1)       !young modulus
     poiss  = sec%mtbas%prope(2)      !poisson ratio
     densi = sec%mtbas%prope(5)*mscal !mass density     
     IF( poiss > 0.45 ) young  = young / (1d0 - 2d0*poiss) !modify
     fac   = SQRT(densi/young)        !inverse of sound velocity
     !associated curves for thermal dependant properties
     ncurv(1) = INT(sec%mtbas%props(1))        !associated curve for Young modulus
     ncurv(2) = INT(sec%mtbas%props(2))        !       "             Poisson ration
     IF( ncurv(1) > 0 ) CALL cur_point(ncurv(1),sec%mtbas%chead,crvy) !pointer to Young curve
     IF( ncurv(2) > 0 ) CALL cur_point(ncurv(2),sec%mtbas%chead,crvn) !pointer to Poisson curve
     osec  = isec                     !update previous material
   END IF
   numct = SUM(ncurv)
   IF( numct > 0 )THEN !compute properties as a function of temperature
     val = SUM(tempe(1,e%lnods(1:3)))/3d0 ! temperature at element center
     i = 1
     IF(ncurv(1) > 0) young = inte_cr(crvy%val,crvy%np,val,i)
     IF(ncurv(2) > 0) poiss = inte_cr(crvn%val,crvn%np,val,i)
     IF( poiss > 0.45 ) young  = young / (1d0 - 2d0*poiss) !modify
     fac = SQRT(densi/young)        !inverse of sound velocity
   END IF   
   
   x = coord(:,e%lnods(1:3))                !nodal coordinates
   DO n=1,nnode                    !for each side
     d(:) = x(:,h(1,n)) - x(:,h(2,n))      ! direction x length
     xl(n) = DOT_PRODUCT(d,d)              ! length squared
   END DO
   !xlong = MINVAL(xl(1:nnode), mask= xl(1:nnode) > 0d0)
   !xlong = SQRT(xlong)
   !-----------------
   xlong = MAX(xl(1),xl(2),xl(3))
   xlong = SQRT(xlong)
   a(1:2) = x(1:2,2)-x(1:2,1)
   b(1:2) = x(1:2,3)-x(1:2,1)
   area2 = ABS (a(1)*b(2) - a(2)*b(1))
   xlong = area2 /xlong
   !-----------------
   IF( sms .AND.  xlong < thl ) xlong = xlong/alp
   deltm = xlong*fac
   IF(deltm < deltc) deltc = deltm
   e => e%next                         !point to next element
 END DO

 RETURN
 END SUBROUTINE delt19

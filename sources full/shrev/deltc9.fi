 SUBROUTINE deltc9(nelem,nnode,deltc,heade,coord)
 !*********************************************************************
 !
 !     determines the maximum time increment for beam elements
 !
 !*********************************************************************
 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE
 !     arguments
 INTEGER (kind=4) :: nelem,nnode
 REAL (kind=8) :: deltc,coord(:,:)
 TYPE(ele09), POINTER :: heade

 !     local variables
 INTEGER (kind=4) ielem,isec,n1,n2,osec
 REAL (kind=8) long,young,ro,deltm,x(2),fac
 TYPE(ele09), POINTER :: e
 TYPE (section), POINTER :: sec
 ! REAL (kind=8),PARAMETER :: factor = 0.75D0

 osec = -1
 e => heade
 DO ielem=1,nelem
   isec = e%matno
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     young  = sec%mtbas%prope(1)      !young modulus
     ro     = sec%mtbas%prope(5)*mscal      !mass density
     fac    = SQRT(ro/young)          !inverse of sound vel.
     osec   = isec                    !update previous material
   END IF
   n1 = e%lnods(1)
   n2 = e%lnods(nnode)
   x  = coord(:,n1) - coord(:,n2)
   long = DOT_PRODUCT(x,x)
   IF(nnode == 3) long = long/2.5d0
   deltm = SQRT(long)*fac
   IF(deltm < deltc) deltc = deltm
   e => e%next
 END DO
 ! deltc = factor*deltc
 RETURN
 END SUBROUTINE deltc9

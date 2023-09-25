 SUBROUTINE codel6(nelem,deltc,heade,coord)
 !*********************************************************************
 !
 !     determines the maximum time increment
 !
 !*********************************************************************
 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: nelem
 TYPE(ele06), POINTER :: heade
 REAL (kind=8), INTENT(IN) :: coord(:,:)
 REAL (kind=8), INTENT(OUT) :: deltc

 INTEGER (kind=4) j,ielem,isec,n,osec
 TYPE(ele06), POINTER :: e
 TYPE (section), POINTER :: sec

 REAL    (kind=8) xlong,young,ro,xl(nnode),deltm,fac

 osec = -1
 e => heade
 DO ielem=1,nelem
   isec = e%matno                     !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     IF( sec%secty == 14 )THEN
       young  = sec%mtbas%prope(1)      !young modulus
       ro     = sec%mtbas%prope(5)*mscal      !mass density
     ELSE
       young  = sec%rprop(2)      !young modulus
       ro     = sec%rprop(3)*mscal      !mass density
     END IF
     fac    = SQRT(ro/young)          !inverse of sound vel.
     osec   = isec                    !update previous material
   END IF
   xl = 0d0
   DO j = 1,3 !ndime
     DO n=1,nnode-1
       xl(n) = xl(n) + (coord(j,e%lnods(n+1))-coord(j,e%lnods(n)))**2
     END DO
     xl(nnode) = xl(nnode)+(coord(j,e%lnods(1))-coord(j,e%lnods(nnode)))**2
   END DO
   xlong = MINVAL(xl, mask = xl > 0d0 )
   xlong = SQRT(xlong)
   deltm = xlong*fac
   IF(deltm < deltc) deltc = deltm
   e => e%next
 END DO

 RETURN
 END SUBROUTINE codel6

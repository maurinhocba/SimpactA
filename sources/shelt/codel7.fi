 SUBROUTINE codel7(nelem,deltc,heade,coord,ngaus)
 !*********************************************************************
 !
 !     determines the maximum time increment
 !
 !*********************************************************************
 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem,ngaus
 REAL (kind=8) :: deltc,coord(:,:)
 TYPE (ele07), POINTER :: heade

 INTEGER (kind=4) j,ielem,isec,ln(3),nvert,n,osec
 TYPE (ele07), POINTER :: e
 TYPE (section), POINTER :: sec

 REAL    (kind=8) xlong,young,ro,xl(3),deltm,fac

 osec = -1
 nvert = 3
 e => heade
 DO ielem=1,nelem
   isec = e%matno
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
   ln = e%lnods(1:nvert)
   DO j = 1,3 !ndime
     DO n=1,nvert-1
       xl(n) = xl(n) + (coord(j,ln(n+1))-coord(j,ln(n)))**2
     END DO
     xl(nvert) = xl(nvert)+(coord(j,ln(1))-coord(j,ln(nvert)))**2
   END DO
   xlong = SQRT(MINVAL(xl))  !minimum side length
   IF( ngaus /= 1)xlong = xlong / 2.5d0  !TQQL only
   deltm = xlong*fac
   IF(deltm < deltc) deltc = deltm
   e => e%next
 END DO

 RETURN
 END SUBROUTINE codel7

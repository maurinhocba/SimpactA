 SUBROUTINE deltc2(nelem,ndime,deltc,heade,coord)
 !*********************************************************************
 !
 !     determines the maximum time increment for beam elements
 !
 !*********************************************************************
 USE ctrl_db, ONLY : mscal
 IMPLICIT NONE

 !**   arguments
 INTEGER(kind=4), INTENT(IN) :: nelem,ndime
 REAL (kind=8), INTENT(IN) :: coord(:,:)
 REAL (kind=8), INTENT(IN OUT)   :: deltc
 TYPE (ele02), POINTER :: heade
 !**   local variables
 INTEGER (kind=4) j,ielem,lprop,n1,n2,omat
 REAL (kind=8)    long,young,ro,deltm
 TYPE (ele02), POINTER :: e
 TYPE (section), POINTER :: sec
 LOGICAL :: newm

 !      REAL (kind=8), PARAMETER :: factor=0.75d0

 !**   code
 omat = -1
 deltc = 1.d20
 e => heade
 DO ielem=1,nelem
   lprop = e%matno
   newm = lprop /= omat
   IF ( newm )THEN
     sec => psecs(lprop)%p
     omat = lprop
     young = sec%mtbas%prope(1)
     ro    = sec%mtbas%prope(5)*mscal
   END IF
   n1 = e%lnods(1)
   n2 = e%lnods(nnode)
   long = 0d0
   DO j = 1,ndime
     long = long + (coord(j,n1)-coord(j,n2))**2
   END DO
   deltm = SQRT(long*ro/young)
   IF(deltm < deltc) deltc = deltm
   e => e%next
 END DO
 !      deltc = factor*deltc
 RETURN
 END SUBROUTINE deltc2

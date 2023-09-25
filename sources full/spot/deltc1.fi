 SUBROUTINE deltc1(nnode,nelem,deltc,heade,emass)
 !*********************************************************************
 !
 !     determines the maximum time increment for spot elements
 !     present stiffness must be computed in terms of actual length
 !     something must be done about rotational stiffness also
 !
 !*********************************************************************
 IMPLICIT NONE

 !**   arguments
 INTEGER(kind=4), INTENT(IN) :: nnode,nelem
 REAL (kind=8), INTENT(IN) :: emass(:,:)
 REAL (kind=8), INTENT(IN OUT)   :: deltc
 TYPE (ele01), POINTER :: heade
 !**   local variables
 INTEGER (kind=4) j,ielem,isec,n1,n2,osec
 REAL (kind=8)    long,stiff,mass,deltm,tcrit !,penal
 TYPE (ele01), POINTER :: e
 TYPE (section), POINTER :: sec
 LOGICAL :: newm

 !      REAL (kind=8), PARAMETER :: factor=0.75d0

 !**   code
 osec = -1
 tcrit = 1.d20
 deltm = tcrit
 e => heade
 DO ielem=1,nelem
   isec = e%matno
   newm = isec /= osec
   IF ( newm )THEN
     sec => psecs(isec)%p
     osec = isec
     IF( nnode == 2)THEN
       stiff = sec%rprop(1)  !initial stiffness
     ELSE
       stiff = sec%rprop(7)*sec%mtbas%prope(1) ! Area * Young
       !penal= sec%rprop(7)/dtime**2
     END IF
   END IF
   IF( nnode == 2)THEN
     n1 = e%lnods(1)
     !x(:,1) = coord(:,n1)
     mass = emass(1,n1)
     n2 = e%lnods(2)
     IF( n2 > 0 ) THEN
       !x(:,2) = coord(:,n2)
       mass = MAX(mass,emass(1,n2))
     ELSE
       !x(:,2) = e%gausv(2:ndime+1)
     END IF
     !long = 0d0
     !DO j = 1,ndime
     !  long = long + (x(j,1)-x(j,2))**2
     !END DO
     deltm = SQRT(mass/stiff)
     IF(deltm < tcrit) tcrit = deltm
   ELSE
     long = e%gausv(11)
     DO j=1,3
       mass = emass(1,e%lnods(j))+ emass(1,e%lnods(j+3))
       deltm = SQRT(mass/stiff)*long
       IF(deltm < tcrit) tcrit = deltm
     END DO
   END IF
   e => e%next
 END DO
 IF( tcrit < deltc) deltc = tcrit
 RETURN
 END SUBROUTINE deltc1

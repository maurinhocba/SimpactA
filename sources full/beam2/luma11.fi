 SUBROUTINE luma11(nelem,ymass,sumat,heade,ntype,coord)
 !********************************************************************
 !
 !***  evaluates mass matrix for beam/shell elements
 !
 !********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem,ntype
 REAL (kind=8), INTENT(IN OUT) :: ymass(:,:),sumat,coord(:,:)
 TYPE( ele11 ), POINTER :: heade

 !        local variables
 INTEGER (kind=4), PARAMETER :: nvar=4
 INTEGER (kind=4) ielem,l,n,isec,osec !,matty
 REAL    (kind=8) dvolu,a(2),tmass,sume1,r0
 TYPE( ele11 ), POINTER :: e
 REAL (kind=8),PARAMETER :: twopi=6.283185307179586d0
 TYPE(section), POINTER :: sec

 osec = -1
 sume1 = 0d0
 e => heade
 DO ielem=1,nelem
   isec = e%matno
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     !matty = sec%secdef(3)
     IF(ntype == 1)THEN
       tmass = sec%mtbas%prope(5)*sec%rprop(4)  !dens*area
     ELSE
       tmass = sec%mtbas%prope(5)*sec%rprop(1)  !dens*thickness
     END IF
     osec  = isec
   END IF

   IF(ntype == 3) THEN              !shells of revolution
!     dvolu = e%l0(2)*twopi/6d0
     dvolu = e%l0(2)*twopi/4d0
     r0 = coord(1,e%lnods(2)) + coord(1,e%lnods(3))
!     a(1) = dvolu*(r0+coord(1,e%lnods(2)))
!     a(2) = dvolu*(r0+coord(1,e%lnods(3)))
     a = dvolu*r0
   ELSE                             !beams or 2-d Shells
     dvolu = e%l0(2)/6d0
     a(1)  = dvolu*3d0
     a(2)  = dvolu*3d0
   END IF

   DO n=1,2
     l = e%lnods(n+1)
     ymass(1:2,l) = ymass(1:2,l) + a(n)*tmass
   END DO

   sume1 = sume1 + SUM(a)*tmass

   e => e%next
 END DO
 sumat = sumat + sume1

 RETURN
 END SUBROUTINE luma11

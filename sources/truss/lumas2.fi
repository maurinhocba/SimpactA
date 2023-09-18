 SUBROUTINE lumas2(ndime,nelem,heade,emass,iwrit,sumat)
 !********************************************************************
 !
 !***   evaluation of mass matrix for truss elements
 !
 !********************************************************************
 IMPLICIT NONE

 !***  dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ndime,nelem,iwrit
 REAL (kind=8), INTENT(IN OUT) :: emass(:,:),sumat
 TYPE(ele02), POINTER :: heade

 !***  local variables
 INTEGER (kind=4) ielem,isec,l,n,osec
 REAL    (kind=8) tmass,sume1,area,dens
 TYPE(ele02), POINTER :: e
 TYPE(section), POINTER :: sec
 LOGICAL :: newm

 !***  loop over all the elements

 osec = -1
 e => heade
 sume1 = 0d0
 DO ielem=1,nelem
   isec  = e%matno
   newm = isec /= osec
   IF ( newm )THEN
     sec => psecs(isec)%p
     area = sec%rprop(1)
     dens = sec%mtbas%prope(5)
     osec = isec
   END IF
   tmass = e%gausv(1)*area*dens/2d0
   DO n=1,nnode
     l = e%lnods(n)
     emass(1:ndime,l) = emass(1:ndime,l) + tmass
   END DO
   sume1 = sume1 + 2d0*tmass
   e => e%next
 END DO
 IF(iwrit == 1) &
 WRITE(lures,"(//'  total mass for TRUSS element :', e15.7)",ERR=9999) sume1
 sumat = sumat + sume1

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lumas2

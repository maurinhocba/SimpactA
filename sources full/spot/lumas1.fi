 SUBROUTINE lumas1(ndime,ndofn,nelem,heade,emass,iwrit,sumat)
 !********************************************************************
 !
 !***   evaluation of mass matrix for truss elements
 !
 !********************************************************************
 IMPLICIT NONE

 !***  dummy arguments
 INTEGER (kind=4), INTENT(IN) :: ndime,ndofn,nelem,iwrit
 REAL (kind=8), INTENT(IN OUT) :: emass(:,:),sumat
 TYPE(ele01), POINTER :: heade

 !***  local variables
 INTEGER (kind=4) ielem,isec,l,n,osec
 REAL    (kind=8) sume1,tmass,rmass
 TYPE(ele01), POINTER :: e
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
     tmass = sec%rprop(5)/2d0
     rmass = sec%rprop(6)/2d0
     osec = isec
   END IF
   DO n=1,2
     l = e%lnods(n)
     IF( l > 0 )THEN
       emass(1:ndime,l) = emass(1:ndime,l) + tmass
       IF( e%rot )emass(ndime+1:ndofn,l) = emass(1+ndime:ndofn,l) + rmass
     END IF
   END DO
   sume1 = sume1 + 2d0*tmass
   e => e%next
 END DO
 IF(iwrit == 1) &
 WRITE(lures,"(//'  total mass for SPOT  element :', e15.7)",ERR=9999) sume1
 sumat = sumat + sume1

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lumas1

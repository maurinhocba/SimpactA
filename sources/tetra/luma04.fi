 SUBROUTINE luma04(heade,emass,sumat)

 IMPLICIT NONE

 INTEGER (kind=4), PARAMETER :: ndime=3
 REAL (kind=8), INTENT(IN OUT) :: emass(:,:),sumat
 TYPE( ele04 ), POINTER :: heade

 INTEGER (kind=4) :: isec,inode,lnode,idofn,osec
 REAL (kind=8) :: diag,suma,voli,rhoel,sume1
 TYPE( ele04 ), POINTER :: e
 TYPE (section), POINTER :: sec

 osec = -1

 sume1=0.
 e => heade
 DO
   IF( .NOT.ASSOCIATED (e) )EXIT
   suma=0.
   isec = e%matno
   IF( osec /= isec )THEN
     sec => psecs(isec)%p
     rhoel = sec%mtbas%prope(5)                !material density
     osec = isec
   END IF
   voli = e%dvol*rhoel
   suma = suma + voli
   diag = voli/4d0

   !---generates lumped mass matrix proportional to diagonal

   DO inode=1,nnode
     lnode = e%lnods(inode)
     DO idofn=1,ndime
       emass(idofn,lnode)=emass(idofn,lnode)+diag
     END DO
   END DO
   sume1 = sume1 + suma
   e => e%next
 END DO
 !IF(iwrit == 1) &
 WRITE(lures,"(//'  total mass for TETRA element :', e15.7)")sume1
 sumat = sumat + sume1
 RETURN
 END SUBROUTINE luma04

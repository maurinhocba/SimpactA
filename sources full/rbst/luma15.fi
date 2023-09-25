 SUBROUTINE luma15(heade,ymass,sumat,nrf,rhead)
 IMPLICIT NONE
 INTEGER (kind=4), PARAMETER :: nnr = 2, ndime = 3,  nnode = 3
 INTEGER (kind=4), INTENT(IN) :: nrf
 REAL (kind=8), INTENT(IN OUT) :: ymass(:,:),sumat
 TYPE( ele15 ), POINTER :: heade
 TYPE( reinf ), POINTER :: rhead

 INTEGER (kind=4) :: inode,lnode,isec,osec
 REAL (kind=8) :: diag,suma,voli,rhoel,thick,sume1
 TYPE( ele15 ), POINTER :: e
 TYPE( reinf ), POINTER :: rf
 TYPE (section), POINTER :: sec

 sume1=0.
 osec = -1
 e => heade
 DO
   IF( .NOT.ASSOCIATED (e) )EXIT
   isec = e%matno                    !associated material
   IF( isec /= osec )THEN
     sec => psecs(isec)%p
     IF( sec%secty == 12 )THEN
       thick = sec%rprop(1)                !shell thickness
       rhoel = sec%mtbas%prope(5)*thick    !surface density
     ELSE
       rhoel = sec%rprop(5)                !surface density
     END IF
     osec = isec
   END IF
   voli = e%area1*rhoel
   diag = voli/3d0

   !---generates lumped mass matrix proportional to diagonal

   DO inode=1,nnode
     lnode = e%lnods(inode)
     ymass(1:ndime,lnode)=ymass(1:ndime,lnode)+diag
   END DO
   sume1 = sume1 + voli
   e => e%next
 END DO
 suma = 0d0
 IF( nrf > 0 )THEN
   osec = -1
   rf => rhead
   DO
     IF( .NOT.ASSOCIATED (rf) )EXIT
     isec = rf%secno                    !associated material
     IF( isec /= osec )THEN
       sec => psecs(isec)%p
       rhoel = sec%rprop(1)             !linear density
       osec = isec
     END IF
     voli = rf%l0*rhoel
     diag = voli/2d0

     !---generates lumped mass matrix proportional to diagonal

     DO inode=1,nnr
       lnode = rf%lnods(inode)
       ymass(1:ndime,lnode)=ymass(1:ndime,lnode)+diag
     END DO
     suma = suma + voli
     rf => rf%next
   END DO
 END IF
 sumat = sumat + sume1 + suma
 !IF(iwrit == 1) &
 WRITE(lures,"(//'  total mass for shell element :', e15.7)")sume1
 !IF(iwrit == 1) &
 WRITE(lures,"(//'  total mass for stiffners ... :', e15.7)")suma
 RETURN
 END SUBROUTINE luma15

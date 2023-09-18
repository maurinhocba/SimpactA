 SUBROUTINE luma25(heade,ymass,sumat)
 IMPLICIT NONE
 REAL (kind=8), INTENT(IN OUT) :: ymass(:,:),  & !nodal mass in each direction
                                  sumat          !total mass of elements
 TYPE( ele25 ), POINTER :: heade                 !pointer to first element

 INTEGER (kind=4), PARAMETER :: ndime=3,nnode=4
 INTEGER (kind=4) :: i,l,isec,osec
 REAL (kind=8) :: voli,rhoel,thick,fac
 TYPE( ele25 ), POINTER :: e
 TYPE (section), POINTER :: sec


 osec = -1
 e => heade      !point to first element
 DO
   IF( .NOT.ASSOCIATED (e) )EXIT       !if all elements processed, exit
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
   voli = e%area(0)*rhoel            !element mass

   !---generates lumped mass matrix proportional to diagonal

   DO i=1,nnode
     l = e%lnods(i)
     fac = (e%area(i)+e%area(nextn(i+2)))/2d0*voli  !nodal factor
     ymass(1:3,l) = ymass(1:3,l) + fac
   END DO
   sumat = sumat + voli
   e => e%next
 END DO
 RETURN
 END SUBROUTINE luma25

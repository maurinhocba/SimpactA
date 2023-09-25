 SUBROUTINE lumas6(ndime,heade,emass,coord,sumat,iwrit,zigzag)
 !********************************************************************
 !
 !***  evaluates mass matrix for shell/quad elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: ndime,iwrit
 REAL (kind=8), INTENT(IN) :: coord(:,:)
 LOGICAL, INTENT(IN) :: zigzag
 REAL (kind=8), INTENT(IN OUT) :: emass(:,:),sumat
 TYPE( ele06 ), POINTER :: heade

 !        local variables
 TYPE( ele06 ), POINTER :: e
 TYPE(section), POINTER :: sec
 REAL (kind=8) :: thick
 INTEGER (kind=4) g,l,n,isec,osec
 REAL (kind=8) sume1,tmass,rhoh,facto,auxvl,x(3,nnode), &
               long,xl(nnode),trmas(nnode),ix(0:2),im(2)

 osec = -1
 sume1 = 0.0
 e =>heade
 DO
   IF( .NOT.ASSOCIATED(e) )EXIT
   x = coord(:,e%lnods)
   !       evaluates minimun mesh size (long)
   xl = 0.0
   DO g = 1,ndime
     DO n=1,nnode-1
       xl(n) = xl(n) + (x(g,n+1)-x(g,n))**2
     END DO
     xl(nnode) = xl(nnode) + (x(g,1)-x(g,nnode))**2
   END DO
   long = MAXVAL(xl)
   long = (SQRT(long)*nnode)/nnode
   !          lumped mass for each element
   tmass = 0.0
   isec  = e%matno                    !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     thick  = sec%rprop(1)
     rhoh   = sec%rprop(5)
     IF( sec%secty == 15 .AND. zigzag ) ix = sec%rprop(74:76)
     osec   = isec                    !update previous material
   END IF

   IF(long/thick <= SQRT(2d0/3d0)) THEN
     IF (zigzag )THEN
       facto = ix(0)/rhoh
     ELSE
       facto = thick**2/12d0    !standard
     END IF
   ELSE
     facto = long**2/8d0      !Hughes correction
   END IF

   trmas = 0d0
   DO g=1,ngaus
     auxvl = e%dvolu(g)*rhoh
     tmass = tmass+auxvl
     DO n=1,nnode
       trmas(n) = trmas(n) + auxvl*shape(n,g)
     END DO
   END DO
   !       sums over the global vector
   sume1 = sume1 + tmass
   DO n=1,nnode
     l = e%lnods(n)
     emass(1:3,l) = emass(1:3,l) + trmas(n)
     emass(4:5,l) = emass(4:5,l)+trmas(n)*facto
     IF( zigzag )THEN
       im(1) = ix(1)*e%nangl(1,n)**2 + ix(2)*e%nangl(2,n)**2/rhoh
       im(2) = ix(2)*e%nangl(1,n)**2 + ix(1)*e%nangl(2,n)**2/rhoh
       emass(7,l) = emass(7,l)+trmas(n)*im(1)
       emass(8,l) = emass(8,l)+trmas(n)*im(2)
     END IF
   END DO
   e => e%next
 END DO

 IF(iwrit == 1) &
 WRITE(lures,"(//'  total mass for shell element :', e15.7)",ERR=9999) sume1
 sumat = sumat + sume1

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lumas6

 SUBROUTINE lumas7(heade,coord,emass,sumat,iwrit,ngaus,zigzag)
 !********************************************************************
 !
 !***  evaluates mass matrix for shell/triangular elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: iwrit,ngaus
 REAL(kind=8), INTENT(IN) :: coord(:,:)
 REAL (kind=8), INTENT(IN OUT) :: sumat,emass(:,:)
 TYPE( ele07 ), POINTER :: heade
 LOGICAL, INTENT(IN) :: zigzag
!          local variables
 TYPE( ele07 ), POINTER :: e
 TYPE(section), POINTER :: sec
 REAL (kind=8) :: thick
 INTEGER (kind=4) j,l,n,isec,nvert,osec
 REAL    (kind=8) sume1,tmass,rhoh,facto,auxvl,x(3,6),x1(3),x2(3), &
                 x3(3),long,xl(4),trmas(6),cval,ix(0:2),im(2)
 INTEGER(kind=4), SAVE :: i(4) = (/6,4,5,6 /), k(4) = (/ 4,5,6,5 /)

 sume1 = 0.0

 osec = -1
 nvert = 3
 e => heade
 DO
   IF( .NOT.ASSOCIATED (e) )EXIT
   isec = e%matno                  !associated material
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     thick  = sec%rprop(1)
     cval = thick**2
     rhoh   = sec%rprop(5)            !surface density
     IF( sec%secty == 15 .AND. zigzag ) THEN
       ix = sec%rprop(74:76)/rhoh
       IF( ix(1) < 0.25d0 .OR. ix(2) < 0.25d0 )ix(1:2) = 0.25d0
       cval = ix(0)*12d0
     END IF
     osec   = isec                    !update previous material
   END IF

   x = coord(:,e%lnods)            !nodal coordinates

   !       evaluates minimun mesh size (long)
   xl(1:nvert) = 0.0
   DO j = 1,3 !ndime
     DO n=1,nvert-1
       xl(n) = xl(n) + (x(j,n+1)-x(j,n))**2
     END DO
     xl(nvert) = xl(nvert) + (x(j,1)-x(j,nvert))**2
   END DO
   IF( ngaus == 1 )THEN !TLLL
     long = MAXVAL(xl(1:nvert))           !distance between vertex nodes
   ELSE                 !TQQL
     long = MAXVAL(xl(1:nvert))/4d0       !distance between side nodes
   END IF
   !          lumped mass for each element
   IF(long <= cval ) THEN ! for small elements
     IF (zigzag )THEN
       facto = ix(0)*4d0
     ELSE
       facto = cval/3d0    !standard
     END IF
   ELSE
     !IF (zigzag )THEN
     !  facto = ix(0)
     !ELSE
       facto = long/3d0      !Hughes correction (modified)
     !END IF
   END IF

   tmass = 0d0
   trmas = 0d0
   IF( ngaus /= 1 )THEN
     !         for each sub-triangle
     DO n=1,4
       x1 = x(:,n) - x(:,k(n))
       x2 = x(:,n) - x(:,i(n))
       CALL vecpro(x1,x2,x3)
       CALL vecuni(3,x3,auxvl)
       auxvl = auxvl*rhoh/6d0
       trmas(n)    = trmas(n)   +auxvl
       trmas(i(n)) = trmas(i(n))+auxvl
       trmas(k(n)) = trmas(k(n))+auxvl
       tmass = tmass+auxvl*3d0
     END DO
   ELSE
     auxvl = e%dvolu(1)*rhoh/3d0
     trmas(1:3) = auxvl
     tmass = 3d0*auxvl
   END IF
   sume1 = sume1 + tmass
   IF( ngaus /= 1)THEN
     DO n=1,6
       l = e%lnods(n)
       emass(1:3,l) = emass(1:3,l)+ trmas(n)
       emass(4:5,l) = emass(4:5,l)+ trmas(n)*facto
       IF( zigzag )THEN
         im(1) = ix(1)*e%nangl(1,n)**2 + ix(2)*e%nangl(2,n)**2
         im(2) = ix(2)*e%nangl(1,n)**2 + ix(1)*e%nangl(2,n)**2
         emass(7,l) = emass(7,l)+trmas(n)*im(1)
         emass(8,l) = emass(8,l)+trmas(n)*im(2)
       END IF
     END DO
   ELSE               !TLLL
     DO n=1,nvert
       l = e%lnods(n)
       emass(1:3,l) = emass(1:3,l) + trmas(n)
       IF( zigzag )THEN
         im(1) = ix(1)*e%nangl(1,n)**2 + ix(2)*e%nangl(2,n)**2
         im(2) = ix(2)*e%nangl(1,n)**2 + ix(1)*e%nangl(2,n)**2
         emass(7,l) = emass(7,l)+trmas(n)*im(1)
         emass(8,l) = emass(8,l)+trmas(n)*im(2)
       END IF
       l = e%lnods(n+3)
       emass(4:5,l) = emass(4:5,l)+trmas(n)*facto
     END DO
   END IF
   e => e%next
 END DO

 IF( iwrit == 1 ) &
 WRITE(lures,"(//'  total mass for shell element :', e15.7)",ERR=9999) sume1
 sumat = sumat + sume1

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lumas7

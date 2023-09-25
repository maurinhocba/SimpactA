 SUBROUTINE masm17(heade,emass,mass,sumat,lumpd,ngaus)

 ! compute lumped/consistent mass  for element 17 (2-D solid TLF)

 IMPLICIT NONE
 LOGICAL, INTENT(IN) :: lumpd
 INTEGER (kind=4), INTENT(IN) :: ngaus
 REAL (kind=8), INTENT(IN OUT) :: emass(:,:),sumat,mass(:)
 TYPE( ele17 ), POINTER :: heade

 INTEGER (kind=4) :: isec,idofn,osec,g,i,j,k,ii,n,l,ng,ngm
 REAL (kind=8) :: suma,voli,rhoel,ms(36),a(10),m(4),xita(2),dvo(4)
 TYPE( ele17 ), POINTER :: e
 TYPE (section), POINTER :: sec

 REAL (kind=8) :: sh(4,4),der(4,2)
 REAL(kind=8), PARAMETER :: r3 = 0.577350269189626D+00
 INTEGER (kind= 4) posmt
 posmt(i,j) = (16-i)*(i-1)/2+j  !position i,j in mass(nxn)

 ! Note: consistent mass matrix requires at least 2x2 integration
 ng = 2
 IF( ngaus == 1) ng = 1
 IF( ng == 2 .OR. .NOT.lumpd)THEN  !two points or consistent mass matrix
   xita = (/ -r3,r3 /)  !standard 2x2
   ng = 2               !
   ngm = 4              !number of gauss points for the matrix
 ELSE                   !one integration point and lumped matrix
   xita = 0d0
   ngm = 1
 END IF
 suma=0.
 osec = -1

 ms = 0D0
 g = 0
 DO i=1,ng
   DO j=1,ng
     g = g+1
     CALL shape3(der,sh(1,g),xita(j),xita(i),nnode)
   END DO
 END DO

 e => heade
 DO
   IF( .NOT.ASSOCIATED (e) )EXIT
   isec = e%matno
   IF( osec /= isec )THEN
     sec => psecs(isec)%p
     rhoel = sec%mtbas%prope(5)                !material density
     osec = isec
   END IF
   m = 0
   IF(.NOT.lumpd) a = 0d0
   IF( ngm == ngaus )THEN  !same number of points
     dvo(1:ngm) = e%dvol           !use computed values
   ELSE                    !different number of points
     dvo = e%dvol(1)/ngm   !correct volume values
   END IF
   DO g=1,ngm

     voli = dvo(g)*rhoel
     suma = suma + voli

     ii = 0                   !Initializes Position in matrix
     DO n = 1,nnode           !For each node
        m(n) = m(n) + sh(n,g)*voli  !associated (lumped) mass
        IF(.NOT.lumpd) THEN    !Consistent matrix
           DO l = n,nnode       !for each node
             ii = ii+1          !position in matrix
             a(ii) = a(ii) + sh(n,g)*sh(l,g)*voli
           END DO
        END IF     !.NOT.lumpd
     END DO     !n=1,nnode
   END DO     !g=1,ngaus

   !---generates lumped mass matrix proportional to diagonal

   DO n=1,nnode
     l = e%lnods(n)
     DO idofn=1,2
       emass(idofn,l)=emass(idofn,l)+m(n)
     END DO
   END DO

   IF(.NOT.lumpd) THEN
                 !consistent mass matrix  Assembles over MASS(maxa)
     i = 0                            ! Position in matrix A initialization
     g = 0                            ! elemental DOF initialization
     DO n = 1,nnode                   ! For each nodee
       j = (n-1)*2                ! Initializes Column position in M
       DO l = n,nnode                 ! For each node
         i = i+1                      ! Matrix position in A
         k = (l-1)*2              ! Initializes Row position in M
         DO ii = 1,2              ! For each direction in space
           ms(posmt(j+ii,k+ii)) = a(i)
         END DO
       END DO
     END DO
     CALL masaux(4,e%lnods,8,2,ms(1),mass(1) )
   END IF

   e => e%next
 END DO
 !IF(iwrit == 1) &
 !WRITE(lures,"(//'  total mass for QUADL element :', e15.7)",ERR=9999)suma
 sumat = sumat + suma
 RETURN
 9999 CALL runen2(' ')
 END SUBROUTINE masm17

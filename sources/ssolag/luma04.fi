 SUBROUTINE luma04(heade,ymass,sumat)

 IMPLICIT NONE
 REAL (kind=8), INTENT(IN OUT) :: ymass(:,:),sumat
 TYPE( ele04 ), POINTER :: heade

 INTEGER (kind=4), PARAMETER ::  ngaus=3
 INTEGER (kind=4) :: isec,osec,g,n,l
 REAL (kind=8) :: suma,voli,rhoel,m(nnode)
 TYPE( ele04 ), POINTER :: e
 TYPE (section), POINTER :: sec
 REAL (kind=8) :: shape(nnode,ngaus),der(nnode,3), &
                  posgp(ngaus),weigp(ngaus),sh(3,ngaus)

 !  gauss points positions and shape functions
 CALL gaussq (ngaus ,posgp ,weigp )
 weigp = weigp*4d0
 DO g=1,ngaus
   CALL shape5(der(1,1),shape(1,g),0d0,0d0,posgp(g),nnode)
   sh(1,g) = -0.5d0*(1d0-posgp(g))*posgp(g)
   sh(2,g) =  1d0-posgp(g)**2
   sh(3,g) =  0.5d0*(1d0+posgp(g))*posgp(g)
 END DO

 ! initialization
 suma=0.
 osec = -1
 e => heade
 DO
   IF( .NOT.ASSOCIATED (e) )EXIT
   !  material and density
   isec = e%matno
   IF( osec /= isec )THEN
     sec => psecs(isec)%p
     IF( sec%mabas == -1 )THEN
       rhoel = sec%rprop(3)
     ELSE
       rhoel = sec%mtbas%prope(5)  !material density
     END IF
     osec = isec
   END IF

   m = 0            !initializes lumped mass matrix
   DO g=1,ngaus        !integration loop

     voli = DOT_PRODUCT(e%dvol,sh(:,g))*weigp(g)*rhoel  !GP mass
     suma = suma + voli                                 !add to total mass

     DO n = 1,nnode           !For each node
       m(n) = m(n) + shape(n,g)*voli  !associated (lumped) mass
     END DO     !n=1,nnode
   END DO     !g=1,ngaus

   !---generates lumped mass matrix proportional to diagonal

   DO n=1,nnode
     l = e%lnods(n)
     ymass(1:3,l) = ymass(1:3,l) + m(n)
   END DO

   e => e%next
 END DO
 !IF(iwrit == 1) &
 WRITE(lures,"(//'  total mass for SSOLAG element:', e15.7)")suma
 sumat = sumat + suma
 RETURN
 END SUBROUTINE luma04

 SUBROUTINE lumas8 (ndime,nelem,iwrit,nnode,ngaus,axesc,         &
                    heade,weigh,emass,shape,sumat)
 !******************************************************************
 !
 ! *** calculates lumped mass for beam elements
 !
 !******************************************************************
 IMPLICIT NONE
 !          routine PARAMETER
 INTEGER (kind=4), PARAMETER :: ndofe=6
 INTEGER (kind=4) :: ndime,nelem,nnode,ngaus,axesc,iwrit
 REAL   (kind=8) :: emass(:,:),weigh(:),shape(:,:),sumat
 TYPE (ele08), POINTER :: heade
 !          local variables
 INTEGER (kind=4) ielem,i,j,k,g,n,isec,l,osec
 REAL (kind=8) sume1,tmass,auxvl,area,rotms(3,3),trmas(nnode),   &
               iner(3),inert(3,3),        &
               aux1,aux2,aux3,aux4
 TYPE(section), POINTER :: sec
 TYPE (ele08), POINTER :: e
 LOGICAL isotro

 sume1 = 0.0
 osec = -1
 e => heade
 DO ielem=1,nelem
   ! lumped mass for each element
   isec = e%matno
   ! evaluates jacobians
   ! section properties
   IF( isec /= osec )THEN             !present and previous material diff
     sec => psecs(isec)%p
     isotro = sec%iprop(1) == 0       !isotropic material
     IF( isotro )THEN
       area = sec%rprop(1)/sec%rprop(16)*sec%rprop(15)   !A*density
       aux1 = SQRT(sec%rprop(1)/2d0/sec%rprop(5))        ![A/I2/2]^1/2
       aux2 = sec%rprop(5)/sec%rprop(16)*sec%rprop(15)   !I2*density
       aux3 = SQRT(sec%rprop(1)/2d0/sec%rprop(6))        ![A/I3/2]^1/2
       aux4 = sec%rprop(6)/sec%rprop(16)*sec%rprop(15)   !I3*density
     ELSE
       area = sec%rprop(1)                               !A*density
       aux1 = SQRT(sec%rprop(1)/2d0/sec%rprop(3))        ![A/I2/2]^1/2
       aux2 = sec%rprop(3)                               !I2*density
       aux3 = SQRT(sec%rprop(1)/2d0/sec%rprop(4))        ![A/I3/2]^1/2
       aux4 = sec%rprop(4)                               !I3*density
     END IF
     osec = isec                    !update previous material
   END IF
   ! correct this values according to Hughes
   auxvl = e%jac(1)*aux1
   IF(auxvl < 1d0) THEN
     iner(2) = aux2                                         !standard
   ELSE
     iner(2) = area*e%jac(1)**2/2d0                         !Hughes
   END IF
   auxvl = e%jac(1)*aux3
   IF(auxvl < 1d0) THEN
     iner(3) = aux4                                         !standard
   ELSE
     iner(3) = area*e%jac(1)**2/2d0                         !Hughes
   END IF
   iner(1) = iner(2)+iner(3)
   IF( .NOT.isotro )THEN
     IF( iner(1) < sec%rprop(2) ) iner(1) = sec%rprop(2)
   END IF

   ! numerical integration of nodal shape functions
   tmass = 0d0                     !initializes length
   trmas(1:nnode) = 0d0            !lumped lengths
   DO g=1,ngaus                       !for each Gauss point
     auxvl = e%jac(g)*weigh(g)        !associated length
     tmass = tmass + auxvl            !addt to total length
     DO n=1,nnode
       trmas(n) = trmas(n) + auxvl*shape(n,g)   !lumped lengths
     END DO
   END DO
   !auxvl = 0d0
   !it doesn't seem to improve the convergence
   DO n=1,nnode
     rotms(1,n) = trmas(n)*iner(1)
     rotms(2,n) = trmas(n)*iner(2) * 4d0  !increase bending mass
     rotms(3,n) = trmas(n)*iner(3) * 4d0  !increase bending mass
   END DO
   !IF(nnode == 3) THEN
   !  rotms(2,2) = rotms(2,2) + auxvl
   !  rotms(3,2) = rotms(3,2) + auxvl
   !END IF
   !transform inertia properties to local system & performs lumping
   IF(axesc > 0) THEN
     DO n=1,nnode            !for each node
       l = MIN(axesc,n)        !position in local lambda llbd(:,:,l)
       DO i=1,3                   !for each direction
         inert(i,n) = 0d0            !initializes inertia
         DO j=1,3                        !double loop on both index k,j
           DO k=1,3
             inert(i,n) = inert(i,n) + e%llbd(i,k,l)*e%llbd(j,k,l)*rotms(k,n)
           END DO
         END DO
       END DO
     END DO
   ELSE
     CALL vecasi(3*nnode,rotms,inert)
   END IF
   ! sums over the global vector
   sume1 = sume1 + tmass*area
   DO n=1,nnode
     l = e%lnods(n)
     emass(1:ndime,l) = emass(1:ndime,l) + trmas(n)*area
     DO j=ndime+1,ndofe
       emass(j,l) = emass(j,l) + inert(j-ndime,n)
     END DO
   END DO
   e => e%next
 END DO

 IF(iwrit == 1) WRITE(lures,"(//,'  Total mass for Beam element :', &
                        &    e15.7)",ERR=9999) sume1
 sumat = sumat + sume1

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE lumas8

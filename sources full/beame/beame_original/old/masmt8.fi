 SUBROUTINE masmt8 (ndime,nelem,iwrit,nnode,ngaus,axesc,         &
                    heade,weigh,emass,shape,sumat,mass,lumpd)
 !******************************************************************
 !
 ! *** calculates lumped mass for beam elements
 !
 !******************************************************************
 IMPLICIT NONE
 !          routine PARAMETER
 INTEGER (kind=4), PARAMETER :: ndofe=6
 LOGICAL, INTENT(IN) :: lumpd
 INTEGER(kind=4), INTENT(IN) :: ndime,nelem,nnode,ngaus,axesc,iwrit
 REAL   (kind=8), INTENT(IN) :: weigh(:),shape(:,:)
 REAL   (kind=8), INTENT(IN OUT) :: emass(:,:),sumat,mass(:)
 TYPE (ele08), POINTER :: heade
 !          local variables
 INTEGER (kind=4) ielem,i,j,k,g,n,m,isec,l,osec,nv,ii,jj,kk,li
 REAL (kind=8) sume1,tmass,auxvl,area,rotms(3,3),trmas(nnode),iner(3),inert(3,3),   &
               aux1,aux2,aux3,aux4
 REAL (kind=8) lti(3,3),ltil(3,3),mm(nnode,nnode),mc(171),pg(nnode),w(nnode),jc(nnode),sh(nnode,nnode)
 TYPE(section), POINTER :: sec
 TYPE (ele08), POINTER :: e
 LOGICAL isotro
 INTEGER (kind= 4) posmt
 posmt(i,j,nv) = (2*nv-i)*(i-1)/2+j  !position i,j in mass(nxn)

 INTERFACE
   INCLUDE 'shape9.h'
 END INTERFACE

 IF( .NOT.lumpd )  CALL shape9(w,pg,sh,mm,nnode,nnode)

 sume1 = 0.0
 osec = -1
 nv = nnode*ndofe  !number of total dof in the element
 mc = 0d0           !initializes just once
 e => heade
 DO ielem=1,nelem
   isec = e%matno                     !section number
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

   ! correct this values according to Hughes  ***************
   auxvl = e%jac(1)*aux1
   IF(auxvl < 1d0) THEN
     iner(2) = aux2                                         !standard
   ELSE
     iner(2) = area*e%jac(1)**2 /2d0                         !Hughes
   END IF
   auxvl = e%jac(1)*aux3
   IF(auxvl < 1d0) THEN
     iner(3) = aux4                                         !standard
   ELSE
     iner(3) = area*e%jac(1)**2 /2d0                         !Hughes
   END IF
   iner(1) = iner(2)+iner(3)
   IF( .NOT.isotro )THEN
     IF( iner(1) < sec%rprop(2) ) iner(1) = sec%rprop(2)
   END IF
   ! *******************

   ! numerical integration of nodal shape functions
   tmass = 0d0                     !initializes length
   trmas(1:nnode) = 0d0            !lumped lengths
   mm = 0d0
   DO g=1,ngaus                       !for each Gauss point
     auxvl = e%jac(g)*weigh(g)        !associated length
     tmass = tmass + auxvl            !addt to total length
     DO n=1,nnode
       trmas(n) = trmas(n) + auxvl*shape(n,g)   !lumped lengths
     END DO
   END DO
   DO n=1,nnode
     rotms(1,n) = trmas(n)*iner(1)
     rotms(2,n) = trmas(n)*iner(2)
     rotms(3,n) = trmas(n)*iner(3)
   END DO
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


   IF( .NOT.lumpd )THEN                !consistent mass matrix  Assembles over MASS(maxa)

     jc = SUM(e%jac)/ngaus
     mm = 0d0
     DO g=1,nnode                        !for each Gauss point
       DO n=1,nnode                        !for each Gauss point
         auxvl = jc(g)*w(g)                !associated length
         DO m=n,nnode
           mm(n,m) = mm(n,m) + auxvl*sh(n,g)*sh(m,g)
         END DO
       END DO
     END DO
     IF( axesc == 0 )THEN
       ltil = 0d0
       ltil(1,1) = iner(1)
       ltil(2,2) = iner(2)
       ltil(3,3) = iner(3)
     END IF

     DO n = 1,nnode                   ! For each nodee
       j = (n-1)*ndofe                ! Initializes Column position in M
       jj= j+3
       IF(axesc > 0) THEN
         l = MIN(axesc,n)               !position in local lambda llbd(:,:,l)
         DO ii = 1,3                  ! For each direction in space
           !lti(ii,:) = e%llbd(:,ii,l)*iner(ii)
           lti(:,ii) = e%llbd(:,ii,l)*iner(ii)
         END DO
       END IF
       DO m = n,nnode                 ! For each node
         k = (m-1)*ndofe              ! Initializes Row position in M
         kk= k+3                      ! Initializes Row position in M
         DO ii = 1,3                  ! For each direction in space
           mc(posmt(j+ii,k+ii,nv)) = mm(n,m)*area   !traslational
           !PRINT *,posmt(j+ii,k+ii,nv)
         END DO
         IF(axesc > 0) THEN
           l = MIN(axesc,m)             !position in local lambda llbd(:,:,l)
           !ltil = MATMUL(lti,e%llbd(:,:,l))
           ltil = MATMUL(lti,TRANSPOSE(e%llbd(:,:,l)))
         END IF
         DO i=1,3
           li = 1
           IF( n == m ) li = i
           DO ii = li,3                  ! For each direction in space
             mc(posmt(jj+i,kk+ii,nv)) = mm(n,m)*ltil(i,ii)  !rotational
             !PRINT *,posmt(jj+i,kk+ii,nv)
           END DO
         END DO
       END DO
     END DO
     CALL masaux(nnode,e%lnods,nv,ndofe,mc(1),mass(1))
   END IF
   e => e%next
 END DO

 IF(iwrit == 1) WRITE(lures,"(//,'  Total mass for Beam element :', &
                        &    e15.7)",ERR=9999) sume1
 sumat = sumat + sume1

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE masmt8

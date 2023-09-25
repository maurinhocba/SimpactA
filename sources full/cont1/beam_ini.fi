 SUBROUTINE beam_ini(nev,nec,lv,lc,iwrit)
 ! this routine computes initial displacements
 ! to match the axes of tube and skin
 IMPLICIT NONE
 ! dummy arguments
 INTEGER(KIND=4), INTENT(IN) :: nev,  & !number of beam elements
                                nec,  & !number of skin elements
                                iwrit   !flag to print
 INTEGER(KIND=4), POINTER ::  lv(:,:),  &  !beam connectivities
                              lc(:,:)      !skin connectivities
 ! parameters
 REAL(KIND=8) :: toler = 1e-8 , &
                 facto = 57.295779513082320876798154814105      ! facto= 180d0/pi
 ! local variables
 INTEGER(KIND=4) :: i,j,k,l,m,npv,npc
 LOGICAL :: rcv,rcc
 REAL(KIND=8) :: dx(3),l0,li,lj,ti(3),tj(3),ca,sa2,a

 npv = nev + 1       !number of beam points
 npc = nec + 1       !number of skin points

 !find first node (coincident nodes)
 rcv = .FALSE.        !modify beam connectivities
 rcc = .FALSE.        !modify skin connectivities
 IF( dist(coora(:,lv(1,1)),coora(:,lc(1,1))) < toler )THEN
   !right
 ELSE IF( dist(coora(:,lv(1,1)),coora(:,lc(2,nec))) < toler )THEN
   ! reverse skin connectivities
   rcv = .TRUE.
 ELSE IF( dist(coora(:,lv(2,nev)),coora(:,lc(1,1))) < toler )THEN
   ! reverse tube connectivities
   rcc = .TRUE.
 ELSE IF( dist(coora(:,lv(2,nev)),coora(:,lc(2,nec))) < toler )THEN
   ! reverse both connectivities
   rcv = .TRUE.
   rcc = .TRUE.
 END IF

 IF( rcv )THEN ! reverse tube  connectivities
   DO i=1,nev/2
     j = nev+1-i
     k = lv(1,i)
     l = lv(2,i)
     lv(1,i) = lv(2,j)
     lv(2,i) = lv(1,j)
     lv(1,j) = l
     lv(2,j) = k
   END DO
   IF( MOD(nev,2) == 1 )THEN
     j = nev/2+1
     k = lv(1,j)
     lv(1,j) = lv(2,j)
     lv(2,j) = k
   END IF
 END IF
 IF( rcc )THEN ! reverse skin  connectivities
   DO i=1,nec/2
     j = nec+1-i
     k = lc(1,i)
     l = lc(2,i)
     lc(1,i) = lc(2,j)
     lc(2,i) = lc(1,j)
     lc(1,j) = l
     lc(2,j) = k
   END DO
   IF( MOD(nec,2) == 1 )THEN
     j = nec/2+1
     k = lc(1,j)
     lc(1,j) = lc(2,j)
     lc(2,j) = k
   END IF
 END IF
 ! COMPUTE PRESENT COORDINATES
 ! firt check initial displacement array
 IF( .NOT.initial_displacements )THEN
   initial_displacements = .TRUE.
   coorc = 0d0
 END IF
 ! compute local system at first node
 k = lv(1,1)                !first node of first segment
 dx = euler(1:3,k)          !euler angles at first node
 CALL inrotm(dx,euler(:,k)) !first local system
 ! go computing new positinos
 k = 1                      !initializes position in skin
 ti = coora(:,lc(2,k)) - coora(:,lc(1,k))  !skin tangent vector
 li = SQRT(DOT_PRODUCT(ti,ti))             !skin length
 ti = ti/li                                !unit vector
 l = lv(1,1)                               !first beam node
 coorc(:,l) = coora(:,l)                   !update coordinates
 DO i=1,nev                 !loop over each beam node (but first)
   ! compute initial lengths
   m = lv(2,i)                             !second node beam
   l0 = dist(coora(:,m), coora(:,l))       !beam length
   IF( l0 <= li )THEN     !beam node is inside the skin element
     coorc(:,m) = coorc(:,l) + l0*ti       !new node coordinates
     li = li - l0                          !remaining skin length
   ELSE  !use quadratic equation to determine the position
     k = k+1                               !new skin segment
     tj = coora(:,lc(2,k)) - coora(:,lc(1,k))  !skin tangent vector
     lj = SQRT(DOT_PRODUCT(tj,tj))             !skin length
     tj = tj/lj                                !unit vector
     ca = DOT_PRODUCT(ti,tj)                   !cos between adjacent segments
     sa2 = 1d0 - ca**2                         !sin^2
     a = SQRT(l0**2+li**2*sa2) - li*ca         !distance along skin
     IF( a > lj )WRITE(*,"(' error SKIN discretization is finer than BEAM discretization')")
     coorc(:,m) = coorc(:,l) + li*ti + a*tj    !new node coordinates
     li = lj - a                               !remaining skin length
     ti = tj                                   !new unit vector
   END IF
   ! compute Euler angles
   dx = coorc(:,m) - coorc(:,l)  !present tangent vector (local x-axis) in beam
   a = SQRT(DOT_PRODUCT(dx,dx))  !compute length
   dx = dx/a                     !unit vector
   CALL vecpro(euler(1:3,l),dx,tj) !rotation vector  t1_i-1 x t1_i
   a = SQRT(DOT_PRODUCT(tj,tj))  !angle SIN
   IF( a > 0d0 )   tj = tj*ASIN(a)/a   !corrected rotation vector
   euler(:,m) = euler(:,l)             !previous local system
   dx(1) = DOT_PRODUCT(euler(1:3,m),tj)  !local rotation vector components
   dx(2) = DOT_PRODUCT(euler(4:6,m),tj)
   dx(3) = DOT_PRODUCT(euler(7:9,m),tj)
   CALL actrot(euler(:,m),dx)        !local system for new node
   l = m                             !update first beam node
 END DO
 ! store information
 IF(iwrit == 1) WRITE(lures,"('INITIAL_CONDITIONS',/,'DISPLACEMENT')")
 k = lv(1,1)           !first beam node
 coorc(:,k) = 0d0      !no displacements
 DO i=1,nev            !for each second node
   k = lv(2,i)         !second node
   coorc(:,k) = coorc(:,k)-coora(:,k) !displacements
   IF(iwrit == 1)WRITE(lures,"(i7,3e15.6)")label(k),coorc(:,k) !export to output file
 END DO
 IF(iwrit == 1)WRITE(lures,"('END_DISPLACEMENT',/,'LOCAL_SYSTEMS')")
 k = lv(1,1)                             !first beam node
 CALL angeul(euler(:,k), dx,.TRUE. )     !compute euler angles in rad
 euler(1:3,k) = dx                       !keep Euler angles
 euler(4:9,k) = 0d0                      !initializes as flag
 DO i=1,nev                              !for each second node
   k = lv(2,i)                           !second node
   CALL angeul(euler(:,k), dx, .TRUE. )  !compute euler angles in rad
   euler(1:3,k) = dx                     !keep Euler angles
   euler(4:9,k) = 0d0                    !initializes as flag
   IF(iwrit == 1)WRITE(lures,"(i7,3e15.6)")label(k),dx*facto  !print Euler angles
 END DO
 IF(iwrit == 1)WRITE(lures,"('END_LOCAL_SYSTEMS',/,'END_INITIAL_CONDITIONS')")
 RETURN
 END SUBROUTINE beam_ini

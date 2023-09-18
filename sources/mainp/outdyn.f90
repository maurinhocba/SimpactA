 SUBROUTINE outdyn(inicia, foutp, ecdis )
 !***********************************************************************
 !
 !** output routine
 !
 !***********************************************************************
 USE ctrl_db,ONLY: ifunc, ifixd, istep, numct, ndime, ndofn, nrotd, &
                   neulr, nload, npoin, dtime, ttime, xbase, &
                   begtm, neq , therm , tscal , itemp, ndoft
 USE kinc_db, ONLY : nvelr,nn,velor,npsdf,nesdf,ftsdf,nsdof !,sldvl
 USE lispa0
 USE npo_db, ONLY : label,coord,coora,ifpre,resid,fcont,force,euler,velnp,  &
                    veloc,acelr,ymass,naeul,ddisp,psia,tempe
 USE outp_db
 USE c_input, ONLY: del_old_files
 USE gvar_db, ONLY : static
 IMPLICIT NONE

   !--- Dummy variables
   LOGICAL,INTENT(IN):: foutp, & ! .TRUE. write compulsory
                        inicia   ! .TRUE. first call, write initial values
   INTEGER(kind=4),INTENT(IN):: ecdis  !internal equation number controlling output

   !--- Function
   REAL(kind=8):: functs  !time function
   !--- Local varibles
   INTEGER(kind=4):: ipoin,ireq,idofn,ieq,iq,ib,ie,i,j,k,nv1,npt
   REAL(kind=8):: velcr,cinet,poter,toler,time1,cputm,auxil,enkint,   &
                  enkinr,fract,tdisp,tstra,value,angls(3),rm(3,3),fac
   LOGICAL:: b1,b2
   ! auxiliar arrays to write down information
   REAL(kind=8),ALLOCATABLE :: disp1(:,:),velo1(:,:),acce1(:,:), disp3(:,:)
   REAL(kind=8),SAVE:: oldt1=0d0
   TYPE (slave_list), POINTER :: sl_d

   INTERFACE
     INCLUDE 'angeul.h'
     INCLUDE 'elemnt.h'
     INCLUDE 'contac.h'
   END INTERFACE


   IF( lastst .AND. .NOT.foutp )THEN
     lastst = .FALSE.
     RETURN
   END IF
   ! first determine if OUTPUT must be done at this step

   nv1 = nvelr+1            !last position in array of prescribed velocities
   CALL timuse(cputm)       !present CPU time
   cputm = cputm-cpui       !elapsed CPU time since the beginning of the process
   npt = INT(toutp(1))      !number of points defining output period
   ! initializes flags
   b1 = .FALSE.       !output at selected points
   b2 = .FALSE.       !global output

   IF (ecdis < 0) THEN
     !POSTYPE = 'C' curve control (curve number = -ECDIS)
     value = functs(-ecdis,ttime)                  !present associated function value
     auxil = functs(-ecdis,ttime+dtime)            !next step function value
     toler = (auxil-value)*(1.0000000001d0)/2d0    !tolerance
     !for Selected points (History)
     time1 = MODULO(value,toutd)                   !rest of the divition
     IF(time1+toler > toutd .OR. time1 < toler)  b1 = .TRUE.  !output at selected points
     !for global output (GiD or TecPlot)
     IF (npt == 1) THEN                !only one value
       auxil = toutp(2)                != PERIOD
       time1 = MODULO(value,auxil)     !rest of the divition
       IF (time1+toler > auxil .OR. time1 < toler)  b2 = .TRUE. !global output
     ELSE                       !multiple value definition
       DO i=1,npt                      !compare each value
         auxil = toutp(i+1)            !time
         time1 = value-auxil           !difference
         IF (ABS(time1) < toler)  b2 = .TRUE.  !Global OUTPUT selected
       END DO
     END IF

   ELSE IF(ecdis > 10)THEN
     !POSTYPE = 'D' displacemente control
     ipoin = ecdis/10             !node (internal)
     idofn = MOD(ecdis,10)        !DOFs
     ieq   = ifpre(idofn,ipoin)   !associated equation
     tdisp = ABS( coora(idofn,ipoin) - coord(idofn,ipoin) ) !total displacement
     IF(ieq > 0)THEN                                        !an active DOF
       IF( static )THEN
         toler = ABS(ddisp(ieq))*(1.2d0)/2d0           !tolerance
       ELSE
         toler = dtime*ABS(veloc(ieq))*(1.050000d0)/2d0     !tolerance
       END IF
     ELSE IF(ieq < -nn )THEN                                !a prescribed DOF
       toler = dtime*ABS(velor(-ieq-nn,nv1))*(1.010000d0)/2d0
     ELSE IF(ieq < 0 .AND. ieq > -nn)THEN                   !a slave DOF
       ib = npsdf(-ieq)                                   !first Master
       ie = npsdf(-ieq+1)-1                               !last Master
       auxil = 0d0                                        !initializes
       DO i = ib,ie                                       !loop to compute velocity
         iq = nesdf(i)                                  !associated equation of master DOF
         IF (iq > 0) THEN                               !if a Active DOF
           IF( static )THEN
             auxil = auxil + ddisp(iq)*ftsdf(i)
           ELSE
             auxil = auxil + veloc(iq)*ftsdf(i)
           END IF
         ELSE IF (iq < -nn) THEN                        !if a Prescribed DOF
           auxil = auxil + velor(-iq-nn,nv1)*ftsdf(i)
         END IF
       END DO
       toler = dtime*auxil*(1.050000d0)/2d0               !tolerance
     ELSE                                                 !ieq == 0 (no associated DOF)
       toler = 0d-0 !this should not happen
     END IF
     !for selected points (History)
     time1 = MODULO(tdisp,toutd)   !rest of the division
     IF (time1+toler > toutd .OR. time1 < toler )THEN
       IF( static ) THEN
         b1 = .TRUE.
       ELSE
         IF( ABS((tdisp-oldt1)/toutd) > 0.5d0) b1 = .TRUE.
       END IF
     END IF
     IF( b1 ) oldt1 = tdisp                               !keep last value
     !for Global Output (GiD or TecPlot)
     IF (npt == 1) THEN
       auxil = toutp(2)
       time1 = MODULO(tdisp,auxil)
       IF(time1+toler > auxil .OR. time1 < toler)  b2 = .TRUE.
     ELSE
       DO i=1,npt
         auxil = toutp(i+1)
         time1 = tdisp-auxil
         IF(ABS(time1) < toler)  b2 = .TRUE.
       END DO
     END IF
     value = tdisp
   ELSE
     !POSTYPE = 'T'  for time control
     toler = dtime*(1.0001d0)/2d0  !half the incremental time
     !for selected points (History)
     time1 = MODULO(ttime,toutd) !rest of the division
     IF (time1+toler > toutd .OR. time1 < toler )THEN !
       IF( static )THEN
         b1 = .TRUE.
       ELSE
         IF( ABS((ttime-oldt1)/toutd) > 0.5d0) b1 = .TRUE.
       END IF
     END IF
     IF( b1 ) oldt1 = ttime
     !for global output (GiD or TecPlot)
     IF(npt == 1)THEN
       auxil = toutp(2)
       time1 = MODULO(ttime,auxil)
       IF(time1+toler > auxil .OR. time1 < toler)  b2 = .TRUE.
     ELSE
       DO i=1,npt
         auxil = toutp(i+1)
         time1 = ttime-auxil
         IF(ABS(time1) < toler)  b2 = .TRUE.
       END DO
     END IF
     value = ttime
   END IF

   IF( foutp .OR. inicia) b1 = .TRUE.  !output compulsory of initial values
   IF( foutp ) b2 = .TRUE.             !output compulsory

   !***  OUTPUT FOR SELECTED VALUES

   IF(b1)THEN
     !         D I S P L A C E M E N T S
     IF(nreqd > 0) THEN
       ALLOCATE ( disp1(ndofn,nreqd) )
       DO ireq = 1,nreqd
         ipoin = nprqd(ireq)
         DO idofn = 1,ndime
           disp1(idofn,ireq) = (coora(idofn,ipoin)-coord(idofn,ipoin))
         END DO
         ! for local systems
         IF(neulr > 0) THEN
           IF(ndime == 2) THEN                   !2-D problems (1 DOF)
             disp1(ndofn,ireq) = euler(1,ipoin)   !present angle
           ELSE                                  !3-D problems (3 DOFs)
             rm = RESHAPE( euler(1:9,ipoin),(/3,3/))          !rotation matrix
             angls = 0d0                    !angles
             CALL angeul(rm,angls,.TRUE.)   !returns Euler angles (in rads)
             disp1(ndime+1:nrotd,ireq) =   angls(1:nrotd-ndime)
           END IF
         END IF
         IF( ndofn == 8 )disp1(7:8,ireq) = psia(:,ipoin)
       END DO
       WRITE(11,ERR=9999) value,((disp1(i,ireq),i=1,ndofn),ireq=1,nreqd)
       DEALLOCATE ( disp1 )
     END IF
     IF( .NOT.static )THEN
       !    V E L O C I T I E S
       IF(nreqv > 0) THEN
         ALLOCATE ( velo1(ndofn,nreqv) )
         DO ireq = 1,nreqv
           ipoin = nprqv(ireq)
           DO idofn = 1,ndofn
             ieq = ifpre(idofn,ipoin)
             IF(ieq > 0 )THEN                        !active DOFs
               velo1(idofn,ireq) = veloc(ieq)
             ELSE IF(ieq < -nn)THEN                  !fixed DOFs
               velo1(idofn,ireq) = velor(-ieq-nn,nv1)
             ELSE IF(ieq < 0 .AND. ieq > -nn)THEN    !slave DOFs
               ib = npsdf(-ieq)
               ie = npsdf(-ieq+1)-1
               auxil = 0d0
               DO i = ib,ie
                 iq = nesdf(i)
                 IF(iq > 0)THEN
                   auxil = auxil + veloc(iq)*ftsdf(i)
                 ELSE IF(iq < -nn)THEN
                   auxil = auxil + velor(-iq-nn,nv1)*ftsdf(i)
                 END IF
               END DO
               velo1(idofn,ireq) = auxil
             ELSE
               velo1(idofn,ireq) = 0d0               !non-existent DOFs
             END IF
           END DO
         END DO
         WRITE(19,ERR=9999) value,((velo1(i,ireq),i=1,ndofn),ireq=1,nreqv)
         DEALLOCATE ( velo1 )
       END IF
       !   A C C E L E R A T I O N S
       IF(nreqa > 0) THEN
         ALLOCATE ( acce1(ndofn,nreqa) )
         DO ireq = 1,nreqa
           ipoin = nprqa(ireq)
           DO idofn = 1,ndofn
             ieq = ifpre(idofn,ipoin)
             IF(ieq > 0 )THEN                        !active DOFs
               acce1(idofn,ireq) = acelr(ieq)
             ELSE
               acce1(idofn,ireq) = 0d0               !otherwise 0
             END IF
           END DO
         END DO
         WRITE(20,ERR=9999) value,((acce1(i,ireq),i=1,ndofn),ireq=1,nreqa)
         DEALLOCATE ( acce1 )
       END IF
     END IF
     !  N O D A L   E Q U I V A L E N T   F O R C E S
     IF(nreql > 0) THEN
       res = resid(:,nprql(1:nreql)) !direct contributions
       !write(58,"(3e15.5)")resid(1:3,365:377)
       !search for slave dependencies
       IF( nsdof > 0 )THEN
         sl_d => sl_head
         DO ipoin=1,nreql
           DO idofn=1,ndofn
             DO i=1,sl_d%nvalues
               j  = sl_d%deps(1,i)
               k  = sl_d%deps(2,i)
               ib = sl_d%deps(3,i)
               res(idofn,ipoin) = res(idofn,ipoin) + ftsdf(ib)*resid(j,k)
               !IF( idofn == 6 ) WRITE(58,"(3i5,3e15.5)")k,j,ib,ftsdf(ib),resid(j,k),ftsdf(ib)*resid(j,k)
             END DO
             !IF( idofn == 6 ) WRITE(58,"(3e15.5)")(ftsdf(sl_d%deps(3,i)),i=1,39)
             sl_d => sl_d%next
           END DO
         END DO

       END IF
       WRITE(12,ERR=9999) value,res(1:ndofn,1:nreql)
     END IF
     !  N O D A L   C O N T A C T   F O R C E S
     IF(numct > 0 .AND. nreqc > 0) WRITE(14,ERR=9999) value,fcont(1:ndime,nprqc(1:nreqc))
     !  P A I R S   C O N T A C T   F O R C E S
     IF(numct > 0) CALL contac('OUTDY1',0, ttime=value, dtcal = dtime)
     !  E N E R G Y   V A L U E S
     IF(iener /=  0 .AND. .NOT.static) THEN
       cinet = 0d0
       IF(nvelr > 0)THEN
         DO ieq = 1,neq
           cinet = cinet + veloc(ieq)**2/ymass(ieq)
         END DO
       END IF
       poter = 0d0
       enkint=0.
       DO ipoin=1,npoin
         DO idofn=1,ndime
           ieq = ifpre(idofn,ipoin)
           IF(ieq > 0)THEN
             tdisp = coora(idofn,ipoin) - coord(idofn,ipoin)
             IF(nload > 0)poter = poter +tdisp*force(ieq,nload+1)
             enkint= enkint+veloc(ieq)**2/ymass(ieq)
           END IF
         END DO
       END DO
       IF(neulr > 0)THEN
         enkinr = cinet - enkint
       ELSE
         enkinr = 0d0
         cinet  = enkint
       END IF
       velcr  = SQRT(enkint/sumat)
       enkint = 0.5*enkint
       enkinr = 0.5*enkinr
       cinet  = 0.5*cinet
       fract = 0.
       IF (energ(2)+cinet-energ(3).ne.0.) fract=1d0-energ(1)/(energ(2)+cinet-energ(3))
       WRITE(15,ERR=9999) value,cinet,velcr,dtime,enkint,enkinr,   &
         energ(1),energ(2),fract,cinet-energ(3),energ(5),poter
       !WRITE(55,"(' Time=',e12.4,' Kin. E=',e12.5,' Pot. E=',e12.5)",ERR=9999) ttime,cinet,poter
     END IF
     ! volume & pressure of volume dependent follower loads
     CALL wrtfl1 (ttime)
     !         T E M P E R A T U R E S
     IF(nreqT > 0) THEN
       ALLOCATE ( disp1(ndoft,nreqt) )
       DO ireq = 1,nreqt
         ipoin = nprqt(ireq)
         DO idofn = 1,ndoft
           disp1(idofn,ireq) = tempe(idofn,ipoin)
         END DO
       END DO
       WRITE(70,ERR=9999) value,((disp1(i,ireq),i=1,ndoft),ireq=1,nreqt)
       DEALLOCATE ( disp1 )
     END IF

   END IF

   !*** GLOBAL OUTPUT
   IF( .NOT. static) lastst = b2
   IF( b2 )THEN
     !message in .RSN file, if IWRIT == 1, results are printed to ASCII file
     WRITE(lures,"(5x,//'Results are reported for postprocess at:',/, &
       10x,'Istep= ',i8,5x,'Ttime= ',e15.7,/)",ERR=9999) istep,ttime
     WRITE(17,ERR=9999) istep,value,postype              !heading
     IF (ifunc==0)WRITE(17,ERR=9999) xbase(1:3)          !ground acc, vel & disp
     !** Total Displacements
     DO ipoin = 1,npoin
       WRITE(17,ERR=9999) (coora(1:ndime,ipoin) - coord(1:ndime,ipoin))
     END DO
     ! for local systems
     IF(neulr > 0) THEN
       fac = 45d0/ATAN(1d0)
       IF(ndime == 2) THEN                   !2-D problems (1 DOF)
         DO ipoin = 1,npoin
           WRITE(17,ERR=9999) euler(1,ipoin) !present angle
         END DO
       ELSE                                  !3-D problems (3 DOFs)
         DO ipoin = 1,npoin
           angls = 0d0                  !angles
           IF(naeul(ipoin) )THEN
             rm = RESHAPE( euler(1:9,ipoin),(/3,3/))          !rotation matrix
             CALL angeul(rm,angls) !returns Euler angles (in degrees)
           END IF
           WRITE (17,ERR=9999) angls(1:3)
         END DO
       END IF
       IF( ndofn == 8 )THEN
         DO ipoin = 1,npoin
           WRITE (17,ERR=9999) psia(1:2,ipoin)
         END DO
       END IF
     END IF
     IF(iwrit == 1) THEN   !write to ASCII file (.RSN)
       WRITE(lures,"(//5x,'Displacements at time step ',i10,5x,'Time ',e20.11/,    &
      &                5x,'Nnode',3x,'X-disp',6x,'Y-disp',6x,'Z-disp'/)",ERR=9999) &
      &                istep,ttime
       DO ipoin = 1,npoin  !displacements
         WRITE(lures,902,ERR=9999) label(ipoin),(coora(1:ndime,ipoin)-coord(1:ndime,ipoin))
       END DO
       IF(neulr > 0) THEN  !for local systmes
         IF(ndime == 2) THEN             !2-D problems
           WRITE(lures,930,ERR=9999)
           DO i=1,npoin
             IF(.NOT.naeul(i) )CYCLE
             WRITE(lures,933,ERR=9999) label(i),euler(1,i)*fac
           END DO
         ELSE                            !3-D problemas
           WRITE(lures,931,ERR=9999)
           DO ipoin = 1,npoin
             IF(.NOT.naeul(ipoin) )CYCLE
             rm = RESHAPE( euler(1:9,ipoin),(/3,3/))          !rotation matrix
             angls = 0d0                  !angles
             CALL angeul(rm,angls,.TRUE.) !returns Euler angles (in rads)
             WRITE(lures,"(i7,3e17.8)",ERR=9999) label(ipoin),angls(1:3)*fac
           END DO
           IF( ndofn == 8 )THEN
             WRITE(lures,"(//5x,'AdditionalDisplacements at time step ',i10,5x,'Time ',e20.11/,    &
      &                      5x,'Nnode',3x,'X-disp',6x,'Y-disp'/)",ERR=9999) &
      &                      istep,ttime
             DO ipoin = 1,npoin
              WRITE(lures,932,ERR=9999)label(ipoin), psia(1:2,ipoin)
             END DO
           END IF
         END IF
       END IF
     END IF
     IF( .NOT.static )THEN
       !**  velocity vector
       DO ipoin=1,npoin
         WRITE(17,ERR=9999) velnp(1:ndofn,ipoin)
       END DO
       IF(iwrit == 1) THEN  !write to ASCII file
         WRITE(lures,950,ERR=9999) istep,ttime
         IF(ndime == 2)  WRITE(lures,952,ERR=9999)
         IF(ndime == 3)  WRITE(lures,953,ERR=9999)
         DO ipoin = 1,npoin
           WRITE(lures,902,ERR=9999) label(ipoin),velnp(1:ndofn,ipoin)
         END DO
       END IF
       !**  aceleration vector
       ALLOCATE( disp3(ndofn,npoin) )   !get memory
       DO ipoin = 1,npoin   ! rearrange first
         DO idofn = 1,ndofn
           ieq = ifpre(idofn,ipoin)
           IF(ieq > 0) THEN
             disp3(idofn,ipoin) = acelr(ieq)
           ELSE
             disp3(idofn,ipoin) = 0d0
           END IF
         END DO
         WRITE(17,ERR=9999) disp3(1:ndofn,ipoin)
       END DO
       IF(iwrit == 1) THEN  ! write to ASCII file
         WRITE(lures,960,ERR=9999) istep,ttime
         IF(ndime == 2)  WRITE(lures,962,ERR=9999)
         IF(ndime == 3)  WRITE(lures,963,ERR=9999)
         DO ipoin=1,npoin
           WRITE(lures,902,ERR=9999) label(ipoin),disp3(1:ndofn,ipoin)
         END DO
       END IF
       DEALLOCATE( disp3 ) ! release memory
     END IF
     !** CONTACT values (gaps, pressures, friction-work)
     IF (numct > 0) CALL contac('OUTDY2',0,ttime=value,dtcal=dtime)
     !** temperatures at  nodal points
     IF (itemp) THEN
       DO ipoin=1,npoin
         WRITE(17,ERR=9999) tempe(:,ipoin)
       END DO
       tstra = (ttime - begtm)*tscal  ! time from the strategy start
       !WRITE(17,ERR=9999) (bqgen(ipoin)/tstra,ipoin=1,npoin) ! average heat generation rate
       ! the above is for wear algorithm
       IF(iwrit==1) THEN
         WRITE(3,961,ERR=9999) istep,ttime
         DO ipoin=1,npoin
           WRITE(3,902,ERR=9999) label(ipoin),(tempe(:,ipoin))
         END DO
       END IF
     END IF
   END IF

   !*** elemental variables (Gauss-points)
   IF(iwrit == 1 .AND. b2) WRITE(lures,"(/,10x,'Stresses ',/)",ERR=9999)
   IF(b1 .OR. b2) THEN
     CALL elemnt ('OUTDYN', deltc=dtime, ttime=value, flag1=b1, flag2=b2)
     IF(iwrit == 1) CALL flushf(lures)
     IF(b2) THEN
       CALL flushf(16)
       CALL flushf(17)
       WRITE(*,"(5x,'Results have been reported for postprocess')")  !screen
     END IF
     !IF (b2) CALL del_old_files(0)   !DO NOT erase GiD old files
   END IF

 RETURN
   902 FORMAT(5x,i5,6e13.5)
   930 FORMAT(//'  Nodal angles '/)
   933 FORMAT(i7,3x,e15.7)
   931 FORMAT(//'  Nodal cartesyan systems (Euler angles in rads) '/)
   932 FORMAT(1x,i10,1x,3f8.4)
   950 FORMAT(//5x,'Velocity at time step ',i10,5x,'Time ',e20.11/)
   952 FORMAT(  5x,'Nnode',3x,'X-vel',7x,'Y-vel',7x,'Omega'/)
   953 FORMAT(  5x,'Nnode',3x,'X-vel',7x,'Y-vel',7x,'Z-vel',6x, &
      &            'Omega-1',5x,'Omega-2',5x,'Omega-3'/)
   960 FORMAT(//5x,'Acceleration at time step ',i10,5x,'Time ',e20.11,/)
   961 FORMAT(//5X,'Temperatures at time step ',I10,5X,'Time ',e20.11,/ &
      &         5X,'Nnode',3X,'Temperatures')
   962 FORMAT(  5x,'Nnode',3x,'X-accel',6x,'Y-accel',6x,'Alpha',/)
   963 FORMAT(  5x,'nnode',3x,'X-accel',6x,'Y-accel',6x,'Z-accel',6x, &
      &            'Alpha-1',6x,'Alpha-2',6x,'Alpha-3'/)
   965 FORMAT(//5x,'Hydrostatic pressure at time step ',i10, &
      &         5x,'Time ',e20.11/)
   966 FORMAT(i10,e13.5)
  9999 CALL runen2('')
 END SUBROUTINE outdyn

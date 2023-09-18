 SUBROUTINE pstatic ( ircon, nrest, istop, actio, ecdis )
!
! Performs a Pseudo-STATIC analysis
!
!  Modifications to be done
!  verify computations when non convergences and step splitting
!  verify stack overflow for large problems (array assigments)
!
 USE param_db,ONLY : mich                       !Global parameters
 USE static_db
 USE ctrl_db, ONLY : istep,dtime,dtold,dtscal,ttime,begtm,endtm,nload,neq,numct,npoin,neulr,ncdlt,ndofn,ndime
 USE kinc_db, ONLY : nvelr,velor,nvfix,lcvel,ndepd,naris,npsdf,nesdf,ftsdf
 USE gvar_db, ONLY : updiv,newiv,static,elastic
 USE outp_db, ONLY : iwrit
 USE npo_db,  ONLY : veloc,ddisp,resid,fcont,coora,coorc,euler,naeul,loass,force,ymass,ifpre,loadv
 USE damp_db, ONLY : damp,updt_dval,gres_damp
 USE rfdp_db, ONLY : nrfdp
 USE sms_db, ONLY  : selective_mass_scaling,fibre_mass
 USE loa_db, ONLY  : ifloa
 USE nsld_db, ONLY : nnsld
 USE c_input, ONLY : openfi
 IMPLICIT NONE
 !dummy variables
 INTEGER(kind=4), INTENT(IN OUT) :: ircon, & !restar file number
                                    istop    !flag for errors in residual force computation
 INTEGER(kind=4), INTENT(IN) :: nrest,     & !frequency (seconds) for dumping
                                ecdis        !DOF to control output
 CHARACTER(len=*) :: actio                   !NEW, NSTRA0, RESTAR
 !--- Functions
 INTEGER(kind=4) :: chnode                   !changes from label to internal numeration
 REAL(kind=8) :: functs                      !time function
 !local variables
 REAL(kind=8), PARAMETER :: ref_val = 0.98d+00  !reference value to anulate velocities (kinetic damping)
 INTEGER(kind=4) :: itera,  & !iteration of present step
                    n1,     & !position of resultant loads
                    i,j,ii, & !index
                    pcdis,  & !control equation for debug output
                    np        !frequency to write control disp
 LOGICAL :: nocvg, & ! if maximum number of iterations (MITER) exceeded
            newid, & ! to generate ID array and arrays for automatic damping
            convg, & ! if convergence achieved
            modmp, & ! if damping modified
            newlm, & ! to re-generate lumped mass
            chkin, & ! if Kinetic Damping is still active (first step only)
            actrst   ! flag for dumping file for restar
 REAL(kind=8) :: rfnor_i,kem_i,     & !initial reference values
                 ratio,ke,          & !force and kinetic ratios for convergence check
                 velfi,             & !pseudo-time used for prescribed velocities increment
                 auxv,auxf,fac,     & !auxiliar value/factor
                 pi,pir,            & !pi and pi/(nides or nperi)
                 dl2,               & !square arc-length of actual DOFs only
                 dlamb                !incremental factor for load or prescribed displacements
  REAL(kind=8) ::  facts(nload), factv(nvelr) !action factors
  REAL(kind=8), ALLOCATABLE :: incdis(:,:)     !incremental displacements

  INTERFACE
    INCLUDE 'elemnt.h'
    INCLUDE 'coract.h'
    INCLUDE 'predic.h'
    INCLUDE 'contac.h'
    INCLUDE 'sxplit.h'
  END INTERFACE


   pi = 4d0*atan(1d0)   !3.14159 ...
   n1 = nload + 1       !position of resultant loads

   IF( TRIM(actio) /= 'RESTAR' )THEN       !new strategy
     dtime = (endtm - begtm )/ sstep       !formal time increment
     IF( ASSOCIATED(disax) )DEALLOCATE(disax)       !check and release
     ALLOCATE( disax(neq+1,3) )                     !get memory for previous incremental displacements
     IF( ASSOCIATED(resid_o) )DEALLOCATE(resid_o)   !check and release
     ALLOCATE( resid_o(ndime,npoin) )               !get memory for previous residual forces
     disax = 0d0                           !initializes incremental displacements
     resid_o = 0d0                         !initializes residual at previous step
     damp = 1d0 - 2d0*damp                 !compute damping factors
     rfnor = 0d0                           !initializes reference norm of residual forces
     kem   = 0d0                           !initializes reference (maximum) kinetic energy
     sitep = 1                             !initializes step for the strategy
     istep = istep + 1                     !Update number of global steps
     !Determine curve number associated to main action (Load or Velocity)
     IF( lovec > 0 .AND. nload > 0 )THEN       !loads are the principal action
       lovec = MIN(nload,lovec)                   !check and correct
     ELSE IF( lovec < 0 .AND. nvelr > 0 )THEN  !displacements are the principal action
       lovec = MAX(-nvelr,lovec)                  !check and correct
     ELSE IF( nload > 0 ) THEN                 !loads are the principal action
       lovec = 1                                  !if not declared by user use loads if exist
     ELSE IF( nvelr > 0 ) THEN                 !displacements are the principal action
       lovec = -1                                 !if not declared by user use displ if exist
     ELSE IF( .NOT.spback )THEN                !no action is present
       CALL runen3(' No solicitation (LOAD OR VEL SET) are present ')  !not possible to continue
     ELSE
       lovec = 0
     END IF
     ! compute load factor according to main action at the start
     IF( lovec > 0 )THEN
       olamb = functs(loass(lovec),ttime)*force(neq+1,lovec)     !load factor
     ELSE IF( lovec < 0 )THEN
       lambd = 0d0  !initializes displacement increment
       olamb = 0d0  !intializes load factor at previous increment
     ELSE
       lambd = ttime
       olamb = ttime
     END IF
     str_type = 0                  !default strategy
     jstep = 1                     !initializes strategy type step
   ELSE
     auxf = 1d0    !factor from previous increment
   END IF
   chkin = sitep == 1 .AND. autod  !kinetic damping for first step
   elastic = spback        !elastic springback

   IF( neulr == 9 )THEN                  !if local systems exist
     IF( ASSOCIATED(locsy) )DEALLOCATE(locsy) !check and release
     ALLOCATE(locsy(neulr,npoin))             !get memory for 'last converged' local systems
     DO i=1,npoin                !for each node
       locsy(:,i) = euler(:,i)   !keep last converged local system
     END DO
   END IF

   ALLOCATE( incdis(ndime,npoin) ) !incremental displacements (contact and damping)
   incdis = 0d0                    !initializes
   ALLOCATE( eqres(neq) )          !get auxiliar memory for non-equilibrated nodal forces
   ALLOCATE( smass(ndofn,npoin) )  !get auxiliar memory for modified mass matrix
   CALL slmass ( )                 !compute modified (artificial) mass matrix
   IF( pcont > 10 .OR. spback )THEN   !update output control equation for output
     i    = chnode(pcont/10)       !node (internal)
     j    =  MOD(pcont,10)         !nodal DOF
     IF( i > 0 .AND. j > 0 )THEN
       pcdis = ifpre(j,i)            !associated equation (must be active)
     ELSE
       pcdis = neq
     END IF
     np = nperi/100                !output frequency
     IF( pcdis > 0 ) CALL openfi(10,rest='disp_cp.prn')  !open file for control disp
   ELSE
     pcdis = 0                     !must be stored for restart
   END IF
   newid = .TRUE.                  !initializes flag to compute ID array

   DO ! L O O P   over each   S T E P
     rfnor_i = rfnor ; kem_i = kem      !keep present reference values in case of non-convergence
     convg = .FALSE. ; nocvg = .FALSE.  !Set exit flags to FALSE
     updiv = .FALSE.                    !do not update internal variables (default)
     ttime = ttime + dtime     !update ttime for present step (used for curve evaluation)

     ! update lumped mass matrix (critical DT computation)
     IF( ncdlt > 0 .AND. MOD(sitep,ncdlt) == 0 )THEN
       CALL slmass( )                   !update mass matrix leading to DDT=1
       newlm = .TRUE.                   !NEW Lumped mass computed
     END IF
     ! update lumped mass matrix when dependent nodes exists or mass recomputed
     IF (((ndepd+naris+nrfdp+nnsld ) > 0) .OR. newlm ) THEN !
       ymass = 0d0                                    !initializes
       CALL ensmal(ndofn*npoin,neq,ifpre(1,1),smass(1,1),ymass(1),npsdf(1),nesdf(1),ftsdf(1))
       IF( selective_mass_scaling )CALL fibre_mass(smass)    !recompute fiber mass
       DO i=1,neq                    !invert lumped mass matrix for optimum computations
         ymass(i) = 1d0/ymass(i)
       END DO
       newlm = .FALSE.                                !forget
     END IF
     ! update reference force vector when dependent nodes exists
     IF ((ndepd+naris+nrfdp ) > 0)  THEN
       IF (nload > 0) force(1:neq,1:nload)=0d0         !initializes force vectors
       DO i=1,nload                                    !recompute force vectors
         CALL ensve1(ndofn*npoin,ifpre(1,1),loadv(1,1,i),force(1,i),npsdf(1),nesdf(1),ftsdf(1))
       END DO
     END IF
     !Compute external equivalent forces (Force)
     IF (nload > 0 ) THEN          !if loads must be recomputed
       IF( str_type /= 1 )THEN
         ! for str_type = 1 (load control) present load vector is assumed proportional to Lambda
         !                  it means that load functions are dissable
         DO i=1,nload                                         !for each load set
           facts(i) = functs(loass(i),ttime)*force(neq+1,i)   !load factor
         END DO
         DO i=1,neq                              !for each active DOF
           force(i,n1) = DOT_PRODUCT(force(i,1:nload),facts(1:nload))  !compound force
         END DO
         IF (ifloa > 0) THEN  ! follower loads (not yet consistent)
           resid = 0d0
           CALL loadfl(ttime)  !Nonconservative loading
           resid = -resid
           CALL ensve1(ndofn*npoin,ifpre(1,1),resid(1,1),force(1,n1),npsdf(1),nesdf(1),ftsdf(1))
         END IF
         IF( lovec > 0 )THEN      !If forces are the DRIVING action
           lambd = facts(lovec)   !load factor
           dlamb = lambd - olamb  !increment in load factor
         END IF
       END IF
     END IF
     !   compute prescribed velocities
     IF (nvelr > 0) THEN               !if prescribed velocities exist
       auxv = ttime - dtime/2d0        !pseudo time at step mid-point
       DO i=1,nvelr                                          !for each VEL SET
         factv(i) = functs(lcvel(i),auxv)*velor(nvfix+1,i)     !velocity factor
       END DO
       DO i=1,nvfix                                         !for each prescribed DOF
         velor(i,nvelr+1) = DOT_PRODUCT(velor(i,1:nvelr),factv(1:nvelr))     !compound velocity
       END DO
       IF( lovec < 0 )THEN      !If displacements are the DRIVING action
         lambd = lambd + factv(-lovec) * dtime   !accumulated prescribed displacement
         dlamb = lambd - olamb                   !increment in prescribed displacement
       END IF
     END IF
     ! compute factor for prescribed velocities
     IF( str_type == 2 ) THEN      !for displacemente strategy
       pir = pi/nperi                !auxiliar value for str_type = 2
     ELSE IF( str_type == 1 ) THEN !for load strategy this have no sense
       pir = 1d0                     !set to 1 but non-consistent
     ELSE IF( jstep  == 1 ) THEN   !first step of the standard strategy
       pir = pi/2d0/nides            !auxiliar value for str_type = 0
     ELSE                          !rest of the steps for standard strategy
       pir = 1d0                     !pseudo time increment at first iteration (standard step)
     END IF
     auxv    = dtime*pir       !auxiliar value
     ! compute prediction
     IF( str_type == 0 )THEN                      !standard strategy
       CALL predic(jstep,neq,disax,ddisp,dlamb,ipred)  !predict initial increment
       DO i=1,neq
         veloc(i) = 0d0                            !initializes velocities
       END DO
     ELSE IF( str_type == 1 ) THEN                !non-standard load strategy (arc-length)
       !no prediction used, use DDISP from previous successfull step
       DO i=1,neq
         ddisp(i) = -disax(i,1)                    !use last increment
       END DO
       ! a new value based on number of iterations could be used
       !CALL predic(jstep+1,neq,disax,ddisp,dlamb,ipred)  !predict initial increment
       DO i=1,neq
         veloc(i) = 0d0                            !initializes velocities
       END DO
     ELSE                                         !non-standard strategy
       fac = auxf*pir
       DO i=1,neq
         veloc(i) = -disax(i,1)*fac                  !initial velocity
         ddisp(i) = veloc(i)                         !initial displacement increment
         disax(i,3) = ddisp(i)                       !initializes incremental displacements
       END DO
     END IF
     velfi = auxv                                   !used as factor for prescribed displacements in the step
     CALL coract(velfi,velor(1:nvfix,nvelr+1),ddisp,coora,euler) !update coordinates => coora - euler
     DO i=1,npoin
       incdis(:,i) = coora(:,i) - coorc(:,i)        ! incremental displacement
     END DO
     ! automatic damping computations   !next IF could be changed to SUCCESSFULL STEP instead
     IF( autod  )THEN     !not for displacement strategy
       IF( jstep > 1 )THEN     !not for displacement strategy
         CALL updt_dval (ndime,pmin,pmax,modmp,fdamp,resid,resid_o,incdis,smass)  !updates damping values
         IF(modmp)THEN
           CALL gres_damp (ndime, ndofn, neq, ifpre)       !regenerates damping
           DO i=1,neq
            damp(i) = 1d0 - 2d0*damp(i)                    !compute damping factors
           END DO
         END IF
       END IF
       IF( sitep == 1 .OR. jstep > 1)THEN
         DO i=1,npoin
           resid_o(:,i) = resid(1:ndime,i)   !keep residual forces
         END DO
       END IF
     END IF
     IF( spback )THEN
       rfnor = 0
       DO i=1,npoin
         resid_o(:,i) = resid(1:ndime,i)   !keep residual forces for springback
         rfnor = rfnor + DOT_PRODUCT(resid(:,i),resid(:,i))
       END DO
       rfnor = SQRT(rfnor)
     END IF
     !
     itera = 1                                        !initialized to first iteration
     DO          ! I T E R A T I O N      L O O P
       newiv = .FALSE.  !initializes flag that indicates if internal variables change during the step
       IF( str_type == 2 .AND. MOD(itera,nides) == 0 )updiv = .TRUE. !update Int-Var for displacement strategy
       !to update internal variables in STR=1, see what to do with arc-length
       !IF( str_type == 1 .AND. MOD(itera,nides) == 0 )updiv = .TRUE. !update Int-Var for displacement strategy
       ! compute equivalent nodal forces
       CALL timing(14,1)
       DO i=1,npoin
         resid(:,i) = 0d0 !initializes residual
       END DO
       CALL elemnt('RESVPL', deltc=1d0, istop=istop, ttime= ttime)
       IF( spback .AND. itera < nperi)THEN
         fac = SIN(pi/2d0/nperi*itera) - 1d0
         DO i=1,npoin
           resid(1:ndime,i) = resid_o(:,i)*fac + resid(1:ndime,i)   !modify residual forces
         END DO
       END IF
       CALL timing(14,2)
       !           Calculation of contact forces
       IF(numct > 0)THEN
         CALL timing(16,1)
         DO i=1,npoin
           fcont(:,i) = 0d0
           incdis(:,i) = coora(:,i) - coorc(:,i)
         END DO
         !OK ttime controls begin-end contact activation
         !OK dtime controls maximum iterative displacement (DISMA) and penalty
         !OK nvfix must be set to 0 to avoid inclusion of velor in DISMA
         !TEST Problems may appear at the beginning of the step
         !OK VELOC is replaced by incremental displacement (INCDIS)
         !there are difference between strategies
         !Use flag UPDIV to update internal variables !!!
         CALL contac('FORCES',iwrit,dtcal=1d0, ttime=ttime, velnp=incdis)    !toutd=toutd,
         CALL timing(16,2)
       END IF
       CALL timing(13,1)
       IF( .NOT.updiv  )THEN       !restore coordinates and local system to converged values
         !coordinates
         DO i=1,npoin
           coora(:,i) = coorc(:,i)
         END DO
         !local systems
         IF(neulr == 9) THEN      !ndime == 3
           DO i=1,npoin
             IF( .NOT. naeul(i) )CYCLE
             euler(:,i) = locsy(:,i)
           END DO
         END IF
       ELSE                     ! update converged coordinates and local systems (str_type = 2)
         !coordinates
         DO i=1,npoin
           coorc(:,i) = coora(:,i)
         END DO
         !local systems
         IF(neulr == 9) THEN      !ndime == 3
           DO i=1,npoin
             IF( .NOT. naeul(i) )CYCLE
             locsy(:,i) = euler(:,i)
           END DO
         END IF
         DO i=1,neq
           disax(i,3) = disax(i,3) + ddisp(i)    !update step displacement
           ddisp(i) = 0d0                              !and initializes
         END DO
         velfi = 0d0                              !initializes displacement factor
         updiv = .FALSE.                          !back to .FALSE.
       END IF
       ! recompute VELFI (factor for prescribed Velocities-Displacements)
       SELECT CASE (str_type)  !according to the strategy type
       CASE (0) !standard strategy
         velfi = dtime  !standard value (total step increment)
         IF(  jstep == 1 .AND. itera < nides  ) velfi = dtime*SIN(pir*itera) !for first step only
       CASE (1) !for load strategy (arc-length) although not consistent
         velfi = dtime
       CASE (2) !for displacement strategy
         IF( itera <= nperi/2 ) THEN         !during prescribed velocities
           velfi = velfi + auxv*COS(pir*itera)   !compute factor
         ELSE                                !after prescribed velocities
           velfi = 0d0                           !set to 0
         END IF
       END SELECT
       CALL sxplit(newid,velfi,eqres,str_type,lambd,chkin,dl2)  ! standard central difference scheme
       CALL timing(13,2)
       ! compute present kinetic energy and update reference value
       ke = 0d0                            !initializes
       DO i=1,neq          ! for each active DOF
         ke = ke + veloc(i)**2/ymass(i)
       END DO
       ! compute norms (kin-energy & forces) and update reference value
       ke = SQRT(ke)                       !kinetic energy ratio (squared)
       IF( ke > kem )  kem = ke            !keep maximum for comparison
       ratio = SQRT(DOT_PRODUCT(eqres,eqres))  !norm of residual forces
       IF( .NOT.spback) THEN
         IF( ratio > rfnor )THEN  !compare with present value
           IF( rfnor == 0d0 ) THEN     !not yet defined
             rfnor = ratio                !update reference value
           ELSE IF( MOD(itera,nides) == 0 ) THEN  !modify only at selected iterations
             rfnor = MIN(ratio,2d0*rfnor) !update reference value n
           END IF
         END IF
       END IF
       ratio = ratio/rfnor                     !ratio for convergence check
       ke    = ke   /kem                       !ratio for convergence check
       IF( chkin )THEN  !IF kinetic damping active
         IF( ke < ref_val )THEN  !maximum arrived
           IF( itera > pmin/4 )THEN
             WRITE(55,*)'maximum kinetic energy for ITERA',itera
             DO i=1,neq
               veloc = 0d0    !zero velocities
             END DO
             ! modify damping
             WRITE(* ,*)'maximum kinetic energy for ITERA',itera
             ii = MIN(ii,pmax/4)
             CALL updt_dval (ndime,pmin,pmax,modmp,fdamp,itera=ii)  !updates damping values
             IF(modmp)THEN
               CALL gres_damp (ndime, ndofn, neq, ifpre)    !regenerates damping
               DO i=1,neq
                 damp(i) = 1d0 - 2d0*damp(i)                !compute damping factors
               END DO
             END IF
             chkin = .FALSE. !not any longer
           END IF
         END IF
       END IF
       ! check convergence
       IF( itera >= nperi )THEN !for standard strategy or after NPERI iterations
         SELECT CASE (ittol)        !convergence check
         CASE (1)                          !residual forces
           convg = ratio < tol
         CASE (2)                          !kinetic energy
           convg = ke    < tol
         CASE (3)
           convg =  ke*ratio < tol*tol     !both residual forces + kinetic
         END SELECT
       END IF
       ! check non-convergence
       nocvg = itera == miter                                !check iteration limit
       ! print present reference values to Screen
       IF(itera == 1) THEN
         IF( spback )THEN
           WRITE(6, "(1X,'Istep=',I6,'  Time=',E11.4,' FR=',f9.6,'  KR=',f9.6)") istep,ttime,ratio,ke
         ELSE
           WRITE(6, "(1X,'Istep=',I6,'  DelT=',E10.3,'  Time=',E11.4,' FR=',f9.6,'  KR=',f9.6)") istep,dtime,ttime,ratio,ke
           WRITE(55,"(1X,'Istep=',I6,'  DelT=',E10.3,'  Time=',E11.4,' FR=',f9.6,'  KR=',f9.6)") istep,dtime,ttime,ratio,ke
         END IF
       END IF
       IF (MOD(itera,nides) == 0 ) THEN
         IF( spback )THEN
           WRITE(6,"('+  Itera=',I6,'  cp_disp=',e13.4,'  FR=',f9.6,' KR=',f9.6)") itera,ddisp(pcdis),ratio,ke !print iterative screen
         ELSE
           WRITE(6,"('+  Itera=',I6,'  Lambda=',e13.4,'  FR=',f9.6,' KR=',f9.6)") itera,lambd,ratio,ke !print iterative screen
           WRITE(55,"('  Itera=',I6,'  Lambda=',e13.4,'  FR=',f9.6,' KR=',f9.6)") itera,lambd,ratio,ke !print iterative screen
         END IF
       END IF
       ! print reference values to control-point file
       IF(pcdis > 0 .AND. MOD(itera,np) == 0)THEN
         SELECT CASE (str_type)
         CASE (0:1)
           WRITE(10,"(i6,3e20.8,i6)") itera,ddisp(pcdis),ke,ratio,istep  !debug value
         CASE (2)
           WRITE(10,"(i6,3e20.8,i4)") itera,disax(pcdis,3)+ddisp(pcdis),ke,ratio,istep  !debug value
         END SELECT
       END IF

       IF (istop /= 0 .OR. nocvg .OR. convg  ) EXIT        !check End of iteration loop
       itera = itera + 1                                   !update iteration

     END DO  ! I T E R A T I O N    L O O P
     WRITE(58,*) rfnor,kem

     !-----------------------------------------------------------------------------------
     !I F   N O   C O N V E R G E N C E   C H E C K   I F  R E L A X I N G
     IF( nocvg .AND. itera == miter) THEN ! if maximum number of iterations exceeded
       ! but convergence is not so far GO AHEAD with present state
       SELECT CASE (ittol)        !convergence check
       CASE (1)                        !residual forces
         convg = ratio < tolm
       CASE (2)                        !kinetic energy
         convg = ke    < tolm
       CASE (3)                        !both residual forces + kinetic
         convg = ke*ratio <  tolm*tol
       END SELECT
       ! for alternative strategies relax convergence and go ahead
       IF( .NOT.convg .AND. str_type /= 0 )THEN
          IF( ratio*ke < tolm )THEN
            convg = .TRUE.
            IF( (jstep+1) / nss  == nds ) nds = 1 !increase number of desp. steps
          END IF
       END IF
     END IF

     !-----------------------------------------------------------------------------------

     IF( convg ) THEN        ! I F   C O N V E R G E N C E   A C H I E V E D
       ! update coordinates, local systems and internal variables (last converged values)
       DO i=1,npoin
         coorc(:,i) = coora(:,i)
       END DO
       IF(neulr == 9) THEN      !ndime == 3
         DO i=1,npoin                !for each node
           IF( .NOT. naeul(i) )CYCLE
           locsy(:,i) = euler(:,i)
         END DO
       END IF
       updiv = .TRUE.
       IF( newiv )THEN   !update internal variables if necessary
         CALL timing(14,1)
         DO i=1,npoin
           resid(:,i) = 0d0
         END DO
         CALL elemnt('RESVPL', istop=istop, ttime= ttime) !, deltc=1d0
         CALL timing(14,2)
       END IF
       IF(numct > 0)THEN  ! update sliding setup (newiv ??)
         CALL timing(16,1)
         DO i=1,npoin
           fcont(:,i) = 0d0
         END DO
         CALL contac('FORCES',iwrit,dtcal=1d0, ttime=ttime, velnp=incdis)  !does not modify incdis
         CALL timing(16,2)
       END IF
       !           Results output
       CALL timing(15,1)
       CALL outdyn(.FALSE.,.FALSE., ecdis)
       CALL timing(15,2)
       ! according to present strategy used update control parameters
       IF( str_type == 0 )THEN      !standard strategy
         WRITE(55,"('Convergence in Istep =',i5,4X,' in Iterations=',I8)") sitep,itera !print to report file
         sitep = sitep + 1          !update step for the loop
       ELSE IF (str_type == 1)THEN  !load control strategy
         WRITE(55,"('Convergence in Istep =',i5,4X,' in Iterations=',I8,' Load Strategy')") sitep,itera
         dlamb = lambd - olamb      !recompute present load increment
         sitep = sitep + 1          !update step for the loop
       ELSE IF (str_type == 2)THEN  !displacement control
         WRITE(55,"(' step ',i3,' successfull with displacement strategy , at STEP ',i4)")jstep,sitep
         IF( MOD(jstep,nss) == 0 )THEN  !according to step
           sitep = sitep + 1        !update step for the loop
           IF( jstep/nss == nds )THEN  !end with this strategy
             str_type = 0           !back to standard
             jstep = sitep -1       !allow prediction
             dtime = dtime * nss    !recompute time increment
           END IF
         END IF
         DO i=1,neq
           ddisp(i) = disax(i,3)
         END DO
         auxf = 1d0                 !keep increment from previous step
       END IF
       istep = istep + 1      !update number of (global) steps
       jstep = jstep + 1      !update number of strategy steps
       ! update incremental displacements for prediction (DISABLE)
       DO i=1,neq
         disax(i,3) = disax(i,2) - ddisp(i)
         disax(i,2) = disax(i,1) - ddisp(i)
         disax(i,1) = - ddisp(i)
       END DO
       ! its difficult to define DLAMB due to multiple LOAD and PRESCRIBED velocities functions
       disax(neq+1,3) = disax(neq+1,2) - dlamb
       disax(neq+1,2) = disax(neq+1,1) - dlamb
       disax(neq+1,1) = - dlamb
       olamb = lambd          !keep old load factor
       IF(sitep > sstep) static = .FALSE.  ! do not dump STATIC_DB at the end

     !----------------------------------------------  N O   C O N V E R G E N C E

     !       CHANGE TO   A L T E R N A T I V E    S T R A T E G I E S   IF POSSIBLE

     ELSE IF( .NOT.sstop .AND. sitep > 1 .AND. nload*nvelr == 0 .AND. str_type == 0 )THEN  !if only one type of actions
       WRITE( 6,"(' No Convergence at Istep=',I5,'  Time=',E11.4,'  f-rat=',f9.6,'  k-rat=',f9.6)") istep,ttime,ratio,ke
       WRITE(55,"(' No Convergence at Istep=',I5,'  Time=',E11.4,'  f-rat=',f9.6,'  k-rat=',f9.6)") istep,ttime,ratio,ke
       rfnor= rfnor_i ; kem  = kem_i   !recover reference values
       llamb  = lambd             !keep failed load factor
       ltime = ttime              !keep failed ttime
       ttime = ttime - dtime      !update ttime back
       IF( nload /= 0 )THEN  !allow global softening
         WRITE(6 ,"(' New strategy with load strategy, at STEP ',i4)") sitep
         WRITE(55,"(' New strategy with load strategy, at STEP ',i4)") sitep
         str_type = 1                     !change strategy to variable load (arc-length)
         dtime = dtime          !keep time increment although not useful any longer
         jstep = 1              !initializes step
         auxf = 1d0             !???
         newid = .TRUE.         !compute auxiliar arrays is sxplit
       ELSE                 !allow degradation at constant displacements
         WRITE(55,"(' New strategy with velocitites only, at STEP ',i4)") sitep
         WRITE(6 ,"(' New strategy with velocitites only, at STEP ',i4)") sitep
         str_type = 2           !change to displacement strategy
         lambd = olamb          !back to previous value
         dtime = dtime/nss      !change time increment into NSS parts
         auxf = 1d0/nss         !modify increment in displacement
         jstep = 1              !initializes strategy step
       END IF
       !      restore coordinates and local system to converged values
       DO i=1,npoin
         coora(:,i) = coorc(:,i)
       END DO
       IF(neulr == 9) THEN      !ndime == 3
         DO i=1,npoin                !for each node
           IF( .NOT. naeul(i) )CYCLE
           euler(:,i) = locsy(:,i)
         END DO
       END IF
       istop = 0
       CYCLE
     ELSE
       ttime = ttime - dtime  !update ttime back (for dumping file)
       istop = -2
     END IF
     !----------------------------------------------------------------------
     !           Restart file dumping check
     IF (nrest > 0) THEN
       CALL timerst(.FALSE.,nrest,actrst)
     ELSE IF (nrest < 0) THEN
       actrst = MOD(sitep,-nrest)==0
     ELSE
       actrst = .FALSE.
     END IF

     IF( str_type /= 0 .AND. jstep == 1 ) actrst = .TRUE. !dumping when strategy change

     IF (actrst .OR. (istop /= 0) .OR. .NOT.convg .OR. sitep > sstep) THEN
       CALL timing(15,1)
       CALL dumpin(ircon, sstep < sitep , .FALSE. ,.FALSE.)
       CALL timing(15,2)
       IF ((nrest > 0)) CALL timerst(.TRUE.,nrest,actrst)   !Initialize time for writting restart file
     END IF
     IF( (sitep > sstep) .OR. (istop == 1) .OR. .NOT.convg ) EXIT   !exit loop

   END DO  ! G L O B A L   S T E P   L O O P
   ! -------------------------------------------------------
   IF (istop /= 0) THEN ! -> there was no convergence
     IF( newiv )THEN   !update internal variables if necessary to show non-converged values
       CALL timing(14,1)
       DO i=1,npoin
         resid(:,i) = 0d0
       END DO
       updiv = .TRUE.
       CALL elemnt('RESVPL', istop=istop, ttime= ttime) !, deltc=1d0
       CALL timing(14,2)
     END IF
     ttime = ttime + dtime  !update ttime again (for output file )
     CALL outdyn(.FALSE.,.TRUE., ecdis)  ! final state
     WRITE(6 ,900)
     WRITE(55,900,err=9999)
   END IF

   DEALLOCATE( disax,eqres,smass)               !release memory
   IF( neulr == 9 ) DEALLOCATE(locsy)      !release memory
   actio = 'NSTRA0'                        !for exit
   IF(pcdis > 0) CLOSE(10)  !close displacement debug file
 RETURN
 9999 CALL runen2('ERROR WHILE WRITING "DEB" FILE.')
 ! ++++FORMATS ++++++++++++++++++
  900 FORMAT(/,12x,' NO CONVERGENCE IN ITERATIVE ALGORITHM : ', &
     & 'PROGRAM STOPPED',/,17x,' FINAL GEOMETRY SAVED FOR POSTPROCESSING!')
 END SUBROUTINE pstatic

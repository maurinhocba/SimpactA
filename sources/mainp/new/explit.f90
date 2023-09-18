SUBROUTINE explit(flag)
!********************************************************************
!*** time stepping routine for DYNAMIC strategy
!********************************************************************
!$ USE omp_lib
USE ctrl_db,ONLY: a0, dtime, dtold, dtrec, ifixd, ifunc, nacce, numct, ndime,   &
                  ndimc, ndofn, neulr, nload, npoin, ttime, tbega, xbase, neq, lumped
USE damp_db
USE kinc_db
USE npo_db
USE rfdp_db, ONLY : nrfdp
USE sms_db
USE nsld_db, ONLY : nnsld
IMPLICIT NONE

  !-- Dummy variables
  LOGICAL,INTENT(IN OUT):: flag
  !--- Function
  REAL(kind=8):: functs
  !--- Local variables
  INTEGER(kind=4):: i,j,k,l, n1, ieq, idofn, ipoin
  INTEGER(kind=4),SAVE:: istep=0
  INTEGER(kind=4),POINTER,SAVE:: id(:,:)
  REAL(kind=8):: dt2, factt, functa, facts(nload), factv(nvelr)
  INTERFACE
    INCLUDE 'dampin.h'
    INCLUDE 'nvdamp.h'
    !INCLUDE 'rigwal.h'
    INCLUDE 'colsol.h'
    INCLUDE 'coract.h'
    INCLUDE 'velnpo.h'
  END INTERFACE

  dt2 = (dtime + dtold)/2d0           !average time increment
  istep = istep + 1                   !update step

  !CALL timingd(1,1)
  !generate array ID to assemble contact forces
  IF (flag .AND. numct > 0) THEN
    IF (ndime == ndofn) THEN
      IF (.NOT.ASSOCIATED(id,ifpre)) id=>ifpre
    ELSE
      IF (ASSOCIATED(id)) DEALLOCATE(id)
      ALLOCATE(id(ndimc,npoin))
      DO i=1,npoin
        id(:,i) = ifpre(1:ndimc,i)
      END DO
    END IF
    flag = .FALSE.
  END IF
  ! update lumped mass matrix and force vector when dependant nodes exists
  IF (((ndepd+naris+nrfdp+nnsld ) > 0) .AND. (MOD(istep,20) == 0)) THEN !
    IF (nload >0) force(1:neq,1:nload)=0d0
    DO i=1,nload
      CALL ensve1(ndofn*npoin,ifpre(1,1),loadv(1,1,i),force(1,i),npsdf(1),nesdf(1),ftsdf(1))
    END DO
    IF( lumped )THEN
      ymass = 0d0
      CALL ensmal(ndofn*npoin,neq,ifpre(1,1),emass(1,1),ymass(1),npsdf(1),nesdf(1),ftsdf(1))
      IF( selective_mass_scaling )CALL fibre_mass(emass)
      DO i=1,neq
        ymass(i) = 1d0/ymass(i)
      END DO
    END IF
  END IF
  !CALL timingd(1,2)


  !CALL timingd(2,1)
  ddisp = 0d0  !initializes forces
  !Add internal equivalente forces (resid)
  CALL ensve1(ndofn*npoin,ifpre(1,1),resid(1,1),ddisp(1),npsdf(1),nesdf(1),ftsdf(1))
  !CALL timingd(2,2)

  !Compute and add external equivalent forces (Force)
  IF (nload > 0) THEN
    !CALL timingd(3,1)
    DO i=1,nload
      facts(i) = functs(loass(i),ttime)*force(neq+1,i)
    END DO
    n1 = nload + 1
    DO ieq=1,neq
      force(ieq,n1) = DOT_PRODUCT(force(ieq,1:nload),facts(1:nload))
      ddisp(ieq) = ddisp(ieq) - force(ieq,n1)
    END DO
    !CALL timingd(3,2)
  END IF
  !   compute prescribed velocities
  IF (nvelr > 0) THEN
    !CALL timingd(4,1)
    DO i=1,nvelr
      factv(i) = functs(lcvel(i),ttime)*velor(nvfix+1,i)
    END DO
    DO i=1,nvfix
      velor(i,nvelr+1) = DOT_PRODUCT(velor(i,1:nvelr),factv(1:nvelr))
    END DO
    !CALL timingd(4,2)
  END IF

  IF (numct > 0) THEN
    !CALL timingd(5,1)
    CALL ensve1(ndimc*npoin,id(1,1),fcont(1,1),ddisp(1),npsdf(1),nesdf(1),ftsdf(1))
    !CALL timingd(5,2)
  END IF
  !Modify previous velocities due to damping
  IF (ndamp > 0) THEN
    !CALL timingd(6,1)
    SELECT CASE (dtype)
    CASE (1) ! standard viscous
      CALL dampin(neq,damp,veloc,dtime)
    CASE (2) ! non-viscous
      CALL nvdamp(neq,damp,veloc,ddisp)
    CASE DEFAULT
    END SELECT
    !CALL timingd(6,2)
  END IF

  !If imposed accelerations modify forces
  IF (ifunc == 0) THEN
    DO idofn=1,ndime
      IF ((ifixd == 0) .OR. (ifixd == idofn)) THEN
        factt = functa(acceg(1:nacce,idofn),dtrec,nacce,ttime-tbega)*a0
        DO ipoin=1,npoin
          ieq = ifpre(idofn,ipoin)
          IF (ieq > 0) ddisp(ieq) = ddisp(ieq) + factt/ymass(ieq)
        END DO
        !Store ground acceleration, velocity & displacement
        xbase(1) = factt                     ! acceleration
        xbase(2) = xbase(2) + factt*dt2      ! velocity
        xbase(3) = xbase(3) + xbase(2)*dtime ! displacement
      END IF
    END DO
  END IF

  !Integrate equations
  ! compute standard accelerations
  !CALL timingd(7,1)
  IF( lumped ) THEN
    !$OMP DO
    DO ieq = 1,neq
      acelr(ieq) = -ddisp(ieq)*ymass(ieq)
    END DO
    !$OMP END DO
  ELSE
    acelr = -ddisp
    CALL colsol(ymass(:),maxav(:),neq,2,58,0,0,v=acelr(:))
  END IF
  !CALL timingd(7,2)
  IF( selective_mass_scaling )THEN    ! loop over fibres
    !CALL timingd(8,1)
    sms_faux = 0d0 !initializes non-equilibrated forces
    DO i=1,sms_nf  !for each fibre
      ! first loop to compute forces on the middle surface
      DO j=1,sms_nnf  !for each node on the fibre
        k = sms_fibre(j,i)           !node number
        IF( k == 0 )EXIT             !if all processed exit loop
        DO l=1,ndime                 !for each space dimension
          ieq = ifpre(l,k)           !global DOF
          IF( ieq < 1 )CYCLE         !if not an active dof go to next
          sms_faux(l,i) = sms_faux(l,i) - ddisp(ieq)  !add force
        END DO
      END DO
      ! compute accelerations on the middle surface
      DO l=1,ndime
        sms_faux(l,i) = sms_faux(l,i) / sms_mf(l,i)   !node accelerations (already scaled by (1-alpha)
      END DO
      ! second loop to compute accelerations on the nodes
      DO j=1,sms_nnf  !for each node on the fibre
        k = sms_fibre(j,i)           !node number
        IF( k == 0 )EXIT             !if all processed exit loop
        DO l=1,ndime                 !for each space dimension
          ieq = ifpre(l,k)           !global DOF
          IF( ieq < 1 )CYCLE         !if not an active DOF go to next
          acelr(ieq) = acelr(ieq) * sms_mf(0,i) + sms_faux(l,i)  !modify acceleration
        END DO
      END DO
    END DO
    !CALL timingd(8,2)
  END IF

  !Integrate velocities and displacements
  !CALL timingd(9,1)
 !$OMP DO
  DO ieq = 1,neq
    veloc(ieq) = veloc(ieq) + acelr(ieq)*dt2
  END DO
 !$OMP END DO
 !$OMP DO
  DO ieq = 1,neq
    ddisp(ieq) = veloc(ieq)*dtime
  END DO
 !$OMP END DO
  !CALL timingd(9,2)

  !Update nodal coordinates
  !CALL timingd(10,1)
 !$OMP DO
  DO i=1,npoin
    coorc(:,i) = coora(:,i)
  END DO
 !$OMP END DO
  CALL coract(dtime,velor(1:nvfix,nvelr+1),ddisp,coora,euler)
  !WRITE(58,"(i5,3e20.12)")(i,coora(:,i)-coorc(:,i),i=1,npoin)
  !Rearrange velocity vector (necessary for some element routines)
  CALL velnpo(ndofn, npoin, nvelr, ifpre, nesdf, npsdf, ftsdf, velor, veloc, velnp)
  !CALL timingd(10,2)

  !IF (nrpla > 0) CALL rigwal(nrpla,ndime,npoin,ifpre,veloc,coora,vplan,pplan)

  dtold = dtime

RETURN
END SUBROUTINE explit

SUBROUTINE ensve1(nvarl,lm,locvc,glovc,npsdf,nesdf,ftsdf)
!*************************************************************************
! Assemble a local vector into a global vector
!*************************************************************************
USE kinc_db, ONLY : nn
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: nvarl, lm(nvarl), npsdf(*), nesdf(*)
  REAL(kind=8),INTENT(IN):: locvc(nvarl), ftsdf(*)
  REAL(kind=8),INTENT(IN OUT):: glovc(*)
  !--- Local variables
  INTEGER(kind=4):: i, j, k, nec, ib, ie

  DO i = 1,nvarl                                           !for each value
    nec = lm(i)                                            !assoc. equation
    SELECT CASE (nec)
    CASE (1:)                                              !if active dof
      glovc(nec) = glovc(nec) + locvc(i)                   !sums on global
    CASE (-nn:-1)                                          !if slave dof
      ib = npsdf(-nec)                                     !first post.
      ie = npsdf(-nec+1)-1                                 !last post.
      DO j=ib,ie                                           !for each master
        k = nesdf(j)                                       !assoc. equation
        IF(k > 0) glovc(k) = glovc(k) + locvc(i)*ftsdf(j)  !sums on global
      END DO
    END SELECT
  END DO

RETURN
END SUBROUTINE ensve1
!!*****Next 2 routines to compute CPU times in EXPLIT *****************
!!*********************************************************************
!      SUBROUTINE timingd (k,ind)
!
!      ! initializes and add CPU times for different tasks
!      ! for development and debug
!
!      USE outp_db, ONLY : timed
!      IMPLICIT NONE
!      INTEGER (kind=4),INTENT(IN) :: k,ind
!
!      REAL (kind=8) :: tcpu
!
!      CALL timuse (tcpu)                          !get clock time
!      IF (ind == 1) THEN                       !start task
!        timed(k+10) = tcpu                         !initial time
!      ELSE                                     !end task
!        timed(k) = timed(k) + tcpu - timed(k+10)     !adds elapsed time
!      END IF
!      RETURN
!
!      END SUBROUTINE timingd
!      SUBROUTINE debtime (time)
!!*********************************************************************
!!     writes cpu times
!!*********************************************************************
!      USE gvar_db, ONLY : actchk
!      IMPLICIT NONE
!      REAL (kind=8),INTENT(IN OUT) :: time(1:20)
!
!      INTEGER (kind=4) :: i
!      REAL (kind=8) :: timet
!
!      IF(actchk)RETURN
!
!      timet = SUM(time(1:10))
!      DO i=1,10
!        time(i+10) = time(i)*100d0/timet
!      END DO
!
!      WRITE (58,4) timet,100.0,(time(i),time(i+10),i=1,10)
!      RETURN
!   4  FORMAT(///,'T I M I N G   I N F O R M A T I O N',//,&
!     & 5X,'ACTION                           CPU(SEC)         %',//,&
!     & 5X,'total time devoted to EXPLIT..',10X,E14.4,2X,F6.2,' %',//,& ! 0
!     &15X,'initial & mass recomputation..',E14.4,2X,F6.2,' %',//,&     ! 1
!     &15X,'Assemblange of internal forces',E14.4,2X,F6.2,' %',//,&     ! 2
!     &15X,'External forces ..............',E14.4,2X,F6.2,' %',//,&     ! 3
!     &15X,'Prescribed velocities.........',E14.4,2X,F6.2,' %',//,&     ! 4
!     &15X,'Contact forces ...............',E14.4,2X,F6.2,' %',//,&     ! 5
!     &15X,'Damping ......................',E14.4,2X,F6.2,' %',//,&     ! 6
!     &15X,'Accel. integration loop ......',E14.4,2X,F6.2,' %',//,&     ! 7
!     &15X,'SMS computations .............',E14.4,2X,F6.2,' %',//,&     ! 8
!     &15X,'Vel-Displac. integration loop.',E14.4,2X,F6.2,' %',//,&     ! 9
!     &15X,'Coord & veloc update .........',E14.4,2X,F6.2,' %',//)      ! 10
!
!      END SUBROUTINE debtime

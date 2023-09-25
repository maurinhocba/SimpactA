SUBROUTINE sxplit(flag,velfi,eqres,str_type,lambd,chkin,dl2)
!********************************************************************
!*** time stepping routine for STATIC strategy
!********************************************************************
!$ USE omp_lib
USE ctrl_db,ONLY: numct, ndime, ndofn, neulr, nload, npoin, neq, ndimc
USE damp_db
USE kinc_db
USE npo_db
USE sms_db
IMPLICIT NONE

  !-- Dummy variables
  REAL(kind=8), INTENT (IN OUT) :: dl2       !norm of incremental displacements (squared)
  REAL(kind=8), INTENT (IN OUT) :: velfi     !factor for prescribed displacement increment
  REAL(kind=8), INTENT (IN OUT) :: eqres(:), &  !non-equilibrated nodal forces
                                   lambd        !present load factor
  INTEGER(kind=4), INTENT(IN) :: str_type
  LOGICAL, INTENT(IN OUT) :: flag,chkin

  !--- Local variables
  INTEGER(kind=4):: i,j,k,l, n1, ieq
  REAL(kind=8) :: factor,aux
  INTEGER(kind=4),POINTER,SAVE:: id(:,:),lind(:)
  INTEGER(kind=4),SAVE :: nf
  REAL(kind=8),POINTER,SAVE:: fn(:)
  REAL(kind=8),SAVE :: lambda

  INTERFACE
    INCLUDE 'dampin.h'
    !INCLUDE 'rigwal.h'
    INCLUDE 'coract.h'
    !INCLUDE 'velnpo.h'
  END INTERFACE


  IF (flag .AND. numct > 0) THEN      !generate array ID to assemble contact forces
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

  eqres = 0d0  !initializes forces
  !Add internal equivalente forces (resid)
  CALL ensve1(ndofn*npoin,ifpre(1,1),resid(1,1),eqres(1),npsdf(1),nesdf(1),ftsdf(1))
  !Add contact forces (FCONT)
  IF (numct > 0) CALL ensve1(ndimc*npoin,id(1,1),fcont(1,1),eqres(1),npsdf(1),nesdf(1),ftsdf(1))
  !add external equivalent forces (Force)
  IF (nload > 0) THEN
    n1 = nload + 1
    IF( str_type == 1 )THEN
      IF( flag )THEN
        lambda = lambd
        ! generation
        i = 0                                 !initializes number of loaded DOFs
        aux   = 0d0                           !initializes NORM
        DO ieq=1,neq                               !for each DOF
          IF( force(ieq,n1) /= 0d0 )THEN             !IF loaded
            i = i + 1                                !increase number of loaded DOFs
            aux   = aux   + ABS(force(ieq,n1))         !add to norm (includes Lambda factor)
          END IF
        END DO
        nf = i                                !number of loaded DOFs
        ALLOCATE(lind(nf),fn(nf))             !get memory for arrays
        aux   = lambda/aux                    ! lambda/|f|
        i = 0                                 !initializes counter
        DO ieq=1,neq                              !for each DOF
          IF( force(ieq,n1) /= 0d0 )THEN             !if loaded DOFs
            i = i + 1                                   !increase counter
            fn(i) =  SIGN(aux  ,force(ieq,n1))             !factor for DOF (includes sign)
            lind(i) = ieq                                  !keep associated equation
          END IF
        END DO
        dl2  =  DOT_PRODUCT(ddisp,ddisp)          !fix arc-length increment
        flag = .FALSE.                            !task done
      END IF
      lambd = 0d0                          !initializes factor
      DO i=1,nf                             ! for each loade DOF
        ieq  = lind(i)                        !associated equation
        lambd = lambd + eqres(ieq)*fn(i)    !average factor
      END DO
      factor = lambd/lambda
      DO i=1,nf                            !for each loaded DOF
        ieq  = lind(i)                        !associated equation
        eqres(ieq) = eqres(ieq) - force(ieq,n1)*factor  !include load
      END DO
    ELSE
      DO ieq=1,neq
        eqres(ieq) = eqres(ieq) - force(ieq,n1)
      END DO
    END IF
  END IF
  !Modify previous velocities due to damping
  IF(.NOT.chkin) THEN
 !$OMP DO
  DO ieq=1,neq
    veloc(ieq) = damp(ieq)*veloc(ieq)
  END DO
 !$OMP END DO
 END IF
  !Integrate equations
  ! compute standard accelerations
 !$OMP DO
  DO ieq = 1,neq
    acelr(ieq) = -eqres(ieq)*ymass(ieq)
  END DO
 !$OMP END DO
  IF( selective_mass_scaling )THEN    ! loop over fibres
    sms_faux = 0d0 !initializes non-equilibrated forces
    DO i=1,sms_nf  !for each fibre
      ! first loop to compute forces on the middle surface
      DO j=1,sms_nnf  !for each node on the fibre
        k = sms_fibre(j,i)           !node number
        IF( k == 0 )EXIT             !if all processed exit loop
        DO l=1,ndime                 !for each space dimension
          ieq = ifpre(l,k)           !global DOF
          IF( ieq < 1 )CYCLE         !if not an active dof go to next
          sms_faux(l,i) = sms_faux(l,i) - eqres(ieq)  !add force
        END DO
      END DO
      ! compute accelerations on the middle surface
      sms_faux(:,i) = sms_faux(:,i) / sms_mf(1,i)   !node accelerations (already scaled by (1-alpha)
      ! second loop to compute accelerations on the nodes
      DO j=1,sms_nnf  !for each node on the fibre
        k = sms_fibre(j,i)           !node number
        IF( k == 0 )EXIT             !if all processed exit loop
        DO l=1,ndime                 !for each space dimension
          ieq = ifpre(l,k)           !global DOF
          IF( ieq < 1 )CYCLE         !if not an active DOF go to next
          acelr(ieq) = acelr(ieq) * sms_mf(2,i) + sms_faux(l,i)  !modify acceleration
        END DO
      END DO
    END DO
  END IF

  !Integrate velocities and displacements
 !$OMP DO
  DO ieq = 1,neq
    veloc(ieq) = veloc(ieq) + acelr(ieq)
  END DO
 !$OMP END DO

  ! for load-strategy (arc-length method)
  IF( str_type == 1 )THEN   !compute incremental displacements
    !first proyect velocities on the hyper-plane normal to displacement increment
    aux  = -DOT_PRODUCT(ddisp,veloc)/dl2        !factor
   !$OMP DO
    DO ieq = 1,neq      !correct velocities
      veloc(ieq) = veloc(ieq) + ddisp(ieq)*aux  !constrained velocities
    END DO
   !$OMP END DO
    !second update displacements normally
   !$OMP DO
    DO ieq = 1,neq      !update displacements
      ddisp(ieq) = veloc(ieq) + ddisp(ieq)
    END DO
   !$OMP END DO
    !third constraint displacement using arc-length
    aux  =  SQRT(dl2/ DOT_PRODUCT(ddisp,ddisp))         !fix arc-length increment
    ddisp = ddisp*aux
  ELSE
   !$OMP DO
    DO ieq = 1,neq
      ddisp(ieq) = ddisp(ieq) + veloc(ieq)
    END DO
   !$OMP END DO
    factor = velfi
  END IF

  CALL coract(factor,velor(1:nvfix,nvelr+1),ddisp,coora,euler)


RETURN
END SUBROUTINE sxplit

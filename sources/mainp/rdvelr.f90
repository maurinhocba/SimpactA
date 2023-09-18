SUBROUTINE rdvelr(nvelr,ndime,ndofn,npoin,iwrit,label,ifpre,nvfix,lcvel,velor)

  !***  Applies (?) rigid body velocities
  USE curv_db, ONLY: getcun
  USE c_input
  USE rve_db
  USE kinc_db, ONLY : nn,nn1
  IMPLICIT NONE
  INTEGER (kind=4), INTENT(IN) :: ndime,ndofn,npoin,iwrit,label(:)
  INTEGER (kind=4), INTENT(IN OUT) :: ifpre(:,:),nvfix
  INTEGER (kind=4), INTENT(IN) :: nvelr
  INTEGER (kind=4), POINTER :: lcvel(:)
  REAL (kind=8), POINTER :: velor(:,:)

  INTEGER (kind=4) g,i,j,l,n,nrve,nposn,flag, chnode
  REAL (kind=8) xg(ndofn),xf(ndofn)
  REAL (kind=8), ALLOCATABLE :: auxil(:,:)
  TYPE (rve_set), POINTER :: rves
  TYPE (rve_nod), POINTER :: rven

  IF (nvelr > 0) THEN  !if prescribed velocities exist

    ALLOCATE( lcvel(nvelr+1) )                !space for associated curves
    ALLOCATE( auxil(0:ndofn*npoin,nvelr) )    !auxiliar space
    auxil = 0d0                               !initializes

    rves => headv                             !point to first set

    DO l=1,nvelr                   !loop on each set

      lcvel(L)=getcun (rves%lc)    ! get associated curve number for this set
      auxil(0,l)=rves%factor       !associated factor for this set
      flag = rves%dspflg
      IF(iwrit == 1) THEN
        WRITE(lures, "(//, &
        & 5X,'Curve scaling this velocity set ....',i10/ &
        & 5X,'Scaling factor .....................',e14.7,/)",ERR=9999)


        IF( flag == 0) THEN
          WRITE(lures, "(//,5X,'RIGID BODY VELOCITIES',//)",ERR=9999)
          IF (ndime==2 .AND. ndofn == 2) WRITE(lures, &
             "(6X,'NODE',5X,'X-VELO.',7X,'Y-VELO.'/)",ERR=9999)
          IF (ndime==2 .AND. ndofn == 3) WRITE(lures, &
             "(6X,'NODE',5X,'X-VELO.',7X,'Y-VELO.',7X,'A-VELO.'/)",ERR=9999)
          IF (ndime==3 .AND. ndofn == 3) WRITE(lures, &
             "(6X,'NODE',5X,'X-VELO.',7X,'Y-VELO.',7X,'Z-VELO.'/)",ERR=9999)
          IF (ndime==3 .AND. ndofn > 3) WRITE(lures, &
             "(6X,'NODE',5X,'X-VELO.',7X,'Y-VELO.',7X,'Z-VELO.', &
             &           7X,'A-VELO.',7X,'B-VELO.',7X,'C-VELO.'/)",ERR=9999)
        ELSE
          IF( flag == 1 )THEN
            WRITE(lures, "(//,5X,'RIGID BODY DISPLACEMENTS',//)",ERR=9999)
          ELSE
            WRITE(lures, "(//,5X,'RIGID BODY INCREMENTAL DISPLACEMENTS',//)",ERR=9999)
          END IF
          IF (ndime==2 .AND. ndofn == 2) WRITE(lures, &
             "(6X,'NODE',5X,'X-DISP.',7X,'Y-DISP.'/)",ERR=9999)
          IF (ndime==2 .AND. ndofn == 3) WRITE(lures, &
             "(6X,'NODE',5X,'X-DISP.',7X,'Y-DISP.',7X,'A-ROTA.'/)",ERR=9999)
          IF (ndime==3 .AND. ndofn == 3) WRITE(lures, &
             "(6X,'NODE',5X,'X-DISP.',7X,'Y-DISP.',7X,'Z-DISP.'/)",ERR=9999)
          IF (ndime==3 .AND. ndofn > 3) WRITE(lures, &
             "(6X,'NODE',5X,'X-DISP.',7X,'Y-DISP.',7X,'Z-DISP.', &
             &           7X,'A-ROTA.',7X,'B-ROTA.',7X,'C-ROTA.'/)",ERR=9999)
        END IF
      END IF

      nrve = rves%nrv         !number of nodes in this set
      rven => rves%head       !point to first node in the set
      DO j=1,nrve             !loop for each node in the set
        g = rven%node         !node label
        n = chnode(g)         !internal number

        xg(1:ndofn) = rven%v(1:ndofn)  !prescribed velocities components
        xf = xg                        !copy onto XF
        DO i=1,ndofn                   !for each DOF
          nposn = ifpre(i,n)           !restriction code
          SELECT CASE (nposn)
          CASE (:-nn1)                        !prescribed values
            auxil(-nposn-nn,l) = xg(i)           !assign value
          CASE (-nn:-1)                       ! slave DOF
            IF(xg(i)/=0d0)WRITE(lures,"(' WARNING: imposed veloc', &
            & e13.4,'on node',i9,' DOF',i3,' not possible')",ERR=9999) xg(i), label(n),i
            xf(i) = 0d0                       !not possible
          CASE (0,1)                          !active DOF or non existent
            IF(xg(i) /= 0d0)THEN                !if a non-zero value
              nvfix = nvfix+1                   !increase number of fixed values
              ifpre(i,n) = -nvfix-nn            !modify restriction code
              auxil(nvfix,l) = xg(i)            !assign value
            END IF
          END SELECT
        END DO
        ! echo effectively assigned values
        IF(iwrit == 1) WRITE(lures,"(i10,6e14.5)",ERR=9999) label(n), xf(1:ndofn)
        rven => rven%next        !point to next node
      END DO
      rves => rves%next       !point to next set
    END DO

    ALLOCATE( velor(nvfix+1,nvelr+1) )     !reserve space for definitive array
    velor(1:nvfix,1:nvelr) = auxil(1:nvfix,1:nvelr)  !transfer velocity data
    velor(nvfix+1,1:nvelr) = auxil(0,1:nvelr)        !transfer factors
    velor(1:nvfix,nvelr+1) = 0.d0  !initializes compounded prescribed velocities
    DEALLOCATE ( auxil )           !release auxiliar array

  ELSE  !no prescribed velocities

    ALLOCATE( velor(MAX(nvfix+1,1),1) ) !to avoid null pointers only
    ALLOCATE( lcvel(nvelr+1) )
    velor = 0d0

  END IF

RETURN
 9999 CALL runen2('')
END SUBROUTINE rdvelr

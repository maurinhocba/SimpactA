SUBROUTINE rdvel0 (iwrit,ndime,ndofn,actio,nvelr)

  !This routine reads prescribed velocities

  ! the convention: the last read overrides the previous data
  ! prescribed velocity will be applied to the nodes previously
  ! fixed o released

  USE param_db,ONLY: mnam
  USE c_input
  USE rve_db
  USE curv_db, ONLY : del_cur
  IMPLICIT NONE
  INTERFACE
     INCLUDE 'rdvels.h'
  END INTERFACE
  CHARACTER(len=*),INTENT(INOUT):: actio
  INTEGER (kind=4) :: iwrit,ndime,ndofn,nvelr
  INTEGER :: nrve, flag

  ! Local
  TYPE (rve_set), POINTER :: rves,anter,posic
  CHARACTER (len=mnam) :: lbl
  LOGICAL :: found

  IF (iwrit == 1) WRITE(lures, "(//,' P R E S C R I B E D   V E L O C I T I E S ',/)",ERR=9999)

  IF (nvelr > 0) THEN              !if previous prescribed values
    CALL listen('RDVEL0')          !read a card
    IF (exists('DELETE')) THEN     !Key word DELETE is present
      IF (iwrit == 1) WRITE(lures,"(' Deleting prescribed velocity sets ',/)",ERR=9999)
      IF (TRIM(actio) == 'NSTRA0') actio='NSTRA1'

      DO                           ! loop over the sets to delete
        CALL listen('RDVEL0')      ! read a card
        IF (exists('ENDDEL') ) EXIT     !kew-word END_DELETE found => exit loop
        lbl = get_name('VELSET',stype='VELS',        &   !assoc curve
              texts='!PRESCRIBED VELOCITY SET ..........')
        CALL srch_rves (headv, anter, rves, lbl, found)  !search for the curve
        IF (.NOT.found) THEN                    !If not found error in data input
          WRITE(lures, "(' Warning! Prescribed Velocity Set namee', &
     &              A,' does not exist')",ERR=9999) TRIM(lbl)
        ELSE
          CALL del_rves (headv, tailv, anter, rves)   !delete set assoc to LC
          CALL del_cur (lbl)          !delete curve data associated LC
          nvelr = nvelr - 1           !correct counter
        END IF
      END DO

    ELSE              !nothing to delete
      backs = .TRUE.                       !one line back
    ENDIF

    CALL listen('RDVEL0')          !read a card
    IF (exists('MODIFY')) THEN     !Key word MODIFY is present
      IF (iwrit == 1) WRITE(lures,"(' Modifying prescribed velocity curves ',/)",ERR=9999)

      DO                           ! loop over the sets to modify
        CALL listen('RDVEL0')      ! read a card
        IF (exists('ENDMOD') ) EXIT     !kew-word END_MODIFY found => exit loop
        lbl = get_name('VELSET',stype='VELS',        &   !assoc curve
              texts='!PRESCRIBED VELOCITY SET ..........')
        CALL srch_rves (headv, anter, rves, lbl, found)  !search for the curve
        IF (.NOT.found) THEN                    !If not found error in data input
          WRITE(lures, "(' ERROR!!! Prescribed Velocity Set name', &
     &              A,' does not exist')",ERR=9999) TRIM(lbl)
          CALL runend('Prescribed velocity set not found')
        ELSE
          CALL del_cur (lbl)          !delete curve data associated LC
          ! copied from below, should work
          rves%factor=getrea('FACTOR',1d0, ' Participation factor of this set .')
          flag = 0 !standard velocity definition
          IF( exists('DISPLA')) flag = 1                 !Displacement definition
          IF( flag == 1 .AND. exists('INCREM')) flag = 2
          rves%dspflg = flag                             !displacement flag
          IF( flag > 0 .AND. exists('SIMPLE')) flag = flag + 2
          IF( flag > 0 .AND. exists('LINACC')) flag = flag + 4
          CALL rdcurv('VELOC',rves%lc)
          ! ------------
        END IF
      END DO

    ELSE              !nothing to modify
      backs = .TRUE.                       !one line back
    END IF

  END IF

  IF (iwrit == 1) WRITE (lures, "(//)",ERR=9999)

  DO
    ! loop over velocity sets - reading
    CALL listen('RDVEL0')              !read a card
    IF (.NOT.exists('VELSET')) THEN    ! if key-word VEL_SET not found
      backs = .TRUE.                       !one line back
      EXIT                             ! exit loop
    END IF
    ALLOCATE (rves)                    !get memory for a list of nodes (SET)
    rves%lc = get_name('VELSET',stype='VELS',        &   !assoc curve
              texts='!PRESCRIBED VELOCITY SET ...........')
    rves%factor=getrea('FACTOR',1d0, ' Participation factor of this set .')
    flag = 0 !standard velocity definition
    IF( exists('DISPLA')) flag = 1                 !Displacement definition
    IF( flag == 1 .AND. exists('INCREM')) flag = 2
    rves%dspflg = flag                             !displacement flag
    IF( flag > 0 .AND. exists('SIMPLE')) flag = flag + 2
    IF( flag > 0 .AND. exists('LINACC')) flag = flag + 4
    !check if associated label LC already used (if TRUE stop)

    CALL srch_rves (headv, anter, posic, rves%lc, found)
    IF (found) CALL runend ('RDVEL0: Set using this name already')
    CALL rdcurv('VELOC',rves%lc)
    IF( flag /= 0) CALL ch_d2v(flag) !'VELOC',rves%lc,
    IF(iwrit == 1) THEN    !echo title
      IF( flag == 0 )THEN
        WRITE(lures, "(//,5X,'RIGID BODY VELOCITIES',//)",ERR=9999)
        IF (ndime==2 .AND. ndofn == 2) WRITE(lures,  &
           "(6X,'NODE',5X,'X-VELO.',7X,'Y-VELO.'/)",ERR=9999)
        IF (ndime==2 .AND. ndofn == 3) WRITE(lures,  &
           "(6X,'NODE',5X,'X-VELO.',7X,'Y-VELO.',7X,'A-VELO.'/)",ERR=9999)
        IF (ndime==3 .AND. ndofn == 3) WRITE(lures,  &
           "(6X,'NODE',5X,'X-VELO.',7X,'Y-VELO.',7X,'Z-VELO.'/)",ERR=9999)
        IF (ndime==3 .AND. ndofn > 3) WRITE(lures,  &
           "(6X,'NODE',5X,'X-VELO.',7X,'Y-VELO.',7X,'Z-VELO.',  &
           &           7X,'A-VELO.',7X,'B-VELO.',7X,'C-VELO.'/)",ERR=9999)
      ELSE
        IF( flag == 1 )THEN
          WRITE(lures, "(//,5X,'RIGID BODY DISPLACEMENT',//)",ERR=9999)
        ELSE
          WRITE(lures, "(//,5X,'RIGID BODY INCREMENTAL DISPLACEMENT',//)",ERR=9999)
        END IF
        IF (ndime==2 .AND. ndofn == 2) WRITE(lures,  &
           "(6X,'NODE',5X,'X-DISP.',7X,'Y-DISP.'/)",ERR=9999)
        IF (ndime==2 .AND. ndofn == 3) WRITE(lures,  &
           "(6X,'NODE',5X,'X-DISP.',7X,'Y-DISP.',7X,'A-ROTA.'/)",ERR=9999)
        IF (ndime==3 .AND. ndofn == 3) WRITE(lures,  &
           "(6X,'NODE',5X,'X-DISP.',7X,'Y-DISP.',7X,'Z-DISP.'/)",ERR=9999)
        IF (ndime==3 .AND. ndofn > 3) WRITE(lures,  &
           "(6X,'NODE',5X,'X-DISP.',7X,'Y-DISP.',7X,'Z-DISP.',  &
           &           7X,'A-ROTA.',7X,'B-ROTA.',7X,'C-ROTA.'/)",ERR=9999)
      END IF
    END IF
    IF (TRIM(actio) == 'NSTRA0') actio='NSTRA1'

    CALL rdvels (iwrit, rves%head, rves%tail, nrve, ndofn)

    rves%nrv=nrve                       !keep number of nodes in the set
    CALL add_rves (rves, headv, tailv)  !add set of nodes to the end of the list
    nvelr = nvelr + 1                   !increment number of sets

  END DO  ! loop for sets of prescribed velocities
  CALL listen('RDVEL0')  ! read a line (END_PRESCRIBED line expected)
  IF (.NOT. exists('ENDPRE'))CALL runend ('RDVEL0:end_prescribed_veloci expc.')

  RETURN
 9999 CALL runen2('')
END SUBROUTINE rdvel0

SUBROUTINE rdtemp (iwrit,actio,flag)

  !***  READ fixed and prescribed temperatures

  USE param_db,ONLY: mnam
  USE ift_db       !temperature prescribed data_base
  USE c_input
  USE nsets_db
  USE curv_db, ONLY : del_cur

  IMPLICIT NONE
  CHARACTER(len=*) :: actio
  INTEGER (kind=4) :: iwrit, & !for echo into .RSN
                      flag     ! 0:fixed or released values  1:curve dependent values

  ! Local
  LOGICAL set, found
  CHARACTER(len=mnam) :: sname,lbl
  TYPE (ift_nod), POINTER :: ift, anter, posic
  INTEGER :: i,j,k,n,ki,kf,icod,nnods,nod,nrve,nrpt
  INTEGER (kind=4), ALLOCATABLE :: nods(:)
  REAL (kind=8) :: v
  TYPE (rpt_set), POINTER :: rves,rpts,sant,spos
  TYPE (rpt_nod), POINTER :: rven
  TYPE (nset), POINTER :: ns

  INTERFACE
     INCLUDE 'rdpret.h'
  END INTERFACE

 !------------------------------------------------------------

 IF( flag == 1 )THEN    !this part reads fixed values and released DOFs

   IF (iwrit == 1) WRITE(lures,"(/,' Fixed Temperature Conditions',/)",ERR=9999)

   !Initialize empty list

   IF ( nift==0 .OR. .NOT.exists('ADD   ')) THEN
     CALL ini_ift(head,tail)
     nift=0                       !initializes for new boundary conditions
   ELSE IF (exists('ADD   ')) THEN
     IF (iwrit==1)WRITE(lures,"(/,' Kept Temperature Conditions from the previous strategy.',/)",ERR=9999)
   END IF

   IF (TRIM(actio) == 'NSTRA0') actio = 'NSTRA1'

   IF (iwrit == 1) THEN
     IF(ndoft == 1) WRITE(lures,"(/,'   Node    TN   value')",ERR=9999)
     IF(ndoft == 2) WRITE(lures,"(/,'   Node    TI TS   val TI  Val TS')",ERR=9999)
     IF(ndoft == 3) WRITE(lures,"(/,'   Node    TN TI TS  val TN  val TI  Val TS')",ERR=9999)
   END IF

   !Loop to read data and add them to the list
   DO
     CALL listen('RDTEMP')          !read a line
     IF (exists('ENDFIX')) EXIT     !end of data =>exit

     IF (exists('SET   ',n))THEN                    !key-word SET found
       sname = get_name(posin=n,stype='NSET')       !set name
       CALL nsdb_search (sname, found, ns)          !search in list of sets
       IF (found) THEN                      !set exists, position is NS
         nnods = get_length (ns)            !number of nodes in the set
         ALLOCATE (nods(nnods))             !get memory
         CALL ns_head (ns)                  !go to top of list
         DO i =1, nnods                     !for each node in the list
           nods(i) = get_label(ns)          !get node label
           CALL ns_next (ns)                !go to next node
         END DO
       ELSE                                 !else SET name not found => error
         WRITE (lures, '(" Set ",a,"  not found")',ERR=9999) sname(1:LEN_TRIM(sname))
         CALL runend('Intime: Set not found ')
       END IF
       ki = 1                               !first node
       kf = nnods                           !last node
       set = .TRUE.                         !set flag to TRUE
     ELSE                             !read fixities for a single node
       ki  = INT(param(1))                  !node label
       kf = ki                              !same as last
       set = .FALSE.                        !set flat to FALSE
     END IF
                                      !store fixities in list
     DO n = ki,kf                           !for each node in the list
       ALLOCATE (ift)                       !get memory
       IF (set) THEN                        !for a whole set
         nod = nods(n)                      !node label
       ELSE
         nod = n                            !only node
       END IF
       ift%ifix(1) = nod                    !store node label
       ift%ifix(2:ndoft+1)= INT (param(2:ndoft+1))   !store fixities
       ift%val(1:ndoft)   = param(ndoft+2:1+2*ndoft) !store values
       IF (iwrit == 1) THEN
         IF( ndoft == 1 ) WRITE(lures,"(i7,2x, i3, f9.2)",ERR=9999)ift%ifix(1:2),ift%val(  1)  !echo
         IF( ndoft == 2 ) WRITE(lures,"(i7,2x,2i3,2f9.2)",ERR=9999)ift%ifix(1:3),ift%val(1:2)  !echo
         IF( ndoft == 3 ) WRITE(lures,"(i7,2x,3i3,3f9.2)",ERR=9999)ift%ifix(1:4),ift%val(1:3)  !echo
       END IF

       !    search in all the prescribed values to change conditions

       rves => headv      !point to the first set of prescribed values

       DO i = 1,npret     !for each prescribed temperature set
         ! the convention: the last read overrides the previous data
         ! if prescribed temperature has been read previously it will
         ! be cancelled at fixed or released degrees of freedom
         nrve = rves%nrv       !number of prescribed values in the set
         rven => rves%head     !point to first node in the set
         DO j=1,nrve                        !for each node in the set
           IF (nod == rven%node) THEN       !compare labels
             DO k = 1,ndoft                    !for each DOF
               v = rven%v(k)                   !prescribed temperature in that DOF
               IF (v /= 0d0)THEN               !if original temperature /= 0
                 icod= INT (param(k+1))        !fixity code
                 rven%v(k) = 0d0               !initializes to 0
                 IF (icod == 0) THEN           !if DOF is free now
                   WRITE (lures,"(' Warning. Node:',i8,', DOF:',i2,' released.'/&
                      & ' Previously prescribed temperature cancelled.')",ERR=9999) nod,k
                 ELSE !icod == 1 (prescribed value)
                   WRITE (lures,"(' Warning. Node:',i8,', DOF:',i2,' fixed.'/ &
                      & ' Previously prescribed temperature cancelled.')",ERR=9999) nod,k
                 END IF
               END IF
             END DO
             EXIT      !node label was found => exit this list of nodes
           END IF
           rven => rven%next   !point to next node in the list
         END DO        !list of nodes of a set of prescribed values
         rves => rves%next     !point to next list of nodes
       END DO   ! list of sets

       ! before adding check if the node has not been fixed previously
       CALL srch_ift (head, anter, posic, ift%ifix(1), found)
       IF (found) THEN
         WRITE (lures, "(' Warning. At node',i8,' previously prescribed ', &
                       & 'boundary conditions overwritten.')",ERR=9999)  ift%ifix(1)
         CALL del_ift (head, tail, anter, posic)
       ELSE
         nift=nift+1       !increase counter of fixed nodes
       END IF
       CALL add_ift( ift, head, tail )    !add node to the list
     END DO   !list of nodes (if set)
     IF (set) DEALLOCATE (nods)  !deallocate temporary array
   END DO  !fixed nodes, end of data of boundaries found

 !------------------------------------------------------------

 ELSE  !flag == 2 This part reads prescribed temperatures

  IF (iwrit == 1) WRITE(lures, "(' PRESCRIBED TEMPERATURES IN TIME ',/)",ERR=9999)

  IF (npret > 0) THEN              !if previous prescribed values
    CALL listen('RDTEMP')          !read a card
    IF (exists('DELETE')) THEN     !Key word DELETE is present
      IF (iwrit == 1) WRITE(lures,"(' Deleting prescribed temperature sets ',/)",ERR=9999)
      IF (TRIM(actio) == 'NSTRA0') actio = 'NSTRA1'

      DO                           ! loop over the sets to delete
        CALL listen('RDTEMP')      ! read a card
        IF (exists('ENDDEL') ) EXIT     !kew-word END_DELETE found => exit loop
        lbl = get_name('TEMSET',stype='CURV',        &   !assoc curpt
              texts='!PRESCRIBED TEMPERATURE SET .......')
        CALL srch_rpts (headv, sant, rpts, lbl, found)  !search for the curpt
        IF (.NOT.found) THEN                    !If not found error in data input
          WRITE(lures, "(' Warning! Temperature set using curpt', &
     &              a,' does not exist')",ERR=9999)TRIM(lbl)
        ELSE
          CALL del_rpts (headv, tailv, sant, rpts)   !delete set assoc to LC
          CALL del_cur (lbl)                         !delete curve data associated LC
          npret = npret - 1                          !correct counter
        END IF
      END DO

    ELSE              !nothing to delete
      backs = .TRUE.                       !one line back
    ENDIF

  END IF

  IF (iwrit == 1) WRITE (lures, "(//)",ERR=9999)

  DO
    ! loop over temperature sets - reading
    CALL listen('RDTEMP')              !read a card
    IF (.NOT.exists('TEMSET')) THEN    ! if key-word TEM_SET not found
      backs = .TRUE.                       !one line back
      EXIT                             ! exit loop
    END IF
    ALLOCATE (rpts)                    !get memory for a list of nodes (SET)
    rpts%lc = get_name('TEMSET',stype='CURV',        &   !assoc curpt
              texts='!PRESCRIBED TEMPERATURE SET ...')
    rpts%factor=getrea('FACTOR',1d0, ' Participation factor of this set .')

    !check if associated function LC already used (if TRUE stop)

    CALL srch_rpts (headv, sant, spos, rpts%lc, found)
    IF (found) CALL runend ('RDTEMP: Set using this name already')
    CALL rdcurv('VELOC',rpts%lc)

    IF(iwrit == 1) THEN    !echo title
      WRITE(lures, "(//,5X,' PRESCRIBED TEMPERATURES',//)",ERR=9999)
      IF ( ndoft == 1) WRITE(lures,  &
         "(6X,'NODE',5X,'TEMP.'/)",ERR=9999)
      IF ( ndoft == 2) WRITE(lures,  &
         "(6X,'NODE',5X,'TEMP-I',9X,'TEMP-S.',/)",ERR=9999)
      IF ( ndoft == 3) WRITE(lures,  &
         "(6X,'NODE',5X,'TEMP-N',8X,'TEMP-I',8X,'TEMP-S',/)",ERR=9999)
    END IF
    IF (TRIM(actio) == 'NSTRA0') actio = 'NSTRA1'

    CALL rdpret (iwrit, rpts%head, rpts%tail,nrpt) !read prescribed temp sets

    rpts%nrv=nrpt                       !keep number of nodes in the set
    CALL add_rpts (rpts, headv, tailv)  !add set of nodes to the end of the list
    npret = npret + 1                   !increment number of sets

  END DO  ! loop for sets of prescribed temperatures
  CALL listen('RDTEMP')  ! read a line (END_PRESCRIBED line expected)
  IF (.NOT. exists('ENDPRE'))CALL runend ('RDTEMP:end_prescribed_temper expc.')
  RETURN
  ! the convention: the last read overrides the previous data
  ! prescribed temperature will be applied to the nodes previously
  ! fixed o released

 END IF

 RETURN
 9999 CALL runen2('')
END SUBROUTINE rdtemp

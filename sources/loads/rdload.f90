 SUBROUTINE rdload (iwrit,ndime,ndofn, nload)
 !This routine reads load sets (nodal loads, edge and surface loads)
 USE param_db,ONLY: mnam
 USE c_input
 USE loa_db !,ONLY: loa_nod, loa_set, ifloa, headls, taills,        & !Module type and variables
            !       add_loa, new_loa,                               & !Module 'loa_nod' routines
            !       add_loas, del_loas, new_loas, srch_loas           !Module routines
 USE esets_db, ONLY : delete_name, add_name
 USE curv_db, ONLY : del_cur
 USE npo_db, ONLY : coora
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4):: iwrit, ndime, ndofn, nload
   !--- Local variables
   CHARACTER(len=mnam):: lbl
   INTEGER(kind=4):: nedge, nsurf, numfl, iplod
   REAL(kind=8):: dummy
   LOGICAL:: found
   TYPE(loa_nod),POINTER:: loa
   TYPE(loa_set),POINTER:: loas, anter

   INTERFACE
     INCLUDE 'dedge0.h'
     INCLUDE 'dsurf0.h'
     INCLUDE 'rdfoll.h'
     INCLUDE 'getnod.h'
   END INTERFACE
   !---------------------------------------------------------------

   !Allows deletion of previous load sets

   IF (nload > 0) THEN       !if loads of previous strategies exist
     CALL listen('RDLOAD')   !read a card
     IF (exists('DELETE')) THEN     !key-word DELETE_LOAD_SET
       ! delete load sets (give a list of load set labels)
       IF (iwrit == 1) WRITE(lures,"(' Deleting load sets ',/)",ERR=9999) !echo

       DO     !loop over the set to delete
         CALL listen('RDLOAD')        ! read a card
         IF (exists('ENDDEL') ) EXIT  ! key word END_DELETE found, exit loop
         lbl = get_name('LOADSE', &
                        texts='!LOAD SET LABEL TO DELETE .........',stype='LOAD') !set to delete
         CALL srch_loas(headls, anter, loas, lbl, found)  !search set
         IF (.NOT.found) THEN                  !set not found, ERROR
           WRITE(lures,"(' Load set ',a,' does not exist')",ERR=9999) TRIM(lbl)
           !CALL runend ('RDLOAD:Specified set does not exist')
         ELSE                                  !set found, delete it
           nload = nload - 1               !updates number of remaining sets
           IF (loas%numfl > 0 .OR.       &  !if follower loads in the set
               ASSOCIATED(loas%hydfl)) ifloa = ifloa - 1 !updates counter
           CALL del_loas(headls, taills, anter, loas)  !delete load set
           CALL del_cur(lbl)          !delete curve data associated LC
           CALL delete_name(lbl,7)
         END IF
       END DO

     ELSE               ! if no set to delete
       backs = .TRUE. !  one card back
     END IF
   END IF

    CALL listen('RDLOAD')          !read a card
    IF (exists('MODIFY')) THEN     !Key word MODIFY is present
      IF (iwrit == 1) WRITE(lures,"(' Modifying applied loads curves ',/)",ERR=9999)

      DO                           ! loop over the sets to modify
        CALL listen('RDLOAD')      ! read a card
        IF (exists('ENDMOD') ) EXIT     !kew-word END_MODIFY found => exit loop
        lbl = get_name('LOADSE',stype='LOAD',        &   !assoc curve
              texts='!APPLIED LOAD SET .................')
        CALL srch_loas(headls, anter, loas, lbl, found)  !search set
        IF (.NOT.found) THEN                    !If not found error in data input
          WRITE(lures, "(' ERROR! Applied LOADS Set name', &
     &              A,' does not exist')",ERR=9999) TRIM(lbl)
        ELSE
          CALL del_cur (lbl)          !delete curve data associated LC
          loas%factor = getrea('FACTOR',1d0,' SCALE FACTOR OF FOR THIS SET .....')
          CALL rdcurv('FORCE',loas%lbl)
        END IF
      END DO

    ELSE              !nothing to modify
      backs = .TRUE.                       !one line back
    END IF
   !Read  new load sets
   DO    !loop over reference load sets
     CALL listen('RDLOAD')              !read a card
     IF (.NOT.exists('LOADSE')) THEN    !NO new set to read, exit loop
       backs = .TRUE.                   !one card back
       EXIT                             !exit loop
     END IF

     WRITE (lures,"(//)",ERR=9999)

     lbl = get_name('LOADSE',texts='!LOAD SEt label ...................',stype='LOAD')
     CALL srch_loas(headls, anter, loas, lbl, found)  !search set
     IF (found) CALL runend ('RDLOAD: LOAD Set using this name already')
     CALL new_loas(loas)                !get memory for the new set
     loas%lbl = lbl
     loas%factor = getrea('FACTOR',1d0,' SCALE FACTOR OF FOR THIS SET .....')
     CALL add_name(loas%lbl,7)

     CALL listen('RDLOAD')                !TITLE
     !header
     IF(iwrit == 1) WRITE (lures,"(/5X,'======== Reference Load Set No. ',A,  &
                    '  ========',//,5X,'Load Case  -',A)",ERR=9999) TRIM(loas%lbl), TRIM(card)

     CALL rdcurv('FORCE',loas%lbl)

     CALL listen('RDLOAD')                !read gravity card
     IF (.NOT.exists('GRAVIT')) backs = .TRUE.

     WRITE (lures,"(/)",ERR=9999)
     igrav=getint('IGRAV ',0,' Gravitational load code ..........')
     loas%igrav = igrav                   !gravity flag

     IF(igrav /= 0) THEN
       !*** READ gravity direction and gravitational constant
       loas%gravy    = getrea('GCON  ',9.81d0,' GRAVitational Constant ...........')
       loas%gvect(1) = getrea('GVECTX',0d0, ' Grav. VECTor X component..........')
       IF(ndime == 2)THEN
         loas%gvect(2) = getrea('GVECTY',-1d0,' Grav. VECTor Y component..........')
       ELSE !IF(ndime == 3)
         loas%gvect(2) = getrea('GVECTY',0d0, ' Grav. VECTor Y component..........')
         loas%gvect(3) = getrea('GVECTZ',-1d0,' Grav. VECTor Z component..........')
       END IF
       CALL vecuni(ndime,loas%gvect,dummy)
       IF( dummy == 0d0) CALL runend('RDLOAD:GRAVITATIONAL DIREC. IS NULL')
     END IF

     !*** READ nodal point loads
     CALL listen('RDLOAD')
     IF (.NOT.exists('POINTL')) THEN       !no point loads
       backs = .TRUE.           !one card back
     ELSE

       IF (iwrit == 1)THEN
         WRITE(lures, "(/' Nodes with concentrated forces ')",ERR=9999)
         IF(ndime == 2) THEN
           WRITE(lures,"(/6x,'node',4x,'px',8x,'py',8x,'mz')",ERR=9999)
         ELSE
           WRITE(lures,"(/6x,'node',4x,'px',8x,'py',8x,'pz',8x,'mx',8x,'my',8x,'mz')",ERR=9999)
         END IF
       END IF

       iplod = 0                        !initializes counter
       DO
         !Loop to read data and add them to the list
         CALL listen('RDLOAD')         !read a card
         IF (exists('ENDPOI')) EXIT    !key-word END_POINT found, exit loop

         CALL new_loa(loa)             !reserve space
         loa%node = INT(param(1))      !nodal label
         loa%forc(1:ndofn) = param(2:ndofn+1) !force components
         iplod = iplod + 1             !increase counter
         IF(iwrit == 1) WRITE(lures,"(5x,i5,6g10.3)",ERR=9999) loa%node, loa%forc(1:ndofn)    !ECHO
         CALL add_loa(loa, loas%headn, loas%tailn)     !add to list
       END DO

       loas%iplod = iplod              !store counter
     END IF

     !*** READ line loads
     CALL dedge0(nedge,ndime,iwrit,loas%heade,loas%taile)
     loas%nedge = nedge                !store counter

     !*** READ surface loads
     CALL dsurf0(nsurf,ndime,iwrit,loas%heads,loas%tails)
     loas%nsurf = nsurf                !store counter

     !*** READ follower loads
     CALL rdfoll(numfl,ndime,iwrit,                   &
                 loas%fltype,loas%flpar,loas%headf,loas%tailf,   &
                 loas%hydfl,loas%lbl,loas%fluid)
     loas%numfl = numfl
     IF(loas%fluid)THEN
       CALL runend('this version does allow fluid coupling')
     END IF

     CALL listen('RDLOAD')             !read last card in the set
     IF (.NOT.exists('ENDSET')) CALL runend('RDLOAD: END_SET CARD EXPECTED') !key-word END_SET expected

     nload = nload + 1          !increase number of load sets
     IF (loas%numfl > 0 .OR. ASSOCIATED(loas%hydfl)) ifloa=ifloa+1  !increase number of sets with follower loads
     CALL add_loas(loas, headls, taills)   ! add set to the list

   END DO  ! over load sets

   IF (ifloa > 0) CALL wrtfl0

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rdload

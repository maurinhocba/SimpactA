 SUBROUTINE rdheat (iwrit,ndoft,nheat,actio)

 !This routine reads heat sets (nodal heats, edge and surface heats)

 USE ctrl_db, ONLY: ntype,ndime
 USE c_input
 USE heat_db
 USE esets_db, ONLY : delete_name, add_name
 USE curv_db, ONLY : del_cur
 IMPLICIT NONE
 INTERFACE
    INCLUDE 'rdhedg.h'
    INCLUDE 'rdhsur.h'
 END INTERFACE
 CHARACTER(len=6 ) :: actio
 INTEGER (kind=4) :: iwrit,ndoft,nheat

 ! Local
 CHARACTER(len=mnam):: lbl
 INTEGER :: nedge,nsurf,ipheat
 LOGICAL :: found
 TYPE (heat_nod), POINTER :: hean
 TYPE (heat_set), POINTER :: heat,anter
 !---------------------------------------------------------------

 !     allows deletion of previous heat sets

 IF (nheat > 0) THEN       !if heats of previous strategies exist
   CALL listen('RDHEAT')   !read a card
   IF (exists('DELETE')) THEN     !key-word DELETE_HEAT_SET
     ! delete heat sets (give a list of heat set labels)
     IF (iwrit == 1) WRITE(lures,"(' Deleting heat sets ',/)",ERR=9999) !echo
     IF (actio == 'NSTRA0') actio = 'NSTRA1'    !??

     DO     !loop over the set to delete
       CALL listen('RDHEAT')        ! read a card
       IF (exists('ENDDEL') ) EXIT  ! key word END_DELETE found, exit loop
       lbl = get_name('HEATSE',texts='!HEAT SET LABEL TO DELETE .........', &
                       stype='HEAT') !set to delete
       CALL srch_heats (headhs, anter, heat, lbl, found)  !search set
       IF (.NOT.found) THEN                  !set not found, ERROR
         WRITE(lures, "(' Heat set ',a,' does not exist')",ERR=9999)lbl(1:LEN_TRIM(lbl))
         CALL runend ('RDHEAT:Specified set does not exist')
       ELSE                                  !set found, delete it
         nheat = nheat - 1               !updates number of remaining sets
         CALL del_heats (headhs, tailhs, anter, heat)  !delete heat set
         CALL del_cur(lbl)          !delete curve data associated LC
         CALL delete_name(lbl,8)
       END IF
     END DO

   ELSE               ! if no set to delete
     backs = .TRUE.   ! one card back
   END IF
 END IF

 heat => headhs       !point to first heat set

 !     Read NEW heat sets in present strategy

 DO    !loop over reference heat sets
   CALL listen('RDHEAT')              !read a card
   IF (.NOT.exists('HEATSE')) THEN    !NO new set to read, exit loop
      backs = .TRUE.                  !one card back
     EXIT                             !exit loop
   END IF

   IF (iwrit == 1) WRITE (lures, &   !Header for HEAT data
       "(/,'  H E A T   R E A D   I N   T H E   ', &
       &   'P R E S E N T   S T R A T E G Y ',/)",ERR=9999)

   WRITE (lures,"(//)",ERR=9999)

   CALL alloc_heats (heat)             !get memory for the new set

   heat%lbl = get_name('HEATSE',texts='!HEAT SEt label ...................',stype='HEAT')
   !additional factor
   heat%factor = getrea('FACTOR',1d0,' SCALE FACTOR FOR THIS SET ........')
   CALL add_name(heat%lbl,8)

   CALL listen('RDHEAT')                !TITLE
   IF(iwrit == 1) WRITE (lures,"(/5X,'======== Reference Heat Set No. ',a,   &
                    &     '  ========',//,5X,'Heat Case  -',a)",ERR=9999)             &
                         TRIM(heat%lbl),TRIM(card)

   CALL rdcurv('FORCE',heat%lbl)
   WRITE (lures,"(/)",ERR=9999)

   !*** READ nodal point heats

   CALL listen('RDHEAT')
   IF (.NOT.exists('POINTS')) THEN       !no point heats
     heat%ipsour = 0
     backs = .TRUE.   !one card back
   ELSE

     IF(iwrit == 1)THEN
       WRITE(lures, "(/' Nodes with concentrated source ')",ERR=9999)
       IF(ndoft == 1) THEN
         WRITE(lures,"(/6x,'node',4x,'f')",ERR=9999)
       ELSE IF(ndoft == 2) THEN
         WRITE(lures,"(/6x,'node',4x,'Fi',8x,'Fs')",ERR=9999)
       ELSE IF(ndoft == 3) THEN
         WRITE(lures,"(/6x,'node',4x,'Fn',8x,'Fi',8x,'Fs')",ERR=9999)
       END IF
     END IF

     ipheat = 0                             !initializes counter
     CALL ini_heat(heat%headn, heat%tailn ) !initializes pointers
     DO
       !Loop to read data and add them to the list
       CALL listen('RDHEAT')         !read a card
       IF (exists('ENDPOI')) EXIT    !key-word END_POINT found, exit loop
       ALLOCATE (heat)               !reserve space
       hean%node = INT(param(1))     !nodal label
       hean%source(1:ndoft)= param(2:ndoft+1) !source components
       ipheat=ipheat+1               !increase counter
       IF(iwrit == 1) WRITE(lures,"(5x,i5,6g10.3)",ERR=9999)hean%node, hean%source(1:ndoft) !ECHO
       CALL add_heat( hean, heat%headn, heat%tailn )     !add to list
     END DO

     heat%ipsour = ipheat              !store counter
   END IF

   !*** READ line heats

   CALL rdhedg(nedge,iwrit,heat%heade,heat%taile)
   heat%nedge = nedge                !store counter
   heat%ntype = ntype                !store NTYPE (2-D problems)

   !*** READ surface heats

   CALL rdhsur(nsurf,iwrit,heat%heads,heat%tails)
   heat%nsurf = nsurf                !store counter

   CALL listen('RDHEAT')             !read last card in the set
   !key-word END_SET expected
   IF (.NOT.exists('ENDSET'))CALL runend('RDHEAT: END_SET CARD EXPECTED  ')

   nheat = nheat + 1          !increase number of heat sets
   CALL add_heats( heat, headhs, tailhs )   ! add set to the list

 END DO  ! over heat sets

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rdheat

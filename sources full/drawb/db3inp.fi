SUBROUTINE db3inp ( ndold)

  !                3d drawbead input routine

  USE param_db,ONLY: mnam
  USE ctrl_db, ONLY: nconp, npoin
  USE c_input
  USE npo_db, ONLY : coora,label
  USE esets_db, ONLY : delete_name
  USE surf_db
  IMPLICIT NONE
  ! dummy arguments
  INTEGER (KIND=4) :: ndold  !number of existing DB
  ! local variables
  LOGICAL :: found    !flag for element set search
  CHARACTER (LEN=mnam) :: name    !sheet name
  INTEGER (KIND=4), ALLOCATABLE :: iwork(:) !auxiliar array
  INTEGER (KIND=4):: nrtm,    & !number of elements to be checked
                     nmn,     & !
                     mst,     & !nothing in this routine
                     ndbpoin, & !number of points in the db line
                     nnode,   & !number of nodes in sheet conns
                     ncrosi,  & !nothing in this routine
                     nelem      !number of elements in the sheet

  INTEGER (KIND=4), POINTER :: icode(:) !
  TYPE (db3pair_db), POINTER :: anter, pair  !auxiliar pointers
  TYPE (cont_srf), POINTER :: surface


  CALL listen ('dbea3d')     !read keyword line
  IF ( exists('DRAWBE') ) THEN       !if drawbead information exists

    !  redefinition of drawbeads
    CALL listen ('dbea3d')   ! read puncht
    IF ( exists('PUNCHT')) THEN
      puncht = getrea ('PUNCHT',1.d10,'!Maximum punch travel .............')
    ELSE
      backs = .TRUE.
      puncht = 1d10          ! this ensures that all the sheet
                             ! will be associated with a dbline
    END IF

    CALL listen('dbea3d')    !read keyword line
    IF (exists('DELETE')) THEN     !if there are pairs to be deleted
      WRITE(lures,"(' Deleted drawbead pairs ',/)",ERR=9999)
      DO
        CALL listen('dbea3d')      !read a line
        IF (exists('ENDDEL') ) EXIT  !end of list found
        name = get_name(posin=1,stype='DBEA')  !db name to be deleted
        WRITE(lures,"('Pair to delete: ',A)",ERR=9999) TRIM(name)
        CALL srch_db3pair (headp, anter, pair, name, found) !search in list
        IF (.NOT.found) THEN !send warning message
          WRITE(*,"('drawbeads: Specified pair does not exis: ',a)")TRIM(name)
          !CALL runend ('dbea3d:Specified pair does not exis')
        ELSE                 !delete db line
          CALL del_db3pair (headp, tailp, anter, pair)
          ndbea = ndbea - 1        !diminish number of DB lines
          CALL delete_name(name,6) !delete name from list
        END IF
      END DO
    ELSE
      backs = .TRUE.
    END IF
  ELSE
    backs = .TRUE.
  ENDIF

  ndbea = ndold ! no initial calculations necessary for nold db lines

  CALL listen('dbea3d')   !read a line
  ! read drawbead control information
  IF (  exists('INTERF')) THEN   !if data exist
    CALL db3con           !read DB control information
  ELSE
    backs = .TRUE.
  END IF
  ! now read drawbead lines connectivities and associated sheets
  CALL listen ('dbea3d') !read a line
  IF ( exists('SEGMEN')) THEN         !keyword to initiate data
    ALLOCATE (iwork(npoin))        !auxiliar array
    DO
      CALL listen ('dbea3d')          !read a lines
      IF (exists('ENDSEG'))   EXIT       !keyword to end data
      name = get_name(posin=1,stype='DBEA') !get drawbead line name
      WRITE(lures,"('Drawbead pair: ',A)",ERR=9999) TRIM(name)     !echo
      CALL srch_db3pair (headp, anter, pair, name, found)          !look in list
      IF (.NOT.found) CALL runend ('Dbea3d: Unknown drawbead pair !')  !data error
      ndbea = ndbea + 1                  !increase number of drawbead lines
      nnode = pair%iparm (5)             !recover number of nodes in connectivities
      CALL new_srf(surface)              !allocate surface for sheet
      ! read db line and surface connectivities
      CALL inidb0 (ndbpoin,nnode,nelem,label, surface,pair%iparm(8))
      pair%iparm(1) = ndbpoin            !number of points along the db line
      pair%sname    = surface%sname      !sheet name
      CALL alloc_db3 (pair, mst, nmn, nrtm, ndbpoin, ncrosi ,1) !lnd,xl,daz
      CALL inidb0a (npoin, nelem, coora, iwork, ndbpoin, nrtm, &
                   icode, pair%xl, surface%head, puncht,pair%lnd)
      CALL alloc_db3 (pair, mst, nmn, nrtm, ndbpoin, ncrosi ,2) !
      pair%iparm(2) = nrtm  !number of elements to be checked
      CALL inidb1 (iwork, nrtm, nmn, &   !generate array to be checked
                   pair%irect, surface%head, icode)
      pair%iparm (3) = nmn  !number of nodes included in checked elements
      CALL alloc_db3 (pair, mst, nmn, nrtm, ndbpoin, ncrosi ,3) !nseg & msr
    END DO
    DEALLOCATE (iwork)  !release auxiliar array
  ELSE
    backs = .TRUE.
  END IF
  CALL listen ('dbea3d')
  IF ( .NOT. exists('ENDDRA')) backs = .TRUE.

  nconp (4) = ndbea ! analogical to contact routines

  RETURN
 9999 CALL runen2('')
END SUBROUTINE db3inp

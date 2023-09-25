SUBROUTINE db3ini (ndold)

  !                3d drawbead input routine

  USE ctrl_db, ONLY: npoin
  USE npo_db, ONLY: coora
  IMPLICIT NONE
  !dummy arguments
  INTEGER (KIND=4) :: ndold  !previous number of DB lines
  !local variables
  INTEGER (KIND=4), ALLOCATABLE :: iwork(:)  !auxiliar
  INTEGER (KIND=4):: n,       &
                     nrtm,    &
                     nmn,     &
                     mst,     &
                     ndbpoin, &
                     nnode,   &
                     ncrosi
  TYPE (db3pair_db), POINTER :: pair

  ALLOCATE (iwork(npoin))     !get memory for auxiliar array

  pair => headp               !point to first pair
  DO n = 1, ndold             !skip old pairs
    pair => pair%next
  END DO

  DO  !loop over new pairs
    ndbpoin = pair%iparm (1)  !number of nodes in the line
    nrtm    = pair%iparm (2)  !number of elements to check
    nmn     = pair%iparm (3)  !number of nodes on elements
    nnode   = pair%iparm (5)  !number of nodes per element
    CALL inidb2 (pair%irect, pair%msr, iwork,       &
                 pair%nseg, nmn, nrtm, mst)
    pair%iparm (4) = mst      !total number of elements to check around 
    CALL alloc_db3 (pair, mst, nmn, nrtm, ndbpoin, ncrosi ,4) !lmsr(mst)
    CALL inidb3 (pair%irect, pair%msr, iwork,          &
                 pair%nseg, pair%lmsr, coora, pair%xl, &
                 pair%ldbs, pair%nstold,pair%fdat, ndbpoin, &
                 nrtm, nmn, nnode, ncrosi)
    pair%iparm (7) = ncrosi
    CALL alloc_db3 (pair, mst, nmn, nrtm, ndbpoin, ncrosi ,5)
    CALL inidb4 ( pair%nstold, pair%nstini, pair%fdat, &
                  pair%fdaini, nrtm)

    IF ( ASSOCIATED (pair%next) ) THEN
      pair => pair%next
    ELSE
      EXIT
    ENDIF
  END DO
  DEALLOCATE (iwork)

  RETURN

END SUBROUTINE db3ini

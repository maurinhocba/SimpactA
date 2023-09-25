SUBROUTINE db3for ( )
!-----------------------------------------------------------------------
!
!     Compute Forces for  3D drawbeads
!
!-----------------------------------------------------------------------

  USE ctrl_db, ONLY: ndofn, npoin
  USE npo_db, ONLY : fcont,coora
  IMPLICIT NONE


  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s

  INTEGER (kind=4) ::  nrt, nmn, ndbpoin, nnode, ihflag
  REAL (kind=8) :: fmax, elmodul
  TYPE (db3pair_db), POINTER :: pair
!======================================================================



  pair => headp

  DO

    IF ( .NOT.ASSOCIATED (pair) ) EXIT
    pair%tdbf  = 0d0             !initializes forces
    ndbpoin = pair%iparm (1)
    nrt     = pair%iparm (2)
    nmn     = pair%iparm (3)
    nnode   = pair%iparm (5)
    ihflag  = pair%iparm (6)     !compute total forces
    fmax    = pair%rparm (1)
    elmodul = pair%rparm (2)

    CALL dbsol3 (pair%xl, nrt,  pair%ldbs, pair%irect,     &
                 fmax, elmodul, pair%nstnew, pair%nstold,   &
                 ndbpoin, pair%fdat,pair%daz,pair%lmsr, pair%msr, pair%nseg, &
                 nmn, ihflag, pair%tdbf, nnode, fcont, coora)


    pair => pair%next
  END DO

  RETURN
 9999 CALL runen2('')
END SUBROUTINE db3for

SUBROUTINE loaini ( )

!*** generate resultant load

  USE ctrl_db, ONLY:  ndime, ndofn, nload, npoin, ntype,neq,mscal
  USE outp_db, ONLY: iwrit
  USE curv_db, ONLY: getcun
  USE kinc_db
  USE lispa0, ONLY :  lures
  USE loa_db
  USE npo_db
  IMPLICIT NONE

  INTEGER (kind=4) :: il, chnode, icont, iplod, &
                      lodpt, nedge,nsurf,ip
  REAL    (kind=8) :: point(ndofn),nmass
  TYPE (loa_set), POINTER :: loas
  TYPE (loa_nod), POINTER :: loa

  INTERFACE
    INCLUDE 'elemnt.h'
!    INCLUDE 'ensvec.h'
    INCLUDE 'dedge1.h'
    INCLUDE 'dsurf1.h'
    INCLUDE 'waterb_ini.h'
  END INTERFACE

  !IF (ASSOCIATED(force)) DEALLOCATE (force, loadv, loass)
  IF( nload == 0 )RETURN
  !ALLOCATE( force(neq+1,nload+1), loadv(ndofn,npoin,nload), loass(MAX(nload,1)) )

  IF (iwrit == 1) WRITE (lures, &
       "(/,'  L O A D   A P P L I E D   I N   T H E   ', &
       &   'P R E S E N T   S T R A T E G Y ',/)",ERR=9999)
  loadv = 0d0
  force = 0d0
  loass = 1
  IF(iwrit == 1) WRITE(lures,"(//,' REFERENCE Load SETS ',//)",ERR=9999)

  loas => headls
  DO il=1,nload               !loop over reference load sets
    force(neq+1,il) = loas%factor
    loass(il) = getcun (loas%lbl) ! get curve number
                                 ! get curve number for follower load
    IF (loas%fltype == 'INFLAT' .OR. loas%fltype == 'TUBHYD')   &
      loas%flpar%lc = getcun (loas%flpar%ltcur)
    IF (loas%numfl > 0) loas%flpar%lc = getcun (loas%flpar%ltcur)
    IF (loas%fltype == 'WATERB')  CALL waterb_ini(loas%numfl,ndime,loas%flpar,loas%headf)
    iplod = loas%iplod
    IF(iwrit == 1) WRITE(lures,"(//,' REFERENCE Load SET ',A,/)",ERR=9999) TRIM(loas%lbl)

    IF(iwrit == 1 .AND. iplod >0) WRITE(lures,"(//,' Nodal loads ',/)",ERR=9999)
    loa => loas%headn
    DO icont = 1,iplod
      lodpt = loa%node
      point(1:ndofn) = loa%forc(1:ndofn)
      IF(iwrit == 1) WRITE(lures,"(5x,i5,6g10.3)",ERR=9999) lodpt,point
      lodpt = chnode(lodpt)
      IF(lodpt >= 1 .AND. lodpt <= npoin) THEN
        loadv(1:ndofn,lodpt,il) = loadv(1:ndofn,lodpt,il) + point
      ELSE
        WRITE(*,"('     ABNORMAL END OF EXECUTION   ')")
        WRITE(lures,"(5x,'ERROR IN LOAD INPUT DATA, non-existent', &
          & ' node =',i5,/,5x,'EXECUTION STOPPED'///)",ERR=9999) label(lodpt)
        CALL runend('LOADPL:POINT LOAD NOT IN THE RANGE ')
      END IF
      loa => loa%next
    END DO

    nedge = loas%nedge

!#if UNIX
!        IF(nedge > 0) CALL dedge1(nedge,ntype,
!     &            ndime,ndofn,npoin,iwrit,coord,loadv(1:,1,il),
!     &    label(1:),loas%heade )
!#else
    IF(nedge > 0) CALL dedge1(nedge,ntype,ndime,ndofn,npoin,iwrit, &
              coord,loadv(:,:,il),loas%heade )
!#endif

    nsurf = loas%nsurf

!#if UNIX
!        IF(nsurf > 0) CALL dsurf1(nsurf, ndime,ndofn,npoin,iwrit,
!     &            coord,loadv(1:,1,il), label(1:),loas%heads)
!#else
    IF(nsurf > 0) CALL dsurf1(nsurf, ndime,ndofn,npoin,iwrit, &
              coord,loadv(:,:,il), loas%heads)
!#endif
!***    Compute element loads (gravitational and temperature effects)

    igrav = loas%igrav
    gravy = loas%gravy
    IF(iwrit == 1) WRITE(lures,"(//, &
          & ' Gravitational loading code   ',i8,/ &
          & ' Gravitational constant       ',e15.7,/)",ERR=9999)igrav, gravy
    IF (igrav /= 0 .AND. gravy /= 0d0) THEN
      gv(1:ndime) = loas%gvect(1:ndime)
      IF(iwrit == 1) WRITE(lures,"(//, &
                &   ' Gravitation direction',/ &
                &   8x,'X-comp.',8x,'Y-comp.',8x,'Z-comp.'/ &
                &   3e15.7)",ERR=9999) gv
      ! compute Gravitational load using Mass matrix
      point(1:ndime) = gv*gravy/mscal
      DO ip=1,npoin
        nmass = SUM(emass(1:ndime,ip))/ndime
        IF( nmass > 0d0)  &
        loadv(1:ndime,ip,il) = loadv(1:ndime,ip,il) + nmass*point(1:ndime)
      END DO
    END IF

!c#if UNIX
!c#if PVM
!      CALL PvmShrMas(ndofn,loadv(1,1,il))  !actualiza la masa
!c#endif
!c#else
!!DEC$ if (PVM==1)
!    CALL PvmShrMas(ndofn,loadv(1,1,il))  !actualiza la masa
!!DEC$ endif
!c#endif

    loas => loas%next
  END DO

  CALL flushf( lures )

  DO il=1,nload               !loop over reference load sets

    IF(iwrit == 1) THEN
      WRITE(lures,"(//,' REFERENCE Loads, for SET ',i5,/)",ERR=9999) il
      IF(ndime == 2) WRITE(lures,"('  Nodal Load Vector '/, &
                                 & ' Node     X',11x,'Y',11x,'Rz')",ERR=9999)
      IF(ndime == 3) WRITE(lures,"('  Nodal Load Vector '/, &
        & ' node     X',11x,'Y',11x,'Z',11x,'Rx',10x,'Ry',10x,'Rz')",ERR=9999)
      DO icont=1,npoin
        IF(ANY(loadv(1:ndofn,icont,il) /= 0d0))WRITE(lures, &
          & "(i5,8e12.4)",ERR=9999) label(icont),loadv(1:ndofn,icont,il)
      END DO
    END IF

!c#if UNIX
!        CALL ensvec (ndofn*npoin,ifpre(1:,1),loadv(1:,1,il),
!     &               force(:,il),npsdf,nesdf,ftsdf)
!c#else
    CALL ensve1 (ndofn*npoin,ifpre(1,1),loadv(1,1,il), &
                force(1,il),npsdf(1),nesdf(1),ftsdf(1))
!c#endif
  END DO

  RETURN
 9999 CALL runen2('')
END SUBROUTINE loaini

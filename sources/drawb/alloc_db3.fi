SUBROUTINE alloc_db3 (pair, mst, nmn, nrtm, ndbpoin, ncrosi,itask)
  ! allocates an array and initializes
  USE ctrl_db, ONLY : npoin
  IMPLICIT NONE
  ! dummy arguments
  INTEGER (kind=4) :: itask, & !task
                      mst, nmn, nrtm, ndbpoin, ncrosi
  TYPE (db3pair_db) :: pair    !pair of the array to be allocated

  INTERFACE
    INCLUDE 'inc_vint.h'
    INCLUDE 'inc_mint.h'
    INCLUDE 'inc_mrea.h'
  END INTERFACE

  SELECT CASE (itask)
  CASE ( 1)
    !x and y coordinates of the drawbead line points
    ALLOCATE (pair%lnd(2,ndbpoin-1),pair%xl(2,ndbpoin))
    pair%lnd(1:2,1:ndbpoin-1) = 0
    pair%xl(1:2,1:ndbpoin)  = 0d0

  CASE ( 2)
    !       nstnew  - codes of the state of the sheet elements with reference
    !                 to the drawbead line at the present step
    CALL inc_vint(pair%nstnew,ninc=nrtm)   !Allocate or resize 'nstnew'
    !       nstold  - codes of the state of the sheet elements with reference
    !                 to the drawbead line at the previous step
    CALL inc_vint(pair%nstold,ninc=nrtm)   !Allocate or resize 'nstold'
    !       ldbs    - drawbead line segments intersecting sheet element sides
    IF (ASSOCIATED(pair%ldbs)) THEN
      CALL inc_mint(pair%ldbs,nfil=0,ncol=nrtm)   !Resize 'ldbs'
    ELSE
      CALL inc_mint(pair%ldbs,nfil=2,ncol=nrtm)   !Allocate 'ldbs'
    END IF
    !       irect   - connectivities of sheet elements
    IF (ASSOCIATED(pair%irect)) THEN
      CALL inc_mint(pair%irect,nfil=0,ncol=nrtm)  !Resize 'irect'
    ELSE
      CALL inc_mint(pair%irect,nfil=4,ncol=nrtm)  !Allocate 'irect'
    END IF
    !       fdat    - real parameters defining the state of the sheet elements
    !                 at the previous step
    IF (ASSOCIATED(pair%fdat)) THEN
      CALL inc_mrea(pair%fdat,nfil=0,ncol=nrtm)   !Resize 'fdat'
    ELSE
      CALL inc_mrea(pair%fdat,nfil=3,ncol=nrtm)   !Allocate 'fdat'
    END IF
    ALLOCATE(pair%daz(nrtm))
    pair%daz(1:nrtm) = .FALSE.

  CASE ( 3)
    !       msr     - sheet node numbers
    !       nseg    - pointers to sheet nodes data in lmsr array
    ALLOCATE (pair%nseg(nmn+1), pair%msr(nmn))
    pair%nseg(1:nmn+1) = 0   ; pair%msr(1:nmn) = 0

  CASE ( 4)
    !       lmsr    - sheet segment (element) numbers surrounding each node
    ALLOCATE (pair%lmsr(mst))
    pair%lmsr(1:mst) = 0

  CASE ( 5)
    !       nstini  - codes of the state of the sheet elements with reference
    !                 to the drawbead line when they are introduced
    !       fdaini  - local natural coordinates (ss,tt) of the midpoint
    !                 of a drawbead segment intersecting a shell element
    !                 at the initial configuration
    ALLOCATE ( pair%nstini(MAX(1,ncrosi)), pair%fdaini(2,MAX(1,ncrosi)) )
    pair%nstini = 0 ; pair%fdaini = 0d0

  END SELECT

  RETURN
END SUBROUTINE

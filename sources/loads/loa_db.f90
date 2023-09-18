 MODULE loa_db
  USE param_db,ONLY: midn, mnam
  USE hydf_db
  IMPLICIT NONE
  SAVE

  INTEGER (kind=4) ::  igrav  !gravity flag
  REAL (kind=8) :: gv(3),gravy !gravity direction and value (auxiliars)

  PRIVATE:: init_flp

  INTEGER (kind=4) ::  ifloa=0      ! Number of sets with follower load

  ! Derived type for a list containing nodal forces
  TYPE loa_nod
    INTEGER (kind=4) :: node           !node label
    REAL (kind=8) :: forc(6)           !point force and moment components
    TYPE (loa_nod), POINTER :: next    !pointer to next node
  END TYPE loa_nod

  ! Derived type for a list containing edge load
  TYPE edg_nod
    CHARACTER (len=1) :: systm2        ! nodal system (local or global)
    INTEGER (kind=4) :: nn2            ! 2 or 3, No of nodes defining the edge
    INTEGER (kind=4) :: lnod2(3)       ! node labels
    REAL (kind=8) :: xa(3)             ! ???
    REAL (kind=8) :: press2(3,3)       ! nodal press components
    TYPE (edg_nod), POINTER :: next    ! pointer to next segment
  END TYPE edg_nod

  ! Derived type for a list containing surface load
  TYPE srf_nod
    CHARACTER (len=1) :: systm3        ! nodal system (local or global)
    INTEGER (kind=4) :: nn3            ! No of nodes defining the surface segm.
    INTEGER (kind=4) :: lnod3(6)       ! node labels
    REAL (kind=8) :: direc3(3)         ! load direction
    REAL (kind=8) :: press3(6)         ! nodal pressure
    TYPE (srf_nod), POINTER :: next    ! pointer to next segment
  END TYPE srf_nod

  ! Derived type for a list containing follower load
  TYPE foll_seg
    INTEGER (kind=4) :: nnode        !No of nodes defining the segment
    INTEGER (kind=4) :: lnofl(4)     !node labels
    REAL (kind=8) :: fload           !press
    TYPE (foll_seg), POINTER :: next !pointer to next segment
  END TYPE foll_seg

  ! Derived type for parameters of:
  ! inflated structures and
  ! tube hydroforming
  ! used to scale follower load: dp = K (Q *dt - dV) / V
  ! both require volume computation
  TYPE flpar
    ! used in both cases
    CHARACTER (len=mnam) :: ltcur    ! curve lbl for flux/temp control
    INTEGER (kind=4) :: lc        ! curve no. for flux/temp control=>facts
    INTEGER (kind=4) :: nodcen      ! node number of a center point for
                                    ! volume computation
    REAL (kind=8) :: timini         ! time to swap press/flux(temp)control
    REAL (kind=8) :: vol            ! volume (from previous step)
    REAL (kind=8) :: press          ! pressure (from previous step)
    ! for tube hydroforming
    REAL (kind=8) :: kvol           ! K - bulk modulus
    REAL (kind=8) :: qref           ! reference flux,  Q = qref * facts(t)
    ! for inflatable structures
    REAL (kind=8) :: pext           ! exterior pressure
    REAL (kind=8) :: volin          ! initial volume (at t=timini)
  END TYPE flpar

  ! FLPAR for water bags
  !  CHARACTER: ltcur     !
  !  INTEGER  : lc        !
  !  INTEGER  : nodcen    !01: symmetric 02: reference  free surface (no gas pressure)
  !                       !11: symmetric 12: reference  full (constant)
  !                       !21: symmetric 22: reference  filling (changes with curve)
  !  REAL : timini        ! x-coordinate to compute volume
  !  REAL : vol           ! present volume
  !  REAL : press         ! y coordinate
  !  REAL : kvol          ! K - bulk modulus
  !  REAL : qref          ! especific weigth
  !  REAL : pext          ! width at free surface
  !  REAL : volin         ! initial volume

  ! Derived type for storing load data
  TYPE loa_set
    CHARACTER (len=mnam) :: lbl ! set label
    !CHARACTER (len=mnam) :: lc  ! load curve label
    REAL (kind=8):: factor  ! scale factor
    ! point load
    INTEGER (kind=4) :: iplod
    TYPE (loa_nod), POINTER :: headn, tailn
    ! data for edge load
    INTEGER (kind=4) :: nedge
    TYPE (edg_nod), POINTER :: heade, taile
    ! data for surface load
    INTEGER (kind=4) :: nsurf
    TYPE (srf_nod), POINTER :: heads, tails
    ! gravitational loading
    INTEGER (kind=4) :: igrav
    REAL (kind=8) :: gravy
    REAL (kind=8) :: gvect(3)
    ! data for follower load
    INTEGER (kind=4)  :: numfl  ! number of segments loaded with fol. load
    CHARACTER (len=midn) :: fltype ! follower load type
                                   ! 'SIMPLE' - pressure controlled
                                   ! 'INFLAT' - inflatable structures,
                                   ! 'TUBHYD' - tube hydroforming
    TYPE (flpar) :: flpar          ! tube hydrofng and infl.structs. params
    LOGICAL :: fluid               !if pressure is provided by other program

    TYPE (foll_seg), POINTER :: headf, tailf
    TYPE (hydf_load), POINTER :: hydfl

    TYPE (loa_set), POINTER :: next
  END TYPE loa_set
  TYPE (loa_set), POINTER :: headls=>NULL(), taills=>NULL()
CONTAINS

  !-----------------------------------------
  !Routines for load sets
  !-----------------------------------------
  SUBROUTINE init_flp(flp)
  !This subroutine initialize parameters for inflated structures and tube hydroforming
  IMPLICIT NONE
  !Dummy variables
  TYPE(flpar):: flp

  !Initialize tube hydroforming and infl. structs. params
  flp%ltcur = ''     !Initialize curve lbl for flux/temp control
  flp%lc = 0         !Initialize curve no. for flux/temp control=>facts
  flp%nodcen = 0     !Initialize node number of a center point for volume computation
  flp%timini = 0d0   !Initialize time to swap press/flux(temp)control
  flp%vol = 0d0      !Initialize volume (from previous step)
  flp%press = 0d0    !Initialize pressure (from previous step)
  !Tube hydroforming variables
  flp%kvol = 0d0     !Initialize K - bulk modulus
  flp%qref = 0d0     !Initialize reference flux,  Q = qref * facts(t)
  !Inflatable structures variables
  flp%pext = 0d0     !Initialize exterior pressure
  flp%volin = 0d0    !Initialize initial volume (at t=timini)


  RETURN
  END SUBROUTINE init_flp

  !-----------------------------------------
  SUBROUTINE new_loas(new)
  !This subroutine allocate and initialize load set data
  IMPLICIT NONE
  !Dummy variables
  TYPE(loa_set),POINTER:: new

  ALLOCATE(new)  !get memory for new "loa_set"
  new%lbl = ''               !Initialize load set label
  new%factor = 0d0           !Initialize scale factor
  !Gravity load
  new%igrav = 0              !Initialize flag for gravity load
  new%gravy = 0d0            !Inititalize value of gravity load
  new%gvect(1:3) = 0d0       !Initialize unitary gravity vector
  !Point force
  new%iplod = 0              !Initialize load node counter
  NULLIFY(new%headn,new%tailn)   !Initialize load node list
  !Edge force
  new%nedge = 0              !Initialize edge load counter
  NULLIFY(new%heade,new%taile)   !Initialize edge load list
  !Surface forces
  new%nsurf = 0              !Intialize surface load counter
  NULLIFY(new%heads,new%tails)   !Initialize surface load list
  !Followed load
  new%numfl = 0              !Initialize number of segments loaded with fol. load
  new%fltype = ''            !Initialize folower load type
  new%fluid = .FALSE.  !Initialize
  NULLIFY(new%headf,new%tailf)   !Initialize follower load list
  NULLIFY(new%hydfl)             !Initialize definition of sheet hydroforming load
  CALL init_flp(new%flpar)   !Initialize tube hydroforming and infl. structs. params

  NULLIFY(new%next)   !Initialize next "loa_set" element in the list

  RETURN
  END SUBROUTINE new_loas

  SUBROUTINE add_loas(new, head, tail)
  !This subroutine adds data to the end of the list
  IMPLICIT NONE

    !--- Dummy arguments
    TYPE(loa_set),POINTER:: new, head, tail

    !Check if a list is empty
    IF (.NOT.ASSOCIATED(head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)
    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new
    ENDIF

  RETURN
  END SUBROUTINE add_loas

  !-----------------------------------------
  SUBROUTINE srch_loas(head, anter, posic, lbl, found)
  !This subroutine searches for a load set identified with a load curve
  IMPLICIT NONE

    !--- Dummy arguments
    LOGICAL:: found
    CHARACTER(len=*):: lbl ! Load set label
    TYPE(loa_set),POINTER:: head, anter, posic

    found = .FALSE.
    NULLIFY(posic,anter)
    !Check if a list is empty
    IF (ASSOCIATED(head)) THEN
      posic => head
      DO
        IF (TRIM(posic%lbl) == TRIM(lbl)) THEN
          found = .TRUE.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next)) THEN
          anter => posic
          posic => posic%next
        ELSE
          EXIT
        END IF
      END DO
    ENDIF
    IF (.NOT.found) NULLIFY(posic,anter)

  RETURN
  END SUBROUTINE srch_loas

  !-----------------------------------------
  SUBROUTINE delete_loas(loas)
  !Deletes a load set pointed with posic
  IMPLICIT NONE

    !--- Dummy variables
    TYPE(loa_set),POINTER:: loas

    CALL del_loa(loas%headn,loas%tailn)   !Deallocate load node list
    CALL del_edg(loas%heade,loas%taile)   !Deallocate edge load list
    CALL del_srf(loas%heads,loas%tails)   !Deallocate surface load list
    CALL del_foll(loas%headf,loas%tailf)  !Deallocate follower load list
    IF( ASSOCIATED (loas%hydfl) ) CALL del_hydfl(loas%hydfl)   !Deallocate variable of type 'hydf_load'
    DEALLOCATE(loas)  !Deallocate variable 'loa_set'

  RETURN
  END SUBROUTINE delete_loas

  !-----------------------------------------
  SUBROUTINE del_loas(head, tail, anter, posic)
  !Deletes a load set pointed with posic
  IMPLICIT NONE

    !--- Dummy variables
    TYPE(loa_set),POINTER:: head, tail, anter, posic

    IF (.NOT.ASSOCIATED(anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    ENDIF
    IF (.NOT.ASSOCIATED(posic%next) ) tail=>anter  ! if posic == tail
    CALL delete_loas(posic)
    !NULLIFY(anter)

  RETURN
  END SUBROUTINE del_loas

  !-----------------------------------------
  !Routines for point forces
  !-----------------------------------------
  SUBROUTINE new_loa(new)
  !This subroutine allocate and initialize variables of 'loa_nod' type
  IMPLICIT NONE
  !Dummy variables
  TYPE(loa_nod),POINTER:: new

  ALLOCATE(new)  !get memory for new "loa_nod"
  new%node = 0         !Initialize node label
  new%forc = 0d0  !Initialize point force and moment components
  NULLIFY(new%next)    !Initialize next "loa_nod" element in the list

  RETURN
  END SUBROUTINE new_loa

  !-----------------------------------------
  SUBROUTINE add_loa(new, head, tail)
  !This subroutine adds data to the end of the list
  IMPLICIT NONE

    !--- Dummy arguments
    TYPE(loa_nod),POINTER:: new, head, tail

    !Check if a list is empty
    IF (.NOT.ASSOCIATED(head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)
    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new
    ENDIF

  RETURN
  END SUBROUTINE add_loa

  !-----------------------------------------
  SUBROUTINE del_loa(head,tail)
  !Deallocates load nodes list
  IMPLICIT NONE

    !--- Dummy variables
    TYPE(loa_nod),POINTER:: head, tail
    !--- Local variables
    TYPE(loa_nod),POINTER:: loa

    NULLIFY(tail)
    DO
      IF (.NOT.ASSOCIATED(head)) EXIT
      loa => head%next
      DEALLOCATE(head)
      head => loa
    END DO

  RETURN
  END SUBROUTINE del_loa

  !-----------------------------------------
  !Routines for edge loads
  !-----------------------------------------
  SUBROUTINE new_edg(new)
  !This subroutine allocate and initialize variables of 'edg_nod' type
  IMPLICIT NONE
  !Dummy variables
  TYPE(edg_nod),POINTER:: new

  ALLOCATE(new)  !get memory for new "edg_nod"
  new%systm2 = ''             !Initialize nodal system (local or global)
  new%nn2 = 0                 !Initialize 2 or 3, No of nodes defining the edge
  new%lnod2(1:3) = 0          !Initialize node labels
  new%xa(1:3) = 0d0           !Initialize coordinates
  new%press2(1:3,1:3) = 0d0   !Intitalize nodal press components
  NULLIFY(new%next)  !Initialize next "edg_nod" element in the list

  RETURN
  END SUBROUTINE new_edg

  !-----------------------------------------
  SUBROUTINE add_edg(new, head, tail)
  !This subroutine adds data to the end of the list
  IMPLICIT NONE

    !--- Dummy arguments
    TYPE(edg_nod),POINTER:: new, head, tail

    !Check if a list is empty
    IF (.NOT.ASSOCIATED(head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)
    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new
    ENDIF

  RETURN
  END SUBROUTINE add_edg

  !-----------------------------------------
  SUBROUTINE del_edg(head,tail)
  !Deallocates edge load list
  IMPLICIT NONE

    !--- Dummy variables
    TYPE(edg_nod),POINTER:: head, tail
    !--- Local variables
    TYPE(edg_nod),POINTER:: edg

    NULLIFY(tail)
    DO
      IF (.NOT.ASSOCIATED(head)) EXIT
      edg => head%next
      DEALLOCATE(head)
      head => edg
    END DO

  RETURN
  END SUBROUTINE del_edg

  !-----------------------------------------
  !Routines for surface loads
  !-----------------------------------------
  SUBROUTINE new_srf(new)
  !This subroutine allocate and initialize variables of 'srf_nod' type
  IMPLICIT NONE
  !Dummy variables
  TYPE(srf_nod),POINTER:: new

  ALLOCATE(new)  !get memory for new "srf_nod"
  new%systm3 = ''         !Initialize nodal system (local or global)
  new%nn3 = 0             !Initialize No of nodes defining the surface segm.
  new%lnod3(1:6) = 0      !Initialize node labels
  new%direc3(1:3) = 0d0   !Initialize load direction
  new%press3(1:6) = 0d0   !Intitalize nodal pressure
  NULLIFY(new%next)  !Initialize next "srf_nod" element in the list

  RETURN
  END SUBROUTINE new_srf

  !-----------------------------------------
  SUBROUTINE add_srf(new, head, tail)
  !This subroutine adds data to the end of the list
  IMPLICIT NONE

    !--- Dummy arguments
    TYPE(srf_nod),POINTER:: new, head, tail

    !Check if a list is empty
    IF (.NOT.ASSOCIATED(head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)
    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new
    ENDIF

  RETURN
  END SUBROUTINE add_srf

  !-----------------------------------------
  SUBROUTINE del_srf(head,tail)
  !Deallocates surface load list
  IMPLICIT NONE

    !--- Dummy variables
    TYPE(srf_nod),POINTER::  head, tail
    !--- Local variables
    TYPE(srf_nod),POINTER:: srf

    NULLIFY(tail)
    DO
      IF (.NOT.ASSOCIATED(head)) EXIT
      srf => head%next
      DEALLOCATE(head)
      head => srf
    END DO

  RETURN
  END SUBROUTINE del_srf

  !-----------------------------------------
  !Rotines for follower loads
  !-----------------------------------------

  SUBROUTINE new_foll(new)
  !This subroutine allocate and initialize variables of 'foll_seg' type
  IMPLICIT NONE
  !Dummy variables
  TYPE(foll_seg),POINTER:: new

  ALLOCATE(new)  !get memory for new "foll_seg"
  new%nnode = 0        !Initialize No of nodes defining the segment
  new%lnofl(1:4) = 0   !Initialize node labels
  new%fload = 0d0      !Initialize press
  NULLIFY(new%next)  !Initialize next "foll_seg" element in the list

  RETURN
  END SUBROUTINE new_foll

  SUBROUTINE add_foll(new, head, tail)
  !This subroutine adds data to the end of the list
  IMPLICIT NONE

    !--- Dummy arguments
    TYPE(foll_seg),POINTER:: new, head, tail

    !Check if a list is empty
    IF (.NOT.ASSOCIATED(head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)
    ELSE
      !add a segment to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new
    ENDIF

  RETURN
  END SUBROUTINE add_foll

  !-----------------------------------------
  SUBROUTINE del_foll(head,tail)
  !Deallocates follower load list
  IMPLICIT NONE

    !--- Dummy variables
    TYPE(foll_seg),POINTER:: head, tail
    !--- Local variables
    TYPE(foll_seg),POINTER:: seg

    NULLIFY(tail)
    DO
      IF (.NOT.ASSOCIATED(head)) EXIT
      seg => head%next
      DEALLOCATE(head)
      head => seg
    END DO

  RETURN
  END SUBROUTINE del_foll

  SUBROUTINE dump_lo (ndime,ndofn,nload)

  !   dumps load database

  IMPLICIT NONE

  INTEGER (kind=4) :: ndime,ndofn,nload

  INTEGER (kind=4) :: i,j
  TYPE (loa_nod), POINTER :: loa
  TYPE (edg_nod), POINTER :: edg
  TYPE (srf_nod), POINTER :: srf
  TYPE (foll_seg),POINTER :: seg
  TYPE (loa_set), POINTER :: loas
  !------------------------------------------------------------

  WRITE (50,ERR=9999) ifloa
  IF (nload > 0) THEN
    loas => headls
    DO i=1,nload
      WRITE (50,ERR=9999) loas%factor, loas%lbl, loas%iplod
      loa => loas%headn
      DO j=1,loas%iplod
        WRITE (50,ERR=9999) loa%node, loa%forc(1:ndofn)
        loa => loa%next
      END DO
      WRITE (50,ERR=9999) loas%nedge
      edg => loas%heade
      DO j=1,loas%nedge
        WRITE (50,ERR=9999) edg%systm2, edg%nn2, edg%xa(1:ndime),     &
                   edg%lnod2(1:3), edg%press2(1:ndime,1:3)
        edg => edg%next
      END DO
      WRITE (50,ERR=9999) loas%nsurf
      srf => loas%heads
      DO j=1,loas%nsurf
        WRITE (50,ERR=9999) srf%systm3, srf%nn3, srf%direc3(1:ndime), &
                   srf%lnod3(1:6), srf%press3(1:6)
        srf => srf%next
      END DO

      WRITE (50,ERR=9999) loas%numfl, loas%fltype, loas%flpar, loas%fluid
      seg => loas%headf
      DO j=1,loas%numfl
        WRITE (50,ERR=9999) seg%nnode, seg%lnofl, seg%fload
        seg => seg%next
      END DO

      IF (loas%fltype == 'HYDROF') CALL dump_hydfl (loas%hydfl)

      WRITE (50,ERR=9999) loas%igrav,loas%gravy,loas%gvect

      loas => loas%next
    END DO
  END IF

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE dump_lo

  SUBROUTINE rest_lo (ndime,ndofn,nload)
  !   restors load database
  IMPLICIT NONE

    !--- Dummy variables
    INTEGER(kind=4):: ndime, ndofn, nload
    !--- Local variables
    INTEGER(kind=4):: i, j
    TYPE(loa_nod),POINTER:: loa
    TYPE(edg_nod),POINTER:: edg
    TYPE(srf_nod),POINTER:: srf
    TYPE(foll_seg),POINTER :: seg
    TYPE(loa_set),POINTER:: loas


    READ(51) ifloa
    NULLIFY(headls)
    DO i=1,nload
      CALL new_loas(loas)
      READ(51) loas%factor, loas%lbl, loas%iplod
      DO j=1,loas%iplod
        ALLOCATE (loa)
        READ(51) loa%node, loa%forc(1:ndofn)
        CALL add_loa(loa, loas%headn, loas%tailn)
      END DO

      READ(51) loas%nedge
      DO j=1,loas%nedge
        CALL new_edg(edg)
        READ(51) edg%systm2, edg%nn2, edg%xa(1:ndime), edg%lnod2(1:3), edg%press2(1:ndime,1:3)
        CALL add_edg(edg, loas%heade, loas%taile)
      END DO

      READ(51) loas%nsurf
      DO j=1,loas%nsurf
        CALL new_srf(srf)
        READ(51) srf%systm3, srf%nn3, srf%direc3(1:ndime), srf%lnod3(1:6), srf%press3(1:6)
        CALL add_srf(srf, loas%heads, loas%tails)
      END DO

      READ(51) loas%numfl, loas%fltype, loas%flpar, loas%fluid
      DO j=1,loas%numfl
        CALL new_foll(seg)
        READ(51) seg%nnode, seg%lnofl, seg%fload
        CALL add_foll(seg, loas%headf, loas%tailf)
      END DO

      IF (loas%fltype == 'HYDROF') THEN
        CALL ini_hydfl(loas%hydfl,ndime)
        CALL rest_hydfl(loas%hydfl)
      END IF

      READ(51) loas%igrav,loas%gravy,loas%gvect

      CALL add_loas(loas, headls, taills)
    END DO

    IF (ifloa > 0) CALL wrtfl0()

  RETURN
  END SUBROUTINE rest_lo

  SUBROUTINE updlon_lo (nload,oldlb)

  !   updates load database

  USE curv_db, ONLY : del_cur
  IMPLICIT NONE

  INTEGER (kind=4) :: nload
  INTEGER (kind=4), POINTER :: oldlb(:)

  INTEGER (kind=4) ::  j,k,lab,chnode,n
  TYPE (loa_nod), POINTER :: loa,aloa
  TYPE (edg_nod), POINTER :: edg,aedg
  TYPE (srf_nod), POINTER :: srf,asrf
  TYPE (foll_seg),POINTER :: seg,aseg
  TYPE (loa_set), POINTER :: loas,anter
  LOGICAL :: keep,keeps
  CHARACTER (len=mnam) :: lbl
  !------------------------------------------------------------
  ! modify internal number and check if load set remains
  IF (nload > 0) THEN  !if exist loads
    NULLIFY(anter)
    loas => headls      !point to firs set
    DO                    !modify internal numeration
      IF( .NOT.ASSOCIATED(loas) )EXIT    !exit if all sets processed
      keeps = loas%igrav /= 0             !initializes flag to keep this set
      !   point loads
      loa => loas%headn
      NULLIFY(aloa)
      n = 0
      DO j=1,loas%iplod
        lab = chnode( loa%node )
        IF( lab > 0 )THEN
          n = n + 1
          aloa => loa
          loa => loa%next
        ELSE
          IF( n > 0 )THEN
            aloa%next => loa%next
            DEALLOCATE(loa)
            loa => aloa%next
          ELSE
            loas%headn => loas%headn%next
            DEALLOCATE(loa)
            loa => loas%headn
          END IF
        END IF
      END DO
      loas%iplod = n
      keeps = n > 0 .OR. keeps
      !   edge loads
      edg => loas%heade
      NULLIFY(aedg)
      n = 0
      DO j=1,loas%nedge
        keep = .TRUE.
        DO k=1,edg%nn2
          lab = chnode( edg%lnod2(k) )
          IF( lab == 0 ) keep = .FALSE.
        END DO
        IF( keep )THEN
          n = n + 1
          aedg => edg
          edg => edg%next
        ELSE
          IF( n > 0 )THEN
            aedg%next => edg%next
            DEALLOCATE(edg)
            edg => aedg%next
          ELSE
            loas%heade => loas%heade%next
            DEALLOCATE(edg)
            edg => loas%heade
          END IF
        END IF
      END DO
      loas%nedge = n
      keeps = n > 0 .OR. keeps
      !   surface loads
      srf => loas%heads
      NULLIFY(asrf)
      n = 0
      DO j=1,loas%nsurf
        keep = .TRUE.
        DO k=1,srf%nn3
          lab = chnode( srf%lnod3(k) )
          IF( lab == 0 ) keep = .FALSE.
        END DO
        IF( keep )THEN
          n = n + 1
          asrf => srf
          srf => srf%next
        ELSE
          IF( n > 0 )THEN
            asrf%next => srf%next
            DEALLOCATE(srf)
            srf => asrf%next
          ELSE
            loas%heads => loas%heads%next
            DEALLOCATE(srf)
            srf => loas%heads
          END IF
        END IF
      END DO
      loas%nsurf = n
      keeps = n > 0 .OR. keeps
      ! follower loads
      seg => loas%headf
      NULLIFY(aseg)
      n = 0
      DO j=1,loas%numfl
        keep = .TRUE.
        DO k=1,seg%nnode
          lab = oldlb( seg%lnofl(k) )
          lab = chnode( lab )
          IF( lab == 0 ) keep = .FALSE.
        END DO
        IF( keep )THEN
          n = n + 1
          aseg => seg
          seg => seg%next
        ELSE
          IF( n > 0 )THEN
            aseg%next => seg%next
            DEALLOCATE(seg)
            seg => aseg%next
          ELSE
            loas%headf => loas%headf%next
            DEALLOCATE(seg)
            seg => loas%headf
          END IF
        END IF
      END DO
      loas%numfl = n
      keeps = n > 0 .OR. keeps
      IF (TRIM(loas%fltype) == 'INFLAT' .OR. TRIM(loas%fltype) == 'TUBHYD' ) THEN
        n = oldlb(loas%flpar%nodcen)     !original label
        loas%flpar%nodcen = chnode(n)         !present label
      END IF

      IF( keeps )THEN            !if load set remains
        anter => loas
        taills => loas
        loas => loas%next
      ELSE                       !delete load set
        lbl = loas%lbl              !set label
        CALL del_loas(headls, taills, anter, loas)
        CALL del_cur (lbl)          !delete curve data associated LC
        nload = nload - 1           !correct counter
        IF( ASSOCIATED(anter) )THEN   !point to next set
          loas => anter%next
        ELSE
          loas => headls
        END IF
        !delete set
      END IF
    END DO
  END IF

  RETURN
  END SUBROUTINE updlon_lo


 END MODULE loa_db

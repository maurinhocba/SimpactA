MODULE cont1_db
   USE param_db,ONLY: mnam,midn
   USE ctrl_db, ONLY : therm,ndoft,npoin,nconp,npoio,initial_displacements
   USE c_input, ONLY : openfi,exists,getrea,backs,get_name,getint,param,words, &
                       listen
   USE esets_db, ONLY : add_name,delete_name
   USE lispa0, ONLY : lures
   USE name_db,ONLY: rsfprj, rsfstg
   USE npo_db, ONLY : label,coora,iftmp,oldlb,coorc,euler,emass,fcont,coord
   USE surf_db, ONLY : surfa,new_srf,store_segs,dalloc_srf
  IMPLICIT NONE

  ! Global Variables for contact 1 Algorithm

  INTEGER (kind=4), PARAMETER :: &
     nisdb = 3, & ! number of integer values for each slave node
     nrsdb = 2    ! number of real values for each slave node
  INTEGER (kind=4) :: &
     nsurf = 0, & ! Number of contact surfaces
     npair = 0    ! Number of contact surfaces pairs
  REAL (kind=8) :: &
     ctime    ! time increment for force computation
  REAL (kind=8), ALLOCATABLE :: &
     surtf(:,:) !(ndime,npair) Total forces on each pair interface
  LOGICAL :: print_data = .FALSE.  !print initial penetrations on surfaces (GiD files)

!** Derived type for contact 1 database **!

  TYPE pair1_db
    ! For each contact pair the following variables are allocated

    CHARACTER (len=mnam) :: &
      pname,  & ! pair name
      master, & ! master surface name
      slave     ! slave surface name
    INTEGER (kind=4)  :: &
      imast,  & ! Master surface internal number
      islav,  & ! Slave surface internal number
      indcon, & ! contact type 0 1 2 forces on both, slave only, master only
      ncnod     ! number of nodes in slave surface
    REAL (kind=8) :: &
      npenal,  & ! Normal penalty coeff
      tpenal,  & ! Tangential penalty coeff
      static,  & ! Static Friction coefficient
      kinet ,  & ! Kinetic Friction coefficient
      start,   & ! Activation time
      end        ! Deactivation time

    LOGICAL :: prev,  & !.TRUE. if pair active in previous step
               press    !.TRUE. if nodal forces are stored

    ! Integer values for slave Data_base
    INTEGER (kind=4), POINTER ::   issdb(:,:)  ! ISSDB(nisdb,ncnod)
                     !(1) = nears    segment with projection
                     !(2) = nearo    previous segment with proy

    ! Real values for slave Data_base
    REAL (kind=8), POINTER ::    rssdb(:,:) !RSSDB(nrsdb,ncnod)
                            !(1) = rp       onset local coord
                            !(2) = gap      penetration in previous step

    REAL (kind=8), POINTER :: presn(:) !presn(ncnod) normal nodal force

    TYPE (pair1_db), POINTER :: next  ! Pointer to next pair

  END TYPE pair1_db

  TYPE (pair1_db), POINTER :: &
      headp,  & ! Pointer to first pair
      tailp     ! Pointer to last pair

  ! Derived type for the surface database
  TYPE surf1_db
    CHARACTER (len=mnam) :: sname    ! surface name
    LOGICAL :: iwrit,  & ! .T. : Save surface for Post Process
               imcod,  & ! Code for master surface
               iscod     ! Code for slave surface
               ! .F., Does not act as a master/slave surf. for any pair
               ! .T., Act as a master/slave surf. for some pair
!    LOGICAL :: curved    !treat surface as faceted or curved

    INTEGER (kind=4)  :: &
       ncnod,    & ! Number of nodes defining a surface
       nsegm       ! Number of segments in the surf

    REAL (kind=8) :: density, & !line density, to compute mass
                     diam       !diameter

    INTEGER (kind=4), POINTER :: &
       lcnod(:),   & !(ncnod) list of nodes in the surface
       lcseg(:,:)    !(2,nsegm) surface connectivities

!    REAL (kind=8), POINTER :: &
!       xc(:,:),  &   !(2,nsegm)  coordinates of the segment center
!       tc(:,:),  &   !(2,nsegm)  normal (outward) at segment center
!       cu(:),    &   !(nsegm)  surface curvatures
!       tn(:,:)       !(2,ncnod)  normal (outward) at nodes

    TYPE (surf1_db), POINTER :: next        ! pointer to next surface
  END TYPE surf1_db

  TYPE (surf1_db), POINTER :: &
      shead,  & ! pointer to first surface
      stail     ! pointer to last surface

CONTAINS

  !************    pair managment routines ************

  SUBROUTINE ini_cont1 (head, tail)
    !Initialize the contact 1 PAIRS database
    IMPLICIT NONE

       !Dummy arguments
    TYPE (pair1_db), POINTER :: head, tail

    NULLIFY (head, tail)
    RETURN
  END SUBROUTINE ini_cont1

  SUBROUTINE add_pair1 (new, head, tail)
    !This subroutine adds a pair dbase to the end of the list

    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair1_db), POINTER :: new, head, tail

       !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE   !add a pair to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    END IF

    RETURN
  END SUBROUTINE add_pair1

  SUBROUTINE srch_pair1 (head, anter, posic, name, found)
    !Searches for a pair named "name"

    IMPLICIT NONE
       !Dummy arguments
    LOGICAL :: found
    CHARACTER(len=*):: name ! set name
    TYPE (pair1_db), POINTER :: head, anter, posic

    found = .FALSE.
    NULLIFY (posic,anter)
       !Check if a list is empty
    IF (ASSOCIATED (head)) THEN
      posic => head
      DO
        IF(posic%pname == name) THEN
          found = .TRUE.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next) ) THEN
          anter => posic
          posic => posic%next
        ELSE
          EXIT
        END IF
      END DO
    END IF
    IF (.NOT.found) NULLIFY (posic,anter)

    RETURN
  END SUBROUTINE srch_pair1

  SUBROUTINE del_pair1 (head, tail, anter, posic)
    !Deletes a pair pointed with posic

    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair1_db), POINTER :: head, tail, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    END IF
       !if posic == tail
    IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
    CALL dalloc_pair1 (posic)
    NULLIFY (anter)

    RETURN
  END SUBROUTINE del_pair1

  SUBROUTINE dalloc_pair1 (pair)
    !Deallocates a pair
    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair1_db), POINTER :: pair

    DEALLOCATE ( pair%issdb, pair%rssdb )
    IF( pair%press ) DEALLOCATE( pair%presn )
    DEALLOCATE (pair)
    RETURN
  END SUBROUTINE dalloc_pair1

  SUBROUTINE new_pair1 (pair)
    !Allocates a pair
    IMPLICIT NONE
       !Dummy arguments
    TYPE (pair1_db), POINTER :: pair

    ALLOCATE (pair)
    NULLIFY ( pair%issdb, pair%rssdb, pair%presn )
    RETURN
  END SUBROUTINE new_pair1

  !************    surface managment routines ************

  SUBROUTINE ini_srf1 (head, tail)
    !Initialize the contact 1 SURFACES database
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf1_db), POINTER :: head, tail

    NULLIFY (head, tail)
    RETURN
  END SUBROUTINE ini_srf1

  SUBROUTINE add1_srf (new, head, tail)
    !Adds a surface to the end of the list

    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf1_db), POINTER :: new, head, tail

       !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN
       !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE   !add a surface to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    END IF
    RETURN
  END SUBROUTINE add1_srf

  SUBROUTINE srch1_srf (head, anter, posic, name, found)
    !This subroutine searches for a surface named "name"

    IMPLICIT NONE
       !Dummy arguments
    LOGICAL, INTENT(OUT) :: found
    CHARACTER(len=*),INTENT(IN):: name ! set name
    TYPE (surf1_db), POINTER :: head, anter, posic
    ! INTENT(IN) :: head  begining of the data base
    ! INTENT(OUT) :: posic  pointer to searched surface
    ! INTENT(OUT) :: anter  pointer to previous surface

    found = .FALSE.
    NULLIFY (posic,anter)
       !Check if a list is empty
    IF (ASSOCIATED (head)) THEN
      posic => head
      DO
        IF(TRIM(posic%sname) == TRIM(name)) THEN
          found = .TRUE.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next) ) THEN
          anter => posic
          posic => posic%next
        ELSE
          EXIT
        END IF
      END DO
    END IF
    IF (.NOT.found) NULLIFY (posic,anter)
    RETURN
  END SUBROUTINE srch1_srf

  SUBROUTINE del1_srf (head, tail, anter, posic)
    !Deletes a surface pointed with posic
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf1_db), POINTER :: head, tail, anter, posic
    ! INTENT(IN OUT) :: head  begining of the data base
    ! INTENT(IN OUT) :: tail  end of the data base
    ! INTENT(IN) :: posic  pointer to surface to delete
    ! INTENT(IN OUT) :: anter  pointer to previous surface

    IF (.NOT.ASSOCIATED (anter)) THEN   !IF deleled surface is the first
      head => posic%next                !head ==> 2nd surface
    ELSE
      anter%next => posic%next          !previous points to next
    END IF
       !if posic == tail                !If deleted surface is the last
    IF (.NOT.ASSOCIATED (posic%next) ) tail => anter  !
    CALL dallo1_srf (posic)   !release memory
    NULLIFY (anter)           !what for ?
    RETURN
  END SUBROUTINE del1_srf

  SUBROUTINE dallo1_srf (surface)
    !Deallocates a surface
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf1_db), POINTER :: surface

    IF( ASSOCIATED ( surface%lcnod )) DEALLOCATE ( surface%lcnod )
    IF( ASSOCIATED ( surface%lcseg )) DEALLOCATE ( surface%lcseg )
!    IF( ASSOCIATED ( surface%xc    )) DEALLOCATE ( surface%xc    )
!    IF( ASSOCIATED ( surface%tc    )) DEALLOCATE ( surface%tc    )
!    IF( ASSOCIATED ( surface%cu    )) DEALLOCATE ( surface%cu    )
!    IF( ASSOCIATED ( surface%tn    )) DEALLOCATE ( surface%tn    )

    DEALLOCATE (surface)
    RETURN
  END SUBROUTINE dallo1_srf

  SUBROUTINE new_surf1 (surf)
    !Allocates a surface
    IMPLICIT NONE
       !Dummy arguments
    TYPE (surf1_db), POINTER :: surf

    ALLOCATE (surf)
    NULLIFY ( surf%lcnod, surf%lcseg )  ! surf%nhseb, surf%xc, surf%tc, surf%cu, surf%tn )
    RETURN
  END SUBROUTINE new_surf1

 FUNCTION dist(x,y)
 REAL(kind=8) :: dist
 REAL(kind=8) :: x(3),y(3)
 REAL(kind=8) :: d(3)
 d = x - y
 dist = SQRT(DOT_PRODUCT(d,d))
 RETURN
 END FUNCTION dist

  INCLUDE 'beam_ini.fi'
  INCLUDE 'cdump1.fi'
  INCLUDE 'celmn1.fi'
  INCLUDE 'celmn1p.fi'
  INCLUDE 'chksu1.fi'
  INCLUDE 'chsur1.fi'
  INCLUDE 'cinpu1.fi'
  INCLUDE 'crest1.fi'
  INCLUDE 'cscf1a.fi'
  INCLUDE 'csdat1.fi'
  INCLUDE 'csinp1.fi'
  INCLUDE 'cspin1.fi'
  INCLUDE 'csrfc1.fi'
  INCLUDE 'cupdl1.fi'
  INCLUDE 'force1.fi'
  INCLUDE 'inpco1.fi'
  INCLUDE 'nearst.fi'
  INCLUDE 'projt1.fi'
  INCLUDE 'prsur1.fi'
  INCLUDE 'surms1.fi'
  INCLUDE 'upndc1.fi'

END MODULE cont1_db

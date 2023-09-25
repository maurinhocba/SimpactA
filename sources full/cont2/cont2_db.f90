MODULE cont2_db
   USE param_db,ONLY: mnam,midn
   USE ctrl_db, ONLY : therm,ndoft,npoin,nconp,npoio,ntype
   USE c_input, ONLY : openfi,exists,getrea,backs,get_name,getint,param,words, &
                       listen
   USE esets_db, ONLY : add_name,delete_name
   USE lispa0, ONLY : lures
   USE name_db,ONLY: rsfprj, rsfstg
   USE npo_db, ONLY : label,coora,oldlb,iftmp
   USE surf_db, ONLY : surfa,new_srf,store_segs,dalloc_srf
  IMPLICIT NONE

  ! Global Variables for contact 2D Algorithm

  INTEGER (kind=4), PARAMETER :: &
     nisdb = 3, & ! number of integer values for each slave node
     nrsdb = 2    ! number of real values for each slave node
  INTEGER (kind=4) :: &
     nsurf = 0, & ! Number of contact surfaces
     npair = 0    ! Number of contact surfaces pairs
  LOGICAL ::  wear     !.TRUE. if friction work is computed
  REAL (kind=8) :: &
     oldis, & ! maximum displacement expected
     ctime    ! time increment for force computation
  REAL (kind=8), ALLOCATABLE :: surtf(:,:) !(ndime,npair) Total forces on each pair interface
  REAL (kind=8), POINTER :: wwear(:)   !wwear(npoin) friction work
  LOGICAL :: print_data = .FALSE.      !print initial penetrations on surfaces (GiD files)

 !** Derived type for contact 2 database **!

   TYPE pair2_db
     ! For each contact pair the following variables are allocated

     CHARACTER (len=mnam) :: &
       pname,  & ! pair name
       master, & ! master surface name
       slave     ! slave surface name
     INTEGER (kind=4)  :: &
       imast,  & ! Master surface internal number
       islav,  & ! Slave surface internal number
       indcon, & ! contact type 0 1 2 forces on both, slave only, master only
       ncnod,  & ! number of nodes in slave surface
       mtsur,  & ! bottom, reversed, central or top surface for master (-2,-1,1,2)
       slsur,  & ! bottom, reversed, central or top surface for slave  (-2,-1,1,2)
       mtdof,  & ! temperature DOF associated to master surface
       sldof,  & ! temperature DOF associated to slave surface
       freq      ! frequency for global search
     REAL (kind=8) :: &
       npenal,  & ! Normal penalty coeff
       tpenal,  & ! Tangential penalty coeff
       static,  & ! Static Friction coefficient
       kinet ,  & ! Kinetic Friction coefficient
       cutoff,  & ! Maximum Addmisible gap
       gapinc,  & ! Maximum incremental gap
       hcont,   & ! h factor for heat conduction
       econt,   & ! e exponent for heat conduction
       hv   ,   & ! Vicker's hardness
       rec  ,   & ! Relative efusivity (slave)
       start,   & ! Activation time
       end        ! Deactivation time

     LOGICAL :: prev,  & !.TRUE. if pair active in previous step
                press, & !.TRUE. if nodal forces are stored
                cpress,& !.TRUE. if master surface is not a blank-holder
                auto,  & !.TRUE. if master and slave are the same surface
                wrink    !.TRUE. if wrinkles control is desired

     ! Integer values for slave Data_base
     INTEGER (kind=4), POINTER ::   issdb(:,:)  ! ISSDB(nisdb,ncnod)
                      !(1) = nears    segment with projection
                      !(2) = nearo    previous segment with proy
                      !(3) < 0  number of steps without contact
                      !    1 penetration and stuck   (Static friccion)
                      !    2 penetration and sliding (Kinetic friccion)

     ! Real values for slave Data_base
     REAL (kind=8), POINTER ::    rssdb(:,:) !RSSDB(nrsdb,ncnod)
                             !(1) = rp       onset local coord
                             !(2) = gap      penetration in previous step

     REAL (kind=8), POINTER :: presn(:) !presn(ncnod) normal nodal force
     REAL (kind=8), POINTER :: mingp(:) !mingp(ncnod) minimum gap

     TYPE (pair2_db), POINTER :: next  ! Pointer to next pair

   END TYPE pair2_db

   TYPE (pair2_db), POINTER :: &
       headp,  & ! Pointer to first pair
       tailp     ! Pointer to last pair

   ! Derived type for the surface database
   TYPE surf2_db
     CHARACTER (len=mnam) :: sname    ! surface name
     LOGICAL :: cxc,    & !.TRUE. compute triangle center coordinates
                bottom, & !.TRUE. bottom surface used by some pair
                confor, & !.TRUE. conforming surface
                iwrit,  & ! .T. : Save surface for Post Process
                press,  & ! .T. : binder pressure computed for some pair
                imcod,  & ! Code for master surface
                iscod     ! Code for slave surface
                ! .F., Does not act as a master/slave surf. for any pair
                ! .T., Act as a master/slave surf. for some pair
     LOGICAL :: curved, & !treat surface as faceted or curved
                auto      !surface will be tested for self contact

     INTEGER (kind=4)  :: &
        ncnod,    & ! Number of nodes defining a surface
        nsegm       ! Number of segments in the surf

     REAL (kind=8) :: density !surface density, to compute mass

     INTEGER (kind=4), POINTER :: &
        lcnod(:),   & !(ncnod) list of nodes in the surface
        lcseg(:,:), & !(2,nsegm) surface connectivities
        nhseg(:,:), & !(2,nsegm) connected segments to each segment
        lcseb(:,:), & !(2,nsegm) inverted surface connectivities (bottom)
 !       nr(:),      & !(npoin) inverted relation between nodes and local nodes
        nhseb(:,:)    !(2,nsegm) inverted connection (bottom)

     REAL (kind=8), POINTER :: &
        xc(:,:),  &   !(2,nsegm)  coordinates of the segment center
        tc(:,:),  &   !(2,nsegm)  normal (outward) at segment center
        cu(:),    &   !(nsegm)  surface curvatures
        tn(:,:),  &   !(2,ncnod)  normal (outward) at nodes
        area(:)       !(ncnod)    area associated to each node for pressure computation

     TYPE (surf2_db), POINTER :: next        ! pointer to next surface
   END TYPE surf2_db

   TYPE (surf2_db), POINTER :: &
       shead,  & ! pointer to first surface
       stail     ! pointer to last surface

 CONTAINS

   !************    pair managment routines ************

   SUBROUTINE ini_cont2 (head, tail)
     !Initialize the contact 2 PAIRS database
     IMPLICIT NONE

        !Dummy arguments
     TYPE (pair2_db), POINTER :: head, tail

     NULLIFY (head, tail)
     RETURN
   END SUBROUTINE ini_cont2

   SUBROUTINE add_pair2 (new, head, tail)
     !This subroutine adds a pair dbase to the end of the list

     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair2_db), POINTER :: new, head, tail

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
   END SUBROUTINE add_pair2

   SUBROUTINE srch_pair2 (head, anter, posic, name, found)
     !Searches for a pair named "name"

     IMPLICIT NONE
        !Dummy arguments
     LOGICAL :: found
     CHARACTER(len=*):: name ! set name
     TYPE (pair2_db), POINTER :: head, anter, posic

     found = .FALSE.
     NULLIFY (posic,anter)
        !Check if a list is empty
     IF (ASSOCIATED (head)) THEN
       posic => head
       DO
         IF(TRIM(posic%pname) == TRIM(name)) THEN
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
   END SUBROUTINE srch_pair2

   SUBROUTINE del_pair2 (head, tail, anter, posic)
     !Deletes a pair pointed with posic

     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair2_db), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
        !if posic == tail
     IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
     CALL dalloc_pair2 (posic)
     !NULLIFY (anter)         !what for

     RETURN
   END SUBROUTINE del_pair2

   SUBROUTINE dalloc_pair2 (pair)
     !Deallocates a pair
     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair2_db), POINTER :: pair

     IF( ASSOCIATED ( pair%issdb) ) DEALLOCATE ( pair%issdb, pair%rssdb )
     IF( pair%press )THEN
       IF( ASSOCIATED ( pair%presn) ) DEALLOCATE( pair%presn )
     END IF
     IF( pair%wrink ) THEN
       IF( ASSOCIATED ( pair%mingp) ) DEALLOCATE( pair%mingp )
     END IF
     DEALLOCATE (pair)
     RETURN
   END SUBROUTINE dalloc_pair2

   SUBROUTINE new_pair2 (pair)
     !Allocates a pair
     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair2_db), POINTER :: pair

     ALLOCATE (pair)
     NULLIFY ( pair%issdb, pair%rssdb, pair%presn )
     RETURN
   END SUBROUTINE new_pair2

   !************    surface managment routines ************

   SUBROUTINE ini_srf2 (head, tail)
     !Initialize the contact 2 SURFACES database
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf2_db), POINTER :: head, tail

     NULLIFY (head, tail)
     RETURN
   END SUBROUTINE ini_srf2

   SUBROUTINE ini2_srf (head, tail)
     !Initialize a list of surfaces
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf2_db), POINTER :: head, tail

     NULLIFY (head, tail)
     RETURN
   END SUBROUTINE ini2_srf

   SUBROUTINE add2_srf (new, head, tail)
     !Adds a surface to the end of the list

     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf2_db), POINTER :: new, head, tail

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
   END SUBROUTINE add2_srf

   SUBROUTINE srch2_srf (head, anter, posic, name, found)
     !This subroutine searches for a surface named "name"

     IMPLICIT NONE
        !Dummy arguments
     LOGICAL, INTENT(OUT) :: found
     CHARACTER(len=*),INTENT(IN):: name ! set name
     TYPE (surf2_db), POINTER :: head, anter, posic
     ! INTENT(IN) :: head  begining of the data base
     ! INTENT(OUT) :: posic  pointer to searched surface
     ! INTENT(OUT) :: anter  pointer to previous surface

     found = .FALSE.     !initializes
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
   END SUBROUTINE srch2_srf

   SUBROUTINE del2_srf (head, tail, anter, posic)
     !Deletes a surface pointed with posic
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf2_db), POINTER :: head, tail, anter, posic
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
     CALL dallo2_srf (posic)   !release memory
     !NULLIFY (anter)           !what for ?
     RETURN
   END SUBROUTINE del2_srf

   SUBROUTINE dallo2_srf (surface)
     !Deallocates a surface
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf2_db), POINTER :: surface

     IF( ASSOCIATED ( surface%lcnod )) DEALLOCATE ( surface%lcnod )
     IF( ASSOCIATED ( surface%lcseg )) DEALLOCATE ( surface%lcseg )
     IF( ASSOCIATED ( surface%nhseg )) DEALLOCATE ( surface%nhseg )
     IF( ASSOCIATED ( surface%lcseb )) DEALLOCATE ( surface%lcseb )
 !    IF( ASSOCIATED ( surface%nr ))    DEALLOCATE ( surface%nr )
     IF( ASSOCIATED ( surface%nhseb )) DEALLOCATE ( surface%nhseb )
     IF( ASSOCIATED ( surface%xc    )) DEALLOCATE ( surface%xc    )
     IF( ASSOCIATED ( surface%tc    )) DEALLOCATE ( surface%tc    )
     IF( ASSOCIATED ( surface%cu    )) DEALLOCATE ( surface%cu    )
     IF( ASSOCIATED ( surface%tn    )) DEALLOCATE ( surface%tn    )

     DEALLOCATE (surface)
     RETURN
   END SUBROUTINE dallo2_srf

   SUBROUTINE new_surf2 (surf)
     !Allocates a surface
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf2_db), POINTER :: surf

     ALLOCATE (surf)
     surf%bottom = .FALSE.
     surf%press  = .FALSE.
     surf%auto   = .FALSE.
     NULLIFY ( surf%lcnod, surf%lcseg, surf%nhseg, surf%lcseb, &
               surf%nhseb, surf%xc, surf%tc, surf%cu, surf%tn, &
               surf%area  )
     !NULLIFY ( surf%nr )
     RETURN
   END SUBROUTINE new_surf2

  INCLUDE 'ccheck2.fi'
  INCLUDE 'cdump2.fi'
  INCLUDE 'celmn2.fi'
  INCLUDE 'celmn2p.fi'
  INCLUDE 'chksu2.fi'
  INCLUDE 'chsur2.fi'
  INCLUDE 'cinpu2.fi'
  INCLUDE 'coffs2.fi'
  INCLUDE 'crest2.fi'
  INCLUDE 'cscf2a.fi'
  INCLUDE 'cscf2aa.fi'
  INCLUDE 'csdat2.fi'
  INCLUDE 'csinp2.fi'
  INCLUDE 'cspin2.fi'
  INCLUDE 'csrfc2.fi'
  INCLUDE 'csrfc2a.fi'
  INCLUDE 'cupdl2.fi'
  INCLUDE 'curve2.fi'
  INCLUDE 'cxctc2.fi'
  INCLUDE 'force2.fi'
  INCLUDE 'inpco2.fi'
  INCLUDE 'mastd2.fi'
  INCLUDE 'mastdb.fi'
  INCLUDE 'nearst.fi'
  INCLUDE 'projt2.fi'
  INCLUDE 'prsur2.fi'
  INCLUDE 'surms2.fi'
  INCLUDE 'upndc2.fi'

END MODULE cont2_db

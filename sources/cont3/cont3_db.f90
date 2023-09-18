 MODULE cont3_db
   USE param_db,ONLY: mnam,midn
   USE ctrl_db, ONLY : npoin,nconp,npoio ,therm,ndoft
   USE c_input, ONLY : openfi,exists,getrea,backs,get_name,getint,param,words, &
                       listen
   USE esets_db, ONLY : add_name,delete_name
   USE lispa0, ONLY : lures
   USE name_db,ONLY: rsfprj, rsfstg
   USE npo_db, ONLY : label,coora,oldlb,iftmp
   USE surf_db, ONLY : surfa,new_srf,store_segs,dalloc_srf
   IMPLICIT NONE

   ! Global Variables for contact 3D Algorithm

   INTEGER (kind=4), PARAMETER :: &
      nisdb = 3, & ! number of integer values for each slave node
      nrsdb = 3    ! number of real values for each slave node
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

 !** Derived type for contact 3 database **!

   TYPE pair3_db
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
       freq      ! Frequency for global search
     REAL (kind=8) :: &
       npenal, & ! Normal penalty coeff
       tpenal, & ! Tangential penalty coeff
       static, & ! Static Friction coefficient
       kinet , & ! Kinetic Friction coefficient
       cutoff, & ! Maximum Addmisible gap
       gapinc, & ! Maximum incremental gap
       hcont,  & ! h factor for heat conduction
       econt,  & ! e exponent for heat conduction
       hv   ,  & ! Vicker's hardness
       rec  ,  & ! Relative efusivity (slave)
       start,  & ! Activation time
       end       ! Deactivation time
     LOGICAL :: prev,  & !.TRUE. if pair active in previous step
                press, & !.TRUE. if nodal forces are stored
                cpress,& !.TRUE. if master surface is not a blank-holder
                auto,  & !.TRUE. if master and slave are the same surface
                wrink    !.TRUE. if wrinkles control is desired

     ! Integer values for slave Data_base
     INTEGER (kind=4), POINTER :: issdb(:,:)  ! ISSDB(nisdb,ncnod)
                !(1) = nears    segment with projection
                !(2) = nearo    previous segment with proy
                !(3) < 0  number of steps without contact
                !    1 penetration and stuck   (Static friccion)
                !    2 penetration and sliding (Kinetic friccion)

     ! Real values for slave Data_base
     REAL (kind=8), POINTER :: rssdb(:,:) !RSSDB(nrsdb,ncnod)
                       !(1) = rp       onset local coord
                       !(2) = sp
                       !(3) = gap      penetration in previous step
     REAL (kind=8), POINTER :: presn(:) !presn(ncnod) normal nodal force
     REAL (kind=8), POINTER :: mingp(:) !mingp(ncnod) minimum gap

     TYPE (pair3_db), POINTER :: next  ! Pointer to next pair

   END TYPE pair3_db

   TYPE (pair3_db), POINTER :: &
       headp,  & ! Pointer to first pair
       tailp     ! Pointer to last pair

   ! Derived type for the surface database
   TYPE surf3_db
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
        ncnod,   & ! Number of nodes defining a surface
        nsegm      ! Number of segments in the surf
     REAL (kind=8) :: density !surface density, to compute mass

     INTEGER (kind=4), POINTER :: &
        lcnod(:),   & !(ncnod) list of nodes in the surface
        lcseg(:,:), & !(3,nsegm) surface connectivities
        nhseg(:,:), & !(3,nsegm) connected segments to each segment
        lcseb(:,:), & !(3,nsegm) inverted surface connectivities (bottom)
        nr(:),      & !(npoin) inverted relation between nodes and local nodes
        nhseb(:,:)    !(3,nsegm) inverted connection (bottom)

     REAL (kind=8), POINTER :: &
        xc(:,:),  & !(3,nsegm)  coordinates of the segment center
        tc(:,:),  & !(3,nsegm)  normal (outward) at segment center
        cu(:,:),  & !(3,nsegm)  surface curvatures
        tn(:,:),  & !(3,ncnod)  normal (outward) at nodes
        area(:)     !(ncnod)    area associated to each node for pressure computation

     TYPE (surf3_db), POINTER :: next        ! pointer to next surface
   END TYPE surf3_db

   TYPE (surf3_db), POINTER :: &
       shead,  & ! pointer to first surface
       stail     ! pointer to last surface

 CONTAINS

   !************    pair managment routines ************

   SUBROUTINE ini_cont3 (head, tail)
     !Initialize the contact 3 PAIRS database
     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair3_db), POINTER :: head, tail

     NULLIFY (head, tail)
     RETURN
   END SUBROUTINE ini_cont3

   SUBROUTINE add_pair3 (new, head, tail)
     !This subroutine adds a pair dbase to the end of the list
     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair3_db), POINTER :: new, head, tail

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
   END SUBROUTINE add_pair3

   SUBROUTINE srch_pair3 (head, anter, posic, name, found)
     !Searches for a pair named "name"
     IMPLICIT NONE
        !Dummy arguments
     LOGICAL :: found
     CHARACTER(len=*):: name ! set name
     TYPE (pair3_db), POINTER :: head, anter, posic

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
   END SUBROUTINE srch_pair3

   SUBROUTINE del_pair3 (head, tail, anter, posic)
     !Deletes a pair pointed with posic
     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair3_db), POINTER :: head, tail, anter, posic

     IF (.NOT.ASSOCIATED (anter)) THEN
       head => posic%next
     ELSE
       anter%next => posic%next
     END IF
        !if posic == tail
     IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
     CALL dalloc_pair3 (posic)
     !NULLIFY (anter)  !what for
     RETURN
   END SUBROUTINE del_pair3

   SUBROUTINE dalloc_pair3 (pair)
     !Deallocates a pair
     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair3_db), POINTER :: pair

     IF( ASSOCIATED ( pair%issdb) ) DEALLOCATE ( pair%issdb, pair%rssdb )
     IF( pair%press )THEN
       IF( ASSOCIATED ( pair%presn) ) DEALLOCATE( pair%presn )
     END IF
     IF( pair%wrink ) THEN
       IF( ASSOCIATED ( pair%mingp) ) DEALLOCATE( pair%mingp )
     END IF
     DEALLOCATE (pair)
     RETURN
   END SUBROUTINE dalloc_pair3

   SUBROUTINE new_pair3 (pair)
     !Allocates a pair
     IMPLICIT NONE
        !Dummy arguments
     TYPE (pair3_db), POINTER :: pair

     ALLOCATE (pair)
     NULLIFY ( pair%issdb, pair%rssdb, pair%presn )
     RETURN
   END SUBROUTINE new_pair3

   !************    surface managment routines ************

   SUBROUTINE ini_srf3 (head, tail)
     !Initialize the contact 3 SURFACES database
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf3_db), POINTER :: head, tail

     NULLIFY (head, tail)
     RETURN
   END SUBROUTINE ini_srf3

   SUBROUTINE ini3_srf (head, tail)
     !Initialize a list of surfaces
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf3_db), POINTER :: head, tail

     NULLIFY (head, tail)
     RETURN
   END SUBROUTINE ini3_srf

   SUBROUTINE add3_srf (new, head, tail)
     !Adds a surface to the end of the list
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf3_db), POINTER :: new, head, tail

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
   END SUBROUTINE add3_srf

   SUBROUTINE srch3_srf (head, anter, posic, name, found)
     !This subroutine searches for a surface named "name"
     IMPLICIT NONE
        !Dummy arguments
     LOGICAL, INTENT(OUT) :: found
     CHARACTER(len=*),INTENT(IN):: name ! set name
     TYPE (surf3_db), POINTER :: head, anter, posic
     ! INTENT(IN) :: head  begining of the data base
     ! INTENT(OUT) :: posic  pointer to searched surface
     ! INTENT(OUT) :: anter  pointer to previous surface

     found = .FALSE.    !initializes
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
   END SUBROUTINE srch3_srf

   SUBROUTINE del3_srf (head, tail, anter, posic)
     !Deletes a surface pointed with posic
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf3_db), POINTER :: head, tail, anter, posic
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
     CALL dallo3_srf (posic)   !release memory
     !NULLIFY (anter)           !what for ?
     RETURN
   END SUBROUTINE del3_srf

   SUBROUTINE dallo3_srf (surface)
     !Deallocates a surface
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf3_db), POINTER :: surface

     IF( ASSOCIATED ( surface%lcnod )) DEALLOCATE ( surface%lcnod )
     IF( ASSOCIATED ( surface%lcseg )) DEALLOCATE ( surface%lcseg )
     IF( ASSOCIATED ( surface%nhseg )) DEALLOCATE ( surface%nhseg )
     IF( ASSOCIATED ( surface%lcseb )) DEALLOCATE ( surface%lcseb )
     IF( ASSOCIATED ( surface%nr ))    DEALLOCATE ( surface%nr )
     IF( ASSOCIATED ( surface%nhseb )) DEALLOCATE ( surface%nhseb )
     IF( ASSOCIATED ( surface%xc    )) DEALLOCATE ( surface%xc    )
     IF( ASSOCIATED ( surface%tc    )) DEALLOCATE ( surface%tc    )
     IF( ASSOCIATED ( surface%cu    )) DEALLOCATE ( surface%cu    )
     IF( ASSOCIATED ( surface%tn    )) DEALLOCATE ( surface%tn    )

     DEALLOCATE (surface)
     RETURN
   END SUBROUTINE dallo3_srf

   SUBROUTINE new_surf3 (surf)
     !Allocates a surface
     IMPLICIT NONE
        !Dummy arguments
     TYPE (surf3_db), POINTER :: surf

     ALLOCATE (surf)
     surf%bottom = .FALSE.
     surf%press  = .FALSE.
     surf%auto   = .FALSE.
     NULLIFY ( surf%lcnod, surf%lcseg, surf%nhseg, surf%lcseb, surf%nr , &
               surf%nhseb, surf%xc, surf%tc, surf%cu, surf%tn, surf%area )
     RETURN
   END SUBROUTINE new_surf3

  INCLUDE 'cdump3.fi'
  INCLUDE 'crest3.fi'
  INCLUDE 'celmn3.fi'
  INCLUDE 'celmn3p.fi'
  INCLUDE 'chckar.fi'
  INCLUDE 'chksu3.fi'
  INCLUDE 'chsur3.fi'
  INCLUDE 'cinpu3.fi'
  INCLUDE 'coffs3.fi'
  INCLUDE 'cupdl3.fi'
  INCLUDE 'cscf3a.fi'
  INCLUDE 'cscf3aa.fi'
  INCLUDE 'csdat3.fi'
  INCLUDE 'csinp3.fi'
  INCLUDE 'csmas3.fi'
  INCLUDE 'cspin3.fi'
  INCLUDE 'csrfc3.fi'
  INCLUDE 'csrfc3a.fi'
  INCLUDE 'curve3.fi'
  INCLUDE 'cxctc3.fi'
  INCLUDE 'force3.fi'
  !INCLUDE 'getno3.fi'
  INCLUDE 'inpco3.fi'
  INCLUDE 'mastd3.fi'
  INCLUDE 'mastda.fi'
  !INCLUDE 'mdsur3.fi'
  INCLUDE 'nears3.fi'
  INCLUDE 'projt3.fi'
  INCLUDE 'projtr.fi'
  INCLUDE 'prsur3.fi'
  INCLUDE 'surms3.fi'
  INCLUDE 'upndc3.fi'

 END MODULE cont3_db

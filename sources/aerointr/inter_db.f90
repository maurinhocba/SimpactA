
MODULE inter_db
  
  USE param_db,     ONLY:   mnam
  USE StructuresA,  ONLY:   tnode, tpanel, tsection,                        & ! types
                            blade1nodes, blade2nodes,  blade3nodes,         & ! variables
                            hubnodes,    nacellenodes, towernodes,          &
                            blade1panels, blade2panels,  blade3panels,      &
                            hubpanels,    nacellepanels, towerpanels,       &
                            blade1sections, blade2sections, blade3sections
  
  IMPLICIT NONE
  
  TYPE pair
    CHARACTER(len=mnam)             :: m_est        ! Setname in the pair
    CHARACTER(len=mnam)             :: m_aer        ! Millpart codename in the pair; may be one of:
                                                    !        blade1, blade2, blade3
                                                    !        tower, nacelle, hub
    ! Structural data -----------------------------------------------------------------------------
    CHARACTER(len=10)               :: eType        ! element type for structural model
    INTEGER (kind=4)                :: numbnodes    ! Number of different nodes in the set
    INTEGER (kind=4),   ALLOCATABLE :: nodes(:)     ! Internal labels of nodes in the set
    INTEGER (kind=4)                :: nelem        ! Number of structural elements
    INTEGER (kind=4),   ALLOCATABLE :: lnods(:,:)   ! Connectivities for structural elements
                                                    !       not necessary for 'RIGID' elements
                                                    !       (2,nelem) for 'BEAM ' elements
                                                    !       (3,nelem) for 'NBST ' elements
    INTEGER (kind=4),   ALLOCATABLE :: elmIndxs(:)  ! Index corresponding to element label in LNODS
                                                    !       for shell model
                                                    !       elmIndxs(lBound:uBound)
                                                    ! allocated in nbstConnects/nbstMstrs/mastrnods_sr
                                                    ! > bed%buse (or bed%blse) is a vector  with elements' external labels
                                                    ! > if connectivities are needed, they are stored in lnods
                                                    ! > but lnods column index do not agree with external labels stored in buse
                                                    ! > the connection between them is elmIndxs, used as follows
                                                    ! > lnodsColumnIndex = elmIndxs(elementExternalLabel)
    TYPE(bed),          POINTER     :: bed          ! Blade Elements Distribution - for shell model
    ! Aerodynamic nodes data ----------------------------------------------------------------------
    TYPE(tnode),        POINTER     :: nodsaer(:)   ! (nan) Aero. nodes information array
    INTEGER                         :: nan          ! Num. of aero. nodes.
    REAL,               ALLOCATABLE :: xita(:,:)    ! Projection coordinate of each aero. node
                                                    !        not necessary for 'RIGID' elements
                                                    !        (nan)    for 'BEAM ' elements
                                                    !        (6,nan)  for 'NBST ' elements
                                                    !                   3 local coords for 1 elem. in each surface - upper first
    INTEGER,            ALLOCATABLE :: mnods(:,:)   ! Master strut. nodes of each aero. node
                                                    !        (1,1)    for 'RIGID' elements
                                                    !        (2,nan)  for 'BEAM ' elements
                                                    !        (6,nan)  for 'NBST ' elements
                                                    !                   3 nodes for 1 elem. in each surface - upper first
    REAL,               ALLOCATABLE :: dists(:,:)   ! Initial distance vectors from master nodes to aero. node
                                                    !        (3,nan)  for 'RIGID' elements
                                                    !        (6,nan)  for 'BEAM ' elements
                                                    !        (2,nan)  for 'NBST ' elements
                                                    !                   normalized distance to each surface - upper first
    LOGICAL,            ALLOCATABLE :: reguDep(:)   ! (nan) Regular dependency .TRUE.
                                                    !       for shell model
    ! Aerodynamic control points data -------------------------------------------------------------
    TYPE(tpanel),       POINTER     :: panel(:)     ! (np) Panels information array
    TYPE(tsection),     POINTER     :: section(:)   ! (nsec) Used only for blade pairs
    INTEGER                         :: np           ! Num. of panels/Control Points
    REAL,               ALLOCATABLE :: cpxita(:,:)  ! Projection coordinate of each aero. node
                                                    !        not necessary for 'RIGID' elements
                                                    !        (np)    for 'BEAM ' elements
                                                    !        (6,np)  for 'NBST ' elements
                                                    !                   3 local coords for 1 elem. in each surface
    INTEGER,            ALLOCATABLE :: cpmnods(:,:) ! Master strut. nodes of each Control Point
                                                    !        (1,1)   for 'RIGID' elements
                                                    !        (2,np)  for 'BEAM ' elements
                                                    !        (6,np)  for 'NBST ' elements
    REAL,               ALLOCATABLE :: cpdists(:,:) ! Initial distance vectors from master nodes to Control Point
                                                    !        (3,np)  for 'RIGID' elements
                                                    !        (6,np)  for 'BEAM ' elements
                                                    !        (2,np)  for 'NBST ' elements
                                                    !                   normalized distance to each surface
    LOGICAL,            ALLOCATABLE :: cpReguDep(:) ! (np) Regular dependency .TRUE.
                                                    !       for shell model
    ! ---------------------------------------------------------------------------------------------
    TYPE(pair),         POINTER     :: next         ! Next pair in the list
  ENDTYPE pair
  
  TYPE bed ! Blade Elements Distribution - for shell model
    CHARACTER(len=mnam)             :: m_est        ! Setname in the pair - structural
    ! ---------------------------------------------------------------------------------------------
    INTEGER (kind=4)                :: numuse       ! NUMber of Upper Surface Elements
    INTEGER (kind=4)                :: numlse       ! NUMber of Lower Surface Elements
    ! ---------------------------------------------------------------------------------------------
    INTEGER (kind=4),   ALLOCATABLE :: buse(:)      ! Blade Upper Surface Elements (external global labels)
    INTEGER (kind=4),   ALLOCATABLE :: blse(:)      ! Blade Lower Surface Elements (external global labels)
    ! ---------------------------------------------------------------------------------------------
    TYPE(bed),          POINTER     :: next         ! Next set in the list
  ENDTYPE bed
  
  CHARACTER(len=10)                 :: bEType                 ! Blades Elements TYPE - used to define pair%eType for blades
  
  INTEGER                           :: n = 0,               & ! aero. step counter
                                       nest,                & ! current struct. step inside an aero. step - between tow aero. calculations
                                       numnest                ! number of struct. steps for present aero. step
  INTEGER( kind = 4 )               :: hub_master             ! Master node internal label
  
  REAL( kind = 8 )                  :: CF2F                   ! Force coefficient
  REAL( kind = 8 ), ALLOCATABLE     :: loa(:,:)               ! Loads from aero. step for adding to resid(ndofn,npoin) in npo_db
                                                              ! allocated in NLUVLM_Calc.f90 (if n==0)
                                                              
  INTEGER                           :: nIrrDCP = 0            ! number of Control Points with IRRegular Dependency
  REAL(8),          ALLOCATABLE     :: nbstH(:,:)             ! NBST velocity transfering factors (for aero control points
                                                              ! with irregular dependency (if shell structural model is used)
                                                              ! allocated in inter_ini/inter_sr
  
  TYPE(pair),       POINTER         :: headpair, tailpair
  TYPE(bed),        POINTER         :: headbed, tailbed
  
CONTAINS ! =============================================
  
  ! ---------- Subroutines to deal with pairs ----------
  ! ====================================================
  
  SUBROUTINE ini_pairs (headpair, tailpair)
    ! Initializes a list of pairs
    
    IMPLICIT NONE
    
    TYPE(pair), POINTER :: headpair, tailpair
    
    NULLIFY (headpair, tailpair)       !initializes first and last pointer

  ENDSUBROUTINE ini_pairs
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE add_pair (new, head, tail)
    !This subroutine adds data to the end of the list
    
    IMPLICIT NONE
    
    !Dummy arguments
    TYPE (pair), POINTER :: new, head, tail
    
    !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN       !list is empty, start it
      head => new                           !first element
      tail => new                           !last element
      NULLIFY (tail%next)                   !last poit to nothing
    
    ELSE                                    !add a segment to the list
      tail%next => new                      !point to the new element
      NULLIFY (new%next)                    !nothing beyond the last
      tail => new                           !new last element
    
    ENDIF
  ENDSUBROUTINE add_pair
  
! ---------------------------------------------------------------------------
  
  ! ----------- Subroutines to deal with BED -----------
  ! ====================================================
  
  SUBROUTINE ini_bed (headbed, tailbed)
    ! Initializes a list of bed
    
    IMPLICIT NONE
    
    TYPE(bed), POINTER :: headbed, tailbed
    
    NULLIFY (headbed, tailbed)       !initializes first and last pointer

  ENDSUBROUTINE ini_bed
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE add_bed (new, head, tail)
    !This subroutine adds data to the end of the list
    
    IMPLICIT NONE
    
    !Dummy arguments
    TYPE (bed), POINTER :: new, head, tail
    
    !Check if a list is empty
    IF (.NOT. ASSOCIATED (head)) THEN       !list is empty, start it
      head => new                           !first element
      tail => new                           !last element
      NULLIFY (tail%next)                   !last poit to nothing
    
    ELSE                                    !add a segment to the list
      tail%next => new                      !point to the new element
      NULLIFY (new%next)                    !nothing beyond the last
      tail => new                           !new last element
    
    ENDIF
  ENDSUBROUTINE add_bed
  
! ---------------------------------------------------------------------------
    
  SUBROUTINE findBed(posic)
  
    IMPLICIT NONE
    
    TYPE(pair),   POINTER ::  posic
    TYPE(bed),    POINTER ::  current
    
    
    IF (ASSOCIATED (headbed)) THEN        ! Check if a list is empty
      current => headbed
      DO ! Loop over bed
    
        IF(TRIM(posic%m_est) == TRIM(current%m_est) )THEN
          posic%bed => current
          exit
        ENDIF
        
        IF (ASSOCIATED (current%next)) THEN    ! Check if it's not the end of the list
          current => current%next
        ELSE
          exit
        ENDIF
        
      ENDDO
    ENDIF
    
  ENDSUBROUTINE findBed
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_inter
    
    ! DUMPINg_INTERaction data
    
    ! Mauro S. Maza - 17/12/2012
    !                 07/01/2016
    
    ! Dumps interaction data in restart file
    
    USE ctrl_db,      ONLY: ndofn, npoin
    
    IMPLICIT NONE
    
    TYPE(pair),     POINTER :: posic
    TYPE(bed),      POINTER :: posicBed
    INTEGER                 :: j, k, lBou, uBou
    
    
    ! Variables in inter_db (inter_db.f90)
    WRITE(50,ERR=9999)  bEType,                             &
                        n, nest, numnest,                   &
                        hub_master, CF2F
    WRITE(50,ERR=9999)  ((loa(j,k),j=1,ndofn),k=1,npoin)
    WRITE(50,ERR=9999)  nIrrDCP
    WRITE(50,ERR=9999)  ((nbstH(j,k),j=1,27),k=1,nIrrDCP) ! needs to be allocated on restart
    
    
    ! User difined variables in inter_db (inter_db.f90)
    ! BED - Blade Elements Distribution
    IF(TRIM(bEType)=='NBST')THEN ! shell model for blades
      IF (ASSOCIATED (headbed)) THEN   ! Check if a list is empty
        posicBed => headbed
        DO ! Loop over pairs
          ! NO POINTERS SAVED
          WRITE(50,ERR=9999) posicBed%m_est, posicBed%numuse, posicBed%numlse
          WRITE(50,ERR=9999) (posicBed%buse(j),j=1,posicBed%numuse)
          WRITE(50,ERR=9999) (posicBed%blse(j),j=1,posicBed%numlse)
          
          
          IF (ASSOCIATED (posicBed%next)) THEN ! Check if it's not the end of the list
            posicBed => posicBed%next
          ELSE
            exit
          ENDIF
          
        ENDDO
      ENDIF
    ENDIF
    
    
    ! PAIRS
    IF (ASSOCIATED (headpair)) THEN ! Check if a list is empty
      posic => headpair
      DO ! Loop over pairs
        ! NO POINTERS SAVED
        ! data not depending on structural element type
        WRITE(50,ERR=9999) posic%m_est, posic%m_aer, posic%eType, posic%numbnodes
        WRITE(50,ERR=9999) (posic%nodes(j),j=1,posic%numbnodes)
        WRITE(50,ERR=9999) posic%nelem, posic%nan, posic%np
        WRITE(50,ERR=9999) posic%eType
        
        ! data depending on structural element type
        SELECT CASE (TRIM(posic%eType))
        CASE('RIGID') ! ---------------------
          WRITE(50,ERR=9999) posic%mnods(1,1)
          WRITE(50,ERR=9999) ((posic%dists(j,k),j=1,3),k=1,posic%nan)
          WRITE(50,ERR=9999) posic%cpmnods(1,1)
          WRITE(50,ERR=9999) ((posic%cpdists(j,k),j=1,3),k=1,posic%np)
          
          
        CASE('BEAM') ! ---------------------
          WRITE(50,ERR=9999) ((posic%lnods(j,k),j=1,2),k=1,posic%nelem)
          WRITE(50,ERR=9999) (posic%xita(1,j),j=1,posic%nan)
          WRITE(50,ERR=9999) ((posic%mnods(j,k),j=1,2),k=1,posic%nan)
          WRITE(50,ERR=9999) ((posic%dists(j,k),j=1,6),k=1,posic%nan)
          WRITE(50,ERR=9999) (posic%cpxita(1,j),j=1,posic%np)
          WRITE(50,ERR=9999) ((posic%cpmnods(j,k),j=1,2),k=1,posic%np)
          WRITE(50,ERR=9999) ((posic%cpdists(j,k),j=1,6),k=1,posic%np)
          
          
        CASE('NBST') ! --------------------- 
          lBou = LBOUND(posic%elmIndxs,1) ! necessary for restart
          uBou = UBOUND(posic%elmIndxs,1) ! necessary for restart
          WRITE(50,ERR=9999) lBou, uBou   ! necessary for restart
          WRITE(50,ERR=9999) ((posic%lnods(j,k),j=1,3),k=1,posic%nelem)
          WRITE(50,ERR=9999) (posic%elmIndxs(j),j=lBou,uBou)
          WRITE(50,ERR=9999) ((posic%xita(j,k),j=1,6),k=1,posic%nan)
          WRITE(50,ERR=9999) ((posic%mnods(j,k),j=1,6),k=1,posic%nan)
          WRITE(50,ERR=9999) ((posic%dists(j,k),j=1,2),k=1,posic%nan)
          WRITE(50,ERR=9999) (posic%reguDep(j),j=1,posic%nan)
          WRITE(50,ERR=9999) ((posic%cpxita(j,k),j=1,6),k=1,posic%np)
          WRITE(50,ERR=9999) ((posic%cpmnods(j,k),j=1,6),k=1,posic%np)
          WRITE(50,ERR=9999) ((posic%cpdists(j,k),j=1,2),k=1,posic%np)
          WRITE(50,ERR=9999) (posic%cpReguDep(j),j=1,posic%np)
          
          
        ENDSELECT
        
        IF (ASSOCIATED (posic%next)) THEN   ! Check if it's not the end of the list
          posic => posic%next
        ELSE
          exit
        ENDIF
        
      ENDDO
    ENDIF
    
    
    RETURN
    
    9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_inter
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_inter
    
    ! RESTARt_INTERaction data
    
    ! Mauro S. Maza - 17/12/2012
    !                 07/01/2016
    
    ! Reads interaction data from restart file
    
    USE ctrl_db,      ONLY: ndofn, npoin
    
    IMPLICIT NONE
    
    TYPE(pair),         POINTER     ::  posic
    TYPE(bed),          POINTER     ::  posicBed
    INTEGER                         ::  i, j, k, lBou, uBou
    
    
    ! Variables in inter_db (inter_db.f90)
    READ (51)           bEType,                             &
                        n, nest, numnest,                   &
                        hub_master, CF2F
    ALLOCATE( loa(ndofn,npoin))
    READ (51)          ((loa(j,k),j=1,ndofn),k=1,npoin)
    READ (51)           nIrrDCP
    ALLOCATE( nbstH(27,nIrrDCP))
    READ (51)          ((nbstH(j,k),j=1,27),k=1,nIrrDCP)
    
    
    ! User difined variables in inter_db (inter_db.f90)
    ! BED - Blade Elements Distribution
    IF(TRIM(bEType)=='NBST')THEN ! shell model for blades
      CALL ini_bed(headbed, tailbed)
      DO  i=1,3 ! Loop over 3 bed
        ALLOCATE(posicBed) ! allocate memory
        CALL add_bed(posicBed, headbed, tailbed)
        READ (51)          posicBed%m_est, posicBed%numuse, posicBed%numlse
        ALLOCATE(posicBed%buse(posicBed%numuse), &
                 posicBed%blse(posicBed%numlse))
        READ (51)          (posicBed%buse(j),j=1,posicBed%numuse)
        READ (51)          (posicBed%blse(j),j=1,posicBed%numlse)
      ENDDO
    ENDIF
    
    
    ! PAIRS
    CALL ini_pairs (headpair, tailpair)
    DO i=1,6 ! Loop over 6 pairs
      ALLOCATE(posic)
      CALL add_pair (posic, headpair, tailpair)
      ! data not depending on structural element type
      READ (51)          posic%m_est, posic%m_aer, posic%eType, posic%numbnodes
      ALLOCATE( posic%nodes(posic%numbnodes) ) 
      READ (51)          (posic%nodes(j),j=1,posic%numbnodes)
      READ (51)          posic%nelem, posic%nan, posic%np
      READ (51)          posic%eType
      
      ! NO POINTERS SAVED
      ! Aero. nodes data
      SELECT CASE (TRIM(posic%m_aer))
      CASE ('blade1')
        posic%nodsaer => blade1nodes    ! Array with nodal info of posic%m_aer
        posic%panel   => blade1panels   ! Array with panels' info of posic%m_aer
        posic%section => blade1sections ! Array with sections info of blade1 lifting surface
      CASE ('blade2')
        posic%nodsaer => blade2nodes    ! Array with nodal info of posic%m_aer
        posic%panel   => blade2panels   ! Array with panels' info of posic%m_aer
        posic%section => blade2sections ! Array with sections info of blade2 lifting surface
      CASE ('blade3')
        posic%nodsaer => blade3nodes    ! Array with nodal info of posic%m_aer
        posic%panel   => blade3panels   ! Array with panels' info of posic%m_aer
        posic%section => blade3sections ! Array with sections info of blade3 lifting surface
      CASE ('hub')
        posic%nodsaer => hubnodes   ! Array with nodal info of posic%m_aer
        posic%panel   => hubpanels  ! Array with panels' info of posic%m_aer
      CASE ('nacelle')
        posic%nodsaer => nacellenodes   ! Array with nodal info of posic%m_aer
        posic%panel   => nacellepanels  ! Array with panels' info of posic%m_aer
      CASE ('tower')
        posic%nodsaer => towernodes     ! Array with nodal info of posic%m_aer
        posic%panel   => towerpanels    ! Array with panels' info of posic%m_aer
      ENDSELECT
      
      ! data depending on structural element type
      SELECT CASE (TRIM(posic%eType))
      CASE('RIGID') ! ---------------------
        ALLOCATE( posic%mnods(1,1) )
        ALLOCATE( posic%dists(3,posic%nan) )
        ALLOCATE( posic%cpmnods(1,1) )
        ALLOCATE( posic%cpdists(3,posic%np) )
        READ (51)          posic%mnods(1,1)
        READ (51)          ((posic%dists(j,k),j=1,3),k=1,posic%nan)
        READ (51)          posic%cpmnods(1,1)
        READ (51)          ((posic%cpdists(j,k),j=1,3),k=1,posic%np)
          
          
      CASE('BEAM') ! ---------------------
        ALLOCATE( posic%lnods(2,posic%nelem) )
        ALLOCATE( posic%xita(1,posic%nan) )
        ALLOCATE( posic%mnods(2,posic%nan) )
        ALLOCATE( posic%dists(6,posic%nan) )
        ALLOCATE( posic%cpxita(1,posic%np) )
        ALLOCATE( posic%cpmnods(2,posic%np) )
        ALLOCATE( posic%cpdists(6,posic%np) )
        READ (51)          ((posic%lnods(j,k),j=1,2),k=1,posic%nelem)
        READ (51)          (posic%xita(1,j),j=1,posic%nan)
        READ (51)          ((posic%mnods(j,k),j=1,2),k=1,posic%nan)
        READ (51)          ((posic%dists(j,k),j=1,6),k=1,posic%nan)
        READ (51)          (posic%cpxita(1,j),j=1,posic%np)
        READ (51)          ((posic%cpmnods(j,k),j=1,2),k=1,posic%np)
        READ (51)          ((posic%cpdists(j,k),j=1,6),k=1,posic%np)
          
          
      CASE('NBST') ! --------------------- 
        READ (51)          lBou, uBou
        ALLOCATE( posic%lnods(3,posic%nelem) )
        ALLOCATE( posic%elmIndxs(lBou:uBou) )
        ALLOCATE( posic%xita(6,posic%nan) )
        ALLOCATE( posic%mnods(6,posic%nan) )
        ALLOCATE( posic%dists(2,posic%nan) )
        ALLOCATE( posic%reguDep(posic%nan) )
        ALLOCATE( posic%cpxita(6,posic%np) )
        ALLOCATE( posic%cpmnods(6,posic%np) )
        ALLOCATE( posic%cpdists(2,posic%np) )
        ALLOCATE( posic%cpReguDep(posic%np) )
        READ (51)          ((posic%lnods(j,k),j=1,3),k=1,posic%nelem)
        READ (51)          (posic%elmIndxs(j),j=lBou,uBou)
        READ (51)          ((posic%xita(j,k),j=1,6),k=1,posic%nan)
        READ (51)          ((posic%mnods(j,k),j=1,6),k=1,posic%nan)
        READ (51)          ((posic%dists(j,k),j=1,2),k=1,posic%nan)
        READ (51)          (posic%reguDep(j),j=1,posic%nan)
        READ (51)          ((posic%cpxita(j,k),j=1,6),k=1,posic%np)
        READ (51)          ((posic%cpmnods(j,k),j=1,6),k=1,posic%np)
        READ (51)          ((posic%cpdists(j,k),j=1,2),k=1,posic%np)
        READ (51)          (posic%cpReguDep(j),j=1,posic%np)
        CALL findBed(posic)             ! Blade Elements Distribution
        
        
      ENDSELECT
      
    ENDDO
    
    
  ENDSUBROUTINE restar_inter
  
  
ENDMODULE inter_db

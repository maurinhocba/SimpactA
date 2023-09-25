
MODULE mastrnods_sr
  
  ! Mauro S. Maza - 15/09/2015
  
  USE mastrnods_db
  USE inter_db,         ONLY:   pair
  USE npo_db,           ONLY:   coord, eule0
  
  IMPLICIT NONE
  
  INCLUDE 'mkl_blas.fi' ! INTERFACE blocks for BLAS/MKL routines
  
  TYPE(pair),   POINTER :: posic

CONTAINS !===================================================================
  
        ! SUBROUTINES
        ! mngMstrs:             manage master nodes determination
        
        ! robMstrs:             RIGID or BEAM elements structural model
        !   proj_ele08:         determine projection BEAM element
        
        ! nbstMstrs:            NBST elements structural model
        !   nbstConnects:       elements' connectivities
        !   nbstMidpoints:      subset midpoints
        !   nbstIntersect:      intersection point of a line and an element plane
        !   nbstLCoords:        local coordinates of a point laying in the plane of an NBST element
        !   nbstIrrDep2:        irregular dependency (for aero nodes out of structural mesh) - algorithm 2
        !   nbstMstrsTecplot:   debug masters determination
        
        ! distancies:           dists. between target point and set of points
  
  !--------------------------------------------------------------------------
  
  SUBROUTINE mngMstrs( posic )
    
    ! Manages master nodes determination
    ! depending on element type for structural model
    
    IMPLICIT NONE
    
    TYPE(pair), POINTER :: posic
    
    
    SELECT CASE (TRIM(posic%eType))
    CASE('RIGID')
      CALL robMstrs( posic )
      
    CASE('BEAM')
      CALL robMstrs( posic )
      
    CASE('NBST')
      CALL nbstMstrs( posic )
    
    ENDSELECT
    
  ENDSUBROUTINE mngMstrs
  
  ! -------------------------------------------------------------------------
  
  ! --------------- RIGID and BEAM elements ---------------
  !       old and unupdated, but (apparently) working
  ! =======================================================
  
  SUBROUTINE robMstrs( posic )
    
    ! Rigid Or Beam elements MaSTeRS
    
    USE inter_db,       ONLY:   hub_master
    USE surf_db,        ONLY:   surfa,                    &    ! variables
                                new_srf, store_segs            ! subroutines
    
    IMPLICIT NONE
    
    ! dummy args
    TYPE(pair), POINTER :: posic
    
    ! internal vars
	LOGICAL							:: found
	INTEGER (kind=4) 				:: tipo, n1, n2
	INTEGER 						:: i
	REAL(kind = 8), 	POINTER		:: midpoints(:,:)
	REAL(kind = 8),		ALLOCATABLE :: xyz_aux(:,:)
	REAL(kind = 8) 					:: r_an(3), d1(3), d1l(3), &
									   lambda1(3,3), lamVec1(9)
    
    
	INTERFACE
      INCLUDE 'elemnt.h'
      INCLUDE 'inrotm.h'
	END INTERFACE
    
    
    ! Structural elements data
    CALL new_srf(surfa) ! Allocate memory for data
    CALL elemnt ('SURFAC', name=posic%m_est, flag2=found, flag1=.TRUE., igrav=tipo) ! Get data
                                                                                    ! flag1=.TRUE. is to find master node of rigid body
    IF (.NOT.found) THEN
      CALL runend('RDSURF:ELEMENT SET NOT FOUND       ')
    ENDIF
    
    SELECT CASE (tipo)
    CASE (8)    !  'BEAM ' elements
      posic%nelem = surfa%nelem    ! Number of elements in the set
      ALLOCATE( posic%lnods(2,posic%nelem) )  ! Conectivities - internal labels
      CALL store_segs(surfa%head,posic%lnods,2,posic%nelem)   ! This assumes that the elements has only 2 nodes
      CALL nodeslabs( posic )
      ! Midpoints of the elements in the set
      ALLOCATE( midpoints(3,posic%nelem) )    ! Elements midpoints array
      DO i=1,posic%nelem
        n1 = posic%lnods(1,i)
        n2 = posic%lnods(2,i)
        midpoints(1:3,i) = 0.5 * ( coord(:,n1) + coord(:,n2) )
      ENDDO
                  
      ! Determine projection element (actually, master struct. nodes),
      ! coordinate for projection and initial distance vectors
      ! Aero. nodes
      ALLOCATE( xyz_aux(posic%nan,3) )
      DO i=1,posic%nan
          xyz_aux(i,1:3) = posic%nodsaer(i)%xyz(:) ! aero. nodes coords. are not stored all together in a single array
      ENDDO
      ALLOCATE( posic%xita(1,posic%nan)       , &
                posic%mnods(2,posic%nan)    , &
                posic%dists(6,posic%nan)    )
      CALL proj_ele08(  posic%xita, posic%mnods, posic%dists, &
                        posic%nelem, posic%lnods, midpoints, posic%nan, xyz_aux )
      DEALLOCATE( xyz_aux )
      ! Control Points
      ALLOCATE( xyz_aux(posic%np,3) )
      DO i=1,posic%np
        xyz_aux(i,1:3) = posic%panel(i)%xyzcp(:)
      ENDDO
      ALLOCATE( posic%cpxita(1,posic%np)       , &
                posic%cpmnods(2,posic%np)    , &
                posic%cpdists(6,posic%np)    )
      CALL proj_ele08(  posic%cpxita, posic%cpmnods, posic%cpdists, &
                        posic%nelem, posic%lnods, midpoints, posic%np, xyz_aux )
      DEALLOCATE( xyz_aux, midpoints ) !posic%lnods,
    
    CASE (10)    !  'RIGID' elements
      n1 = surfa%nelem        ! Master node internal label
      CALL inrotm(eule0(1:3,n1),lamVec1) ! elements of rotation matrix from definition euler angles
      lambda1 = RESHAPE(lamVec1,(/3,3/))    ! Rotation matrix for n1 local system
      ALLOCATE( posic%mnods(1,1) , posic%cpmnods(1,1) )
      posic%mnods(1,1) = n1
      posic%cpmnods(1,1) = n1
      IF( TRIM(posic%m_aer) == 'hub' )THEN
        hub_master = n1
      ENDIF
                  
      ! Determine initial distance vectors
      ! Aero. nodes
      ALLOCATE( posic%dists(3,posic%nan) )
      DO i=1,posic%nan    ! Loop over aero. nodes
        r_an = posic%nodsaer(i)%xyz(:)    ! Position of aero. node
        d1 = r_an - coord(:,n1)            ! Distance vector from master node to aero. node n1-->an
        d1l = MATMUL(d1,lambda1)        ! Initial distance vector in local coords. of master node n1
        posic%dists(1:3,i) = d1l
      ENDDO
      ! Control Points
      ALLOCATE( posic%cpdists(3,posic%np) )
      DO i=1,posic%np    ! Loop over Control Points
        r_an = posic%panel(i)%xyzcp(:)    ! Position of Control Point
        d1 = r_an - coord(:,n1)            ! Distance vector from master node to Control Point n1-->CP
        d1l = MATMUL(d1,lambda1)        ! Initial distance vector in local coords. of master node n1
        posic%cpdists(1:3,i) = d1l
      ENDDO
    
    ENDSELECT
    
  CONTAINS !...............................................
  
      SUBROUTINE nodeslabs( posic )
        
        ! Mauro S. Maza - 17/11/2011
        
        ! Determines the labels of the nodes present in the set
        ! ASSUMES that the elements are of type 'BEAM ' (two nodes)
        
        IMPLICIT NONE
        
        TYPE(pair), 		    POINTER		:: posic
        INTEGER (kind = 4 ),    ALLOCATABLE :: nodes_aux(:)
        INTEGER                             :: i, j, k, try
        LOGICAL                             :: flag
        
        ALLOCATE( nodes_aux( 2*posic%nelem ) )  ! auxiliar variable
        nodes_aux(1) = posic%lnods(1,1)         ! first label in the list
        k = 1   ! counter for number of different labels found
        DO i=1,posic%nelem    ! loop over elements (posic%lnods columns)
          ! labels of first node of each element (first row of posic%lnods)
          try = posic%lnods(1,i)
          flag = .TRUE.
          DO j=1,k  ! loop over labels already found
            IF( try == nodes_aux(j) )THEN
              flag = .FALSE.
            ENDIF
          ENDDO
          IF( flag )THEN
            k = k +1
            nodes_aux(k) = try
          ENDIF
          ! labels of second node of each element (second row of posic%lnods)
          try = posic%lnods(2,i)
          flag = .TRUE.
          DO j=1,k  ! loop over labels already found
            IF( try == nodes_aux(j) )THEN
              flag = .FALSE.
            ENDIF
          ENDDO
          IF( flag )THEN
            k = k +1
            nodes_aux(k) = try
          ENDIF
        ENDDO
        ALLOCATE( posic%nodes(k) )
        posic%nodes = nodes_aux( 1:k )
        posic%numbnodes = k
        DEALLOCATE( nodes_aux )
      
      ENDSUBROUTINE nodeslabs
  
! ---------------------------------------------------------------------------
    
      SUBROUTINE proj_ele08( xita, mnods, dists, nelem, lnods, midpoints, npoints, xyz ) 
                             ! out ------------  ! in --------------------------------
        
        ! Determine projection element (actually, master struct. nodes),
        ! coordinate for projection and initial distance vectors
        
        IMPLICIT NONE
            
        INTEGER (kind=4), ALLOCATABLE   ::  lnods(:,:)
        INTEGER (kind=4)                ::  nelem, n1, n2
        INTEGER,          ALLOCATABLE   ::  proele(:), mnods(:,:)
        INTEGER                         ::  i, j, befor, after, loc, npoints
        
        REAL,             ALLOCATABLE   ::  xita(:,:), dists(:,:)
        REAL(kind = 8),   POINTER       ::  midpoints(:,:)
        REAL(kind = 8),   ALLOCATABLE   ::  delta_aux(:), xyz(:,:)
        REAL(kind = 8)                  ::  length, xita1, &
                                            r_an(3), delta_xyz(3), vj(3), d1(3), d2(3), d1l(3), d2l(3), &
                                            lambda1(3,3), lambda2(3,3), lamVec1(9), lamVec2(9)
        
        INTERFACE
          INCLUDE 'inrotm.h'
        END INTERFACE
        
        ALLOCATE( proele(npoints) , delta_aux(nelem) )  ! Auxiliar variable
        
        DO i=1,npoints    ! Loop over points
          r_an = xyz(i,:)    ! Position of point
          ! Search nearest midpoint for each point
          ! Calculate distancies
          DO j=1,nelem    ! Loop over elements' midpoints
            delta_xyz = r_an - midpoints(:,j)
            !delta_aux(j) = SQRT( DOT_PRODUCT(delta_xyz,delta_xyz) )
            delta_aux(j) = dnrm2(3,delta_xyz,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
          ENDDO
          proele(i) = MINLOC(delta_aux,1)    ! Find nearest element to the ith aero. node
                                             ! First guess
          befor = 0    ! For nodes in dead zones
          search1:DO    ! Until \xita_1 \in [0,1]
            ! All this assumes that every node is sheared ONLY by two elements
            n1 = lnods(1,proele(i))    ! Element first node
            n2 = lnods(2,proele(i))    ! Element second node
            d1 = r_an - coord(:,n1)    ! Distance vector from first  node to aero. node n1-->an
            d2 = r_an - coord(:,n2)    ! Distance vector from second node to aero. node n2-->an
            CALL inrotm(eule0(1:3,n1),lamVec1) ! elements of rotation matrix from definition euler angles
            CALL inrotm(eule0(1:3,n2),lamVec2) ! elements of rotation matrix from definition euler angles
            lambda1 = RESHAPE(lamVec1,(/3,3/))    ! Rotation matrix for n1 local system
            lambda2 = RESHAPE(lamVec2,(/3,3/))    ! Rotation matrix for n2 local system
            d1l = MATMUL(d1,lambda1)    ! Initial distance vector in local coords. of master node n1
            d2l = MATMUL(d2,lambda2)    ! Initial distance vector in local coords. of master node n2
            vj = coord(:,n2) - coord(:,n1)        ! Calculate tangent vector to master segment, n1-->n2
            !length = SQRT(DOT_PRODUCT(vj,vj))    ! Segment length
            length = dnrm2(3,vj,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            vj = vj/length                        ! Versor
            xita1 = DOT_PRODUCT(d1,vj)/length    ! Local coordinate [0,1] means projects in
            IF (xita1 < 0d0) THEN    ! Find neighbour (sharing node n1)
              loc = 0
              search2:DO j=1,nelem  ! determines which element is the neighbour
                IF (lnods(1,j)==n1 .OR. lnods(2,j)==n1 .AND. j/=proele(i)) THEN
                  loc = j
                  exit search2
                ENDIF
              ENDDO search2
              IF (loc/=0) THEN    ! There exists another element (not proele(i)) sharing node n1
                after = loc
                IF( after==befor )THEN ! aero node in dead zone
                  xita(1,i) = 1
                  mnods( 1 ,i) = n1
                  mnods( 2 ,i) = n1
                  dists(1:3,i) = d1l
                  dists(4:6,i) = d1l
                  exit search1
                ELSE
                  befor = proele(i)
                  proele(i) = after
                ENDIF
              ELSE ! no other element sharing n1 - aero node out of str. mesh range
                xita(1,i) = 1
                mnods( 1 ,i) = n1
                mnods( 2 ,i) = n1
                dists(1:3,i) = d1l
                dists(4:6,i) = d1l
                exit search1
              ENDIF
            ELSEIF (xita1 > 1d0) THEN    ! Find neighbour (sharing node n2)
              loc = 0
              search3:DO j=1,nelem  ! determines which element is the neighbour
                IF (lnods(1,j)==n2 .OR. lnods(2,j)==n2 .AND. j/=proele(i)) THEN
                  loc = j
                  exit search3
                ENDIF
              ENDDO search3
              IF (loc/=0) THEN    ! There exists another element (not proele(i)) shearing node n2
                after = loc
                IF (after==befor) THEN ! aero node in dead zone
                  xita(1,i) = 1
                  mnods( 1 ,i) = n2
                  mnods( 2 ,i) = n2
                  dists(1:3,i) = d2l
                  dists(4:6,i) = d2l
                  exit search1
                ELSE
                  befor = proele(i)
                  proele(i) = after
                ENDIF
              ELSE ! no other element sharing n2 - aero node out of str. mesh range
                xita(1,i) = 1
                mnods( 1 ,i) = n2
                mnods( 2 ,i) = n2
                dists(1:3,i) = d2l
                dists(4:6,i) = d2l
                exit search1
              ENDIF
            ELSE ! projects - it may project in other elements
              xita(1,i) = xita1
              mnods( 1 ,i) = n1
              mnods( 2 ,i) = n2
              dists(1:3,i) = d1l
              dists(4:6,i) = d2l
              exit search1
            ENDIF
          ENDDO search1
        ENDDO
        
        DEALLOCATE( proele , delta_aux )    ! Auxiliar variable
        
      ENDSUBROUTINE proj_ele08
    
  ENDSUBROUTINE robMstrs
  
  !--------------------------------------------------------------------------
  
  ! -------------------- NBST elements --------------------
  ! =======================================================
  
  SUBROUTINE nbstMstrs( posic )
    
    ! ABSTRACT
    ! 1. initial computations
    !   1.1 find structural elements set
    !   1.2 set elements' connectivities
    !   1.3 set elements' midpoints
    !   1.4 versor parallel to blade twist axis
    ! 2. loop over aero nodes - NOT in one shot
    !   (because of the way nodes are numbered)
    !   2.1 lifting surface
    !   | 2.1.1 loop over lines between sections (as defined by Cristian G.)
    !   | | 2.1.1.1 thickness vector (from LE to TE)
    !   | | 2.1.1.1 loop over nodes in the line
    !   | |   2.1.1.1.1 upper surface
    !   | |     2.1.1.1.1.1 distance from midpoints
    !   | |     2.1.1.1.1.2 loop over 10 best chances
    !   | |       2.1.1.1.2.1 element plane
    !   | |       2.1.1.1.2.2 intersection of plane with thikness vector
    !   | |       2.1.1.1.2.3 projects?
    !   | |       2.1.1.1.2.4 YES => save local coordinates and distance from plane
    !   | |       2.1.1.1.2.5 NO => grid out of mesh, determine position with local system and save coords and flag, print message
    !   | |   2.1.1.1.2 lower surface
    !   | |             SAME FOR upper surface
    !   | 2.1.2 loop over sections (for the control points)
    !   |   2.1.2.1 loop over CPs in the section
    !   |           SAME FOR nodes
    !   |           (determine thickness vector from nodal points at each side)
    !   |     2.1.2.1.1 upper surface
    !   |     2.1.2.1.2 lower surface
    !   2.2 root-transition
    !       SAME FOR lifting surface
    !       there are more nodes and CPs for each section
    !     2.2.1 loop over lines between sections (analogous to those defined by Cristian G.)
    !         2.2.1.1.1 upper surface
    !         2.2.1.1.2 lower surface
    !     2.2.2 loop over sections (for the control points)
    !         2.2.2.1.1 upper surface
    !         2.2.2.1.2 lower surface
    ! 3. final computations
    !   3.1 normalize distancies from elements
    !   3.2 deallocate internal arrays
    !   3.3 determine number of control points with irregular dependency and data order in nbstH array
    
    ! speed may be improved not looking for lower surface projection when the
    ! upper one has failed (projectsU=.FALSE.) - but this will cause the code
    ! to become less readable, and it is already fast enough
    
    USE ele13_db,       ONLY:   srch_ele13,         & ! routines
                                ele13_set, ele13,   & ! types
                                head                          ! variables
    USE introut_db,     ONLY:   msg
    USE introut_sr,     ONLY:   introut_mng
    USE inter_db,       ONLY:   nIrrDCP
    USE StructuresA,    ONLY:   nsec, npan
    USE param_db,       ONLY:   mich
    USE ifport  !Intel Compiler - what for?
    
    IMPLICIT NONE
    
    ! dummy args
    TYPE(pair), POINTER :: posic
    ! Functions
    CHARACTER(len=mich):: inttoch
    
    ! internal vars
    LOGICAL                                 ::  exists, projectsU, projectsL
    INTEGER                                 ::  i, j, k, ind,                   &
                                                n1, n2, n3, n4, na, ns, eel,    &
                                                nnl, npl, nprt, nsrt
    REAL(8)                                 ::  d, tvn, chord(3), thick(3), tv(3),   &
                                                xa(3), x1(3), x2(3), x3(3), ip(3)
    REAL(8), ALLOCATABLE, DIMENSION(:,:)    ::  UpMidPo, LoMidPo
    REAL(8), ALLOCATABLE, DIMENSION(:)      ::  distsUp, distsLo, total
    LOGICAL                                 ::  found
    TYPE (ele13_set),   POINTER             ::  set, anterSet
    
    
    ! 1. INITIAL COMPUTATIONS
    ! find structural elements set
    CALL srch_ele13(head, anterSet, set, posic%m_est, found) ! outputs SET pointing to set named elmSetName=posic%m_est
                                                             ! in the list of NBST elements (whose first link is HEAD)
    ! set element connectivities
    ALLOCATE( posic%lnods(3,posic%nelem ) )
    posic%lnods = 0
    CALL nbstConnects(posic%m_est, posic%lnods, posic%elmIndxs)
    ! midpoints
    ALLOCATE( UpMidPo(3,posic%bed%numuse) ) ! midpoints of elements whose external labels are listed in posic%bed%buse
    ALLOCATE( LoMidPo(3,posic%bed%numlse) ) ! midpoints of elements whose external labels are listed in posic%bed%blse
    CALL nbstMidpoints(posic%bed%buse, posic%bed%numuse, posic%lnods(:,posic%elmIndxs(posic%bed%buse)), UpMidPo) ! upper surface
    CALL nbstMidpoints(posic%bed%blse, posic%bed%numlse, posic%lnods(:,posic%elmIndxs(posic%bed%blse)), LoMidPo) ! lower surface
    
    ! versor parallel to blade twist axis
    n1 = (nsec+1)*(npan+1)+1 ! root leading edge node - also first root node
    n2 = n1 + npan ! root trailing edge node
    n3 = n1 + FLOOR(REAL(npan)/2.0) ! root middle node in upper surface
    n4 = n2 + FLOOR(REAL(npan)/2.0) ! root middle node in lower surface
    chord = posic%nodsaer(n2)%xyz-posic%nodsaer(n1)%xyz
    thick = posic%nodsaer(n4)%xyz-posic%nodsaer(n3)%xyz
    CALL cross_product(tv ,chord ,thick) ! cross_product(out, in1, in2)
    tvn = dnrm2(3,tv,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
    tv = ( 1 / tvn ) * tv ! twist versor
    
    ALLOCATE( posic%xita (6,posic%nan)    , & ! aero nodes data
              posic%mnods(6,posic%nan)    , &
              posic%dists(2,posic%nan)    , &
              posic%reguDep(posic%nan)    )
    ALLOCATE( posic%cpxita (6,posic%np)   , & ! control points data
              posic%cpmnods(6,posic%np)   , &
              posic%cpdists(2,posic%np)   , &
              posic%cpReguDep(posic%np)   )
    ALLOCATE( distsUp(posic%bed%numuse), & ! distancies from aero node to midpoints of elements whose external labels are listed in posic%bed%buse
              distsLo(posic%bed%numlse)  ) ! distancies from aero node to midpoints of elements whose external labels are listed in posic%bed%blse
    
    
    ! 2. LOOP OVER AERO NODES - NOT IN ONE SHOT
    !   2.1 LIFTING SURFACE
    !     2.1.1 LOOP OVER LINES BETWEEN SECTIONS (AS DEFINED BY CRISTIAN G.)
    DO i=1,nsec+1
      
      ! thickness vector
      n1 = (i-1)*(npan+1)+1 ! leading edge node
      n2 =     i*(npan+1)   ! trailing edge node
      chord = posic%nodsaer(n2)%xyz-posic%nodsaer(n1)%xyz
      CALL cross_product(thick, chord, tv) ! cross_product(out, in1, in2)
      
      DO j=1,npan+1 ! loop over nodes in the line
        na = (i-1)*(npan+1)+j      ! aero node label
        xa = posic%nodsaer(na)%xyz ! aero node coords
        
    !         2.1.1.1.1 UPPER SURFACE
        ! loop over 150 best chances
        CALL distancies(distsUp, UpMidPo, xa) ! distancies between aero node and elems midpoints - upper surface
        d = -1.0D0
        projectsU = .FALSE.
        DO k=1,150
          d   = MINVAL( distsUp, 1, MASK = distsUp .GT. d ) ! distance to nearest element midpoint
          ind = MINLOC( distsUp, 1, MASK = distsUp .EQ. d ) ! index of nearest element midpoint in posic%bed%buse array
          ! find element, its nodes' coords and normal vector and determine intersection point
          ns = posic%bed%buse(ind) ! stru element label
          CALL nbstIntersect(ns,posic,xa,posic%mnods(1,na),posic%mnods(2,na),posic%mnods(3,na),x1,x2,x3,ip,exists)
          IF(exists)THEN
            posic%dists(1,na) = dnrm2(3,xa-ip,1) ! MKL library routine
                                                 ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            ! determine local coords
            CALL nbstLCoords(x1,x2,x3,ip,posic%xita(1:3,na),projectsU)
            
            IF(projectsU)THEN
              EXIT ! do-loop over distancies
            ENDIF
          ENDIF
        ENDDO
        
    !         2.1.1.1.2 LOWER SURFACE
        ! loop over 150 best chances
        CALL distancies(distsLo, LoMidPo, xa) ! distancies between aero node and elems midpoints - lower surface
        d = -1.0D0
        projectsL = .FALSE.
        DO k=1,150
          d   = MINVAL( distsLo, 1, MASK = distsLo .GT. d ) ! distance to nearest element midpoint
          ind = MINLOC( distsLo, 1, MASK = distsLo .EQ. d ) ! index of nearest element midpoint in posic%bed%buse array
          ! find element, its nodes' coords and normal vector and determine intersection point
          ns = posic%bed%blse(ind) ! stru element label
          CALL nbstIntersect(ns,posic,xa,posic%mnods(4,na),posic%mnods(5,na),posic%mnods(6,na),x1,x2,x3,ip,exists)
          IF(exists)THEN
            posic%dists(2,na) = dnrm2(3,xa-ip,1) ! MKL library routine
                                                 ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            ! determine local coords
            CALL nbstLCoords(x1,x2,x3,ip,posic%xita(4:6,na),projectsL)
            
            IF(projectsL)THEN
              EXIT ! do-loop over distancies
            ENDIF
          ENDIF
        ENDDO
        
        IF(projectsU .AND. projectsL)THEN
          posic%reguDep(na) = .TRUE.
        ELSE ! aero node out of near structural mesh range - take nearest struc. node as master
          posic%reguDep(na) = .FALSE.
          
          ! determine special dependency
          CALL nbstIrrDep2(xa,distsUp,distsLo, &
                           eel, &
                           posic%mnods(1,na),posic%mnods(2,na),posic%mnods(3,na), &
                           posic%xita(1:4,na) )
          
          ! nullify useless variables
          posic%mnods(4:6,na) = 0
          posic%xita(5:6,na) = 0
          
          ! write warning
          msg = 'No proj. elem. for aero. node '//TRIM(inttoch(na,6))//' in '//TRIM(posic%m_aer)//' grid - stru. elem '//TRIM(inttoch(eel,6))//' used as master'
          CALL introut_mng('repProbl')
        ENDIF
    
      ENDDO
    ENDDO ! end with grid nodes in lifting surface
    
    !     2.1.2 LOOP OVER SECTIONS (FOR THE CONTROL POINTS)
    DO i=1,nsec
      
      ! thickness vector
      n1 = (i-1)*(npan+1)+1 ! root-side leading edge node
      n2 =     i*(npan+1)   ! root-side trailing edge node
      n3 = n2+1        ! tip-side leading edge node
      n4 = n2+(npan+1) ! tip-side trailing edge node
      x1 = 0.5*( posic%nodsaer(n1)%xyz + posic%nodsaer(n3)%xyz) ! leading edge virtual node coords
      x2 = 0.5*( posic%nodsaer(n2)%xyz + posic%nodsaer(n4)%xyz) ! trailing edge virtual node coords
      chord = x2-x1
      CALL cross_product(thick, chord, tv) ! cross_product(out, in1, in2)
      
      DO j=1,npan ! loop over nodes in the section
        na = (i-1)*npan+j          ! control point (panel) label
        xa = posic%panel(na)%xyzcp ! control point coords
        
    !         2.1.2.1.1 UPPER SURFACE
        ! loop over 150 best chances
        CALL distancies(distsUp, UpMidPo, xa) ! distancies between aero node and elems midpoints - upper surface
        d = -1.0D0
        projectsU = .FALSE.
        DO k=1,150
          d   = MINVAL( distsUp, 1, MASK = distsUp .GT. d ) ! distance to nearest element midpoint
          ind = MINLOC( distsUp, 1, MASK = distsUp .EQ. d ) ! index of nearest element midpoint in posic%bed%buse array
          ! find element, its nodes' coords and normal vector and determine intersection point
          ns = posic%bed%buse(ind) ! stru element label
          CALL nbstIntersect(ns,posic,xa,posic%cpmnods(1,na),posic%cpmnods(2,na),posic%cpmnods(3,na),x1,x2,x3,ip,exists)
          IF(exists)THEN
            posic%cpdists(1,na) = dnrm2(3,xa-ip,1) ! MKL library routine
                                                   ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            ! determine local coords
            CALL nbstLCoords(x1,x2,x3,ip,posic%cpxita(1:3,na),projectsU)
            
            IF(projectsU)THEN
              EXIT ! do-loop over distancies
            ENDIF
          ENDIF
        ENDDO
        
    !         2.1.2.1.2 LOWER SURFACE
        ! loop over 150 best chances
        CALL distancies(distsLo, LoMidPo, xa) ! distancies between aero node and elems midpoints - lower surface
        d = -1.0D0
        projectsL = .FALSE.
        DO k=1,150
          d   = MINVAL( distsLo, 1, MASK = distsLo .GT. d ) ! distance to nearest element midpoint
          ind = MINLOC( distsLo, 1, MASK = distsLo .EQ. d ) ! index of nearest element midpoint in posic%bed%buse array
          ! find element, its nodes' coords and normal vector and determine intersection point
          ns = posic%bed%blse(ind) ! stru element label
          CALL nbstIntersect(ns,posic,xa,posic%cpmnods(4,na),posic%cpmnods(5,na),posic%cpmnods(6,na),x1,x2,x3,ip,exists)
          IF(exists)THEN
            posic%cpdists(2,na) = dnrm2(3,xa-ip,1) ! MKL library routine
                                                   ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            ! determine local coords
            CALL nbstLCoords(x1,x2,x3,ip,posic%cpxita(4:6,na),projectsL)
            
            IF(projectsL)THEN
              EXIT ! do-loop over distancies
            ENDIF
          ENDIF
        ENDDO
        
        IF(projectsU .AND. projectsL)THEN
          posic%cpReguDep(na) = .TRUE.
        ELSE ! aero node out of near structural mesh range - take nearest struc. node as master
          posic%cpReguDep(na) = .FALSE.
          
          ! determine special dependency
          CALL nbstIrrDep2(xa,distsUp,distsLo, &
                           eel, &
                           posic%cpmnods(1,na),posic%cpmnods(2,na),posic%cpmnods(3,na), &
                           posic%cpxita(1:4,na) )
          
          ! nullify useless variables
          posic%cpmnods(4:6,na) = 0
          posic%cpxita(5:6,na) = 0
          
          ! write warning
          msg = 'No proj. elem. for control point '//TRIM(inttoch(na,6))//' in '//TRIM(posic%m_aer)//' grid - stru. elem '//TRIM(inttoch(eel,6))//' used as master'
          CALL introut_mng('repProbl')
        ENDIF
    
      ENDDO
    ENDDO ! end with panel control points in lifting surface
    
    
    !   2.2 ROOT-TRANSITION
    ! here, the direction used to find masters might be other than that of the
    ! airfoil thickness (any diameter, for example)
    ! however, we still use the thickness direction because the aero grid
    ! converges from a circle at the root to a line at the beginning of the
    ! lifting-surface, and this is done wrt the thickness direction
    nnl = (nsec+1)*(npan+1) ! Number of Nodes in Lifting surface
    npl =  nsec   * npan    ! Number of Panels in Lifting surface
    nprt = 2*npan ! Number of Panels per section in Root-Transition
    nsrt = (posic%np - nsec*npan)/nprt ! Number of Sections in Root-Transition
    
    !     2.2.1 LOOP OVER LINES BETWEEN SECTIONS (ANALOGOUS TO THOSE DEFINED BY CRISTIAN G. FOR THE LIFTING SURFACE)
    DO i=1,nsrt+1 ! nodes in last line in root-transition (nsrt+1) are triplicated
      
      ! thickness vector
      n1 = nnl+(i-1)*(nprt+1)+1 ! leading edge node
      n2 = n1 + npan            ! trailing edge node
      chord = posic%nodsaer(n2)%xyz-posic%nodsaer(n1)%xyz
      CALL cross_product(thick, chord, tv) ! cross_product(out, in1, in2)
      
      DO j=1,nprt+1 ! loop over nodes in the line
        na = nnl+(i-1)*(nprt+1)+j  ! aero node label
        xa = posic%nodsaer(na)%xyz ! aero node coords
        
    !         2.2.1.1.1 UPPER SURFACE
        ! loop over 150 best chances
        CALL distancies(distsUp, UpMidPo, xa) ! distancies between aero node and elems midpoints - upper surface
        d = -1.0D0
        projectsU = .FALSE.
        DO k=1,150
          d   = MINVAL( distsUp, 1, MASK = distsUp .GT. d ) ! distance to nearest element midpoint
          ind = MINLOC( distsUp, 1, MASK = distsUp .EQ. d ) ! index of nearest element midpoint in posic%bed%buse array
          ! find element, its nodes' coords and normal vector and determine intersection point
          ns = posic%bed%buse(ind) ! stru element label
          CALL nbstIntersect(ns,posic,xa,posic%mnods(1,na),posic%mnods(2,na),posic%mnods(3,na),x1,x2,x3,ip,exists)
          IF(exists)THEN
            posic%dists(1,na) = dnrm2(3,xa-ip,1) ! MKL library routine
                                                 ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            ! determine local coords
            CALL nbstLCoords(x1,x2,x3,ip,posic%xita(1:3,na),projectsU)
            
            IF(projectsU)THEN
              EXIT ! do-loop over distancies
            ENDIF
          ENDIF
        ENDDO
        
    !         2.2.1.1.2 LOWER SURFACE
        ! loop over 150 best chances
        CALL distancies(distsLo, LoMidPo, xa) ! distancies between aero node and elems midpoints - lower surface
        d = -1.0D0
        projectsL = .FALSE.
        DO k=1,150
          d   = MINVAL( distsLo, 1, MASK = distsLo .GT. d ) ! distance to nearest element midpoint
          ind = MINLOC( distsLo, 1, MASK = distsLo .EQ. d ) ! index of nearest element midpoint in posic%bed%buse array
          ! find element, its nodes' coords and normal vector and determine intersection point
          ns = posic%bed%blse(ind) ! stru element label
          CALL nbstIntersect(ns,posic,xa,posic%mnods(4,na),posic%mnods(5,na),posic%mnods(6,na),x1,x2,x3,ip,exists)
          IF(exists)THEN
            posic%dists(2,na) = dnrm2(3,xa-ip,1) ! MKL library routine
                                                 ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            ! determine local coords
            CALL nbstLCoords(x1,x2,x3,ip,posic%xita(4:6,na),projectsL)
            
            IF(projectsL)THEN
              EXIT ! do-loop over distancies
            ENDIF
          ENDIF
        ENDDO
        
        IF(projectsU .AND. projectsL)THEN
          posic%reguDep(na) = .TRUE.
        ELSE ! aero node out of near structural mesh range - take nearest struc. node as master
          posic%reguDep(na) = .FALSE.
          
          ! determine special dependency
          CALL nbstIrrDep2(xa,distsUp,distsLo, &
                           eel, &
                           posic%mnods(1,na),posic%mnods(2,na),posic%mnods(3,na), &
                           posic%xita(1:4,na) )
          
          ! nullify useless variables
          posic%mnods(4:6,na) = 0
          posic%xita(5:6,na) = 0
          
          ! write warning
          msg = 'No proj. elem. for aero. node '//TRIM(inttoch(na,6))//' in '//TRIM(posic%m_aer)//' grid - stru. elem '//TRIM(inttoch(eel,6))//' used as master'
          CALL introut_mng('repProbl')
        ENDIF
    
      ENDDO
    ENDDO ! end with grid nodes in root-transition
    
    !     2.2.2 LOOP OVER SECTIONS (FOR THE CONTROL POINTS)
    DO i=1,nsrt
      
      ! thickness vector
      n1 = nnl+(i-1)*(nprt+1)+1 ! root-side leading edge node
      n2 = n1 + npan            ! root-side trailing edge node
      n3 = nnl+ i   *(nprt+1)+1 ! tip-side leading edge node
      n4 = n3 + npan            ! tip-side trailing edge node
      x1 = 0.5*( posic%nodsaer(n1)%xyz + posic%nodsaer(n3)%xyz) ! leading edge virtual node coords
      x2 = 0.5*( posic%nodsaer(n2)%xyz + posic%nodsaer(n4)%xyz) ! trailing edge virtual node coords
      chord = x2-x1
      CALL cross_product(thick, chord, tv) ! cross_product(out, in1, in2)
      
      DO j=1,nprt ! loop over nodes in the section
        na = npl+(i-1)*nprt+j      ! control point (panel) label
        xa = posic%panel(na)%xyzcp ! control point coords
        
    !         2.2.2.1.1 UPPER SURFACE
        ! loop over 150 best chances
        CALL distancies(distsUp, UpMidPo, xa) ! distancies between aero node and elems midpoints - upper surface
        d = -1.0D0
        projectsU = .FALSE.
        DO k=1,150
          d   = MINVAL( distsUp, 1, MASK = distsUp .GT. d ) ! distance to nearest element midpoint
          ind = MINLOC( distsUp, 1, MASK = distsUp .EQ. d ) ! index of nearest element midpoint in posic%bed%buse array
          ! find element, its nodes' coords and normal vector and determine intersection point
          ns = posic%bed%buse(ind) ! stru element label
          CALL nbstIntersect(ns,posic,xa,posic%cpmnods(1,na),posic%cpmnods(2,na),posic%cpmnods(3,na),x1,x2,x3,ip,exists)
          IF(exists)THEN
            posic%cpdists(1,na) = dnrm2(3,xa-ip,1) ! MKL library routine
                                                   ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            ! determine local coords
            CALL nbstLCoords(x1,x2,x3,ip,posic%cpxita(1:3,na),projectsU)
            
            IF(projectsU)THEN
              EXIT ! do-loop over distancies
            ENDIF
          ENDIF
        ENDDO
        
    !         2.2.2.1.2 LOWER SURFACE
        ! loop over 150 best chances
        CALL distancies(distsLo, LoMidPo, xa) ! distancies between aero node and elems midpoints - lower surface
        d = -1.0D0
        projectsL = .FALSE.
        DO k=1,150
          d   = MINVAL( distsLo, 1, MASK = distsLo .GT. d ) ! distance to nearest element midpoint
          ind = MINLOC( distsLo, 1, MASK = distsLo .EQ. d ) ! index of nearest element midpoint in posic%bed%buse array
          ! find element, its nodes' coords and normal vector and determine intersection point
          ns = posic%bed%blse(ind) ! stru element label
          CALL nbstIntersect(ns,posic,xa,posic%cpmnods(4,na),posic%cpmnods(5,na),posic%cpmnods(6,na),x1,x2,x3,ip,exists)
          IF(exists)THEN
            posic%cpdists(2,na) = dnrm2(3,xa-ip,1) ! MKL library routine
                                                   ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
            ! determine local coords
            CALL nbstLCoords(x1,x2,x3,ip,posic%cpxita(4:6,na),projectsL)
            
            IF(projectsL)THEN
              EXIT ! do-loop over distancies
            ENDIF
          ENDIF
        ENDDO
        
        IF(projectsU .AND. projectsL)THEN
          posic%cpReguDep(na) = .TRUE.
        ELSE ! aero node out of near structural mesh range - take nearest struc. node as master
          posic%cpReguDep(na) = .FALSE.
          
          ! determine irregular dependency
          CALL nbstIrrDep2(xa,distsUp,distsLo, &
                           eel, &
                           posic%cpmnods(1,na),posic%cpmnods(2,na),posic%cpmnods(3,na), &
                           posic%cpxita(1:4,na) )
          
          ! nullify useless variables
          posic%cpmnods(4:6,na) = 0
          posic%cpxita(5:6,na) = 0
          
          ! write warning
          msg = 'No proj. elem. for control point '//TRIM(inttoch(na,6))//' in '//TRIM(posic%m_aer)//' grid - stru. elem '//TRIM(inttoch(eel,6))//' used as master'
          CALL introut_mng('repProbl')
        ENDIF
    
      ENDDO
    ENDDO ! end with panel control points in root-transition
    
    
    ! 3. FINAL COMPUTATIONS
    !   3.1 NORMALIZE DISTANCIES FROM ELEMENTS
    ! aero grid nodes
    ALLOCATE( total(posic%nan) )
    total = posic%dists(1,:) + posic%dists(2,:)
    posic%dists(1,:) = posic%dists(1,:) / total
    posic%dists(2,:) = posic%dists(2,:) / total
    DEALLOCATE( total )
    ! control points
    ALLOCATE( total(posic%np) )
    total = posic%cpdists(1,:) + posic%cpdists(2,:)
    posic%cpdists(1,:) = posic%cpdists(1,:) / total
    posic%cpdists(2,:) = posic%cpdists(2,:) / total
    DEALLOCATE( total )
    
    !   3.2 DEALLOCATE INTERNAL ARRAYS
    DEALLOCATE(UpMidPo,LoMidPo,distsUp,distsLo)
    
    !   3.3 DETERMINE NUMBER OF CONTROL POINTS WITH IRREGULAR DEPENDENCY AND DATA ORDER IN NBSTH ARRAY
    DO i=1,posic%np
      IF(.NOT.posic%cpReguDep(i))THEN
        nIrrDCP = nIrrDCP+1
        posic%cpmnods(4,i) = nIrrDCP
      ENDIF
    ENDDO
    
    !CALL nbstMstrsTecplot( posic ) ! for debugging purposes only - file has to be deleted manually in order to be rewritten
    
  CONTAINS !...............................................
    
      SUBROUTINE nbstConnects(elmSetName, lnods, indxs)
        
        ! Determine array with connectivities (internal node labels) of a set
        ! of elements that belong to element set named elmSetName
        
        USE param_db,     ONLY: mnam
        USE ele13_db,     ONLY: srch_ele13,         & ! routines
                                ele13_set, ele13,   & ! types
                                head                  ! variables
      
        IMPLICIT NONE
        
        ! dummy args
        CHARACTER(len=mnam)    ,              INTENT(in)    ::  elmSetName  ! set name
        INTEGER, DIMENSION(:,:),              INTENT(out)   ::  lnods       ! elements connectivities
        INTEGER, DIMENSION(:)  , ALLOCATABLE, INTENT(out)   ::  indxs       ! index corresponding to element label in LNODS
        
        ! internal vars
        LOGICAL                                             ::  found
        INTEGER                                             ::  i
        INTEGER,            ALLOCATABLE,    DIMENSION(:)    ::  eLbls
        TYPE (ele13),       POINTER                         ::  ele
        TYPE (ele13_set),   POINTER                         ::  set, anterSet
        
        
        ! find element set
        CALL srch_ele13(head, anterSet, set, elmSetName, found) ! outputs SET pointing to set named elmSetName
                                                                ! in the list of NBST elements (whose first link is HEAD)
        IF(.NOT.found)THEN
          PRINT *, 'Element set named ', TRIM(elmSetName), ' not found - nbstConnects'
          PRINT *, 'program  S T O P P E D'
          STOP
        ENDIF
        
        ! extract and save connectivities of each element
        IF (ASSOCIATED (set%head)) THEN        ! Check if a list is empty
        
          ele => set%head
          ALLOCATE( eLbls(set%nelem) )
          i = 1
          DO ! Loop over elements
          
            lnods(:,i) = ele%lnods(1:3) ! internal labels of element nodes
            eLbls(i) = ele%numel ! external labels of elements
            
            i=i+1
            
            IF (ASSOCIATED (ele%next)) THEN    ! Check if it's not the end of the list
              ele => ele%next
            ELSE
              exit
            ENDIF
            
          ENDDO
          
          ALLOCATE( posic%elmIndxs(MINVAL(eLbls):MAXVAL(eLbls)) )
          DO i=1,set%nelem
            indxs(eLbls(i)) = i
          ENDDO
            
        ELSE
        
          PRINT *, 'Head element in set named ', TRIM(elmSetName), ' not found - nbstConnects'
          PRINT *, 'program  S T O P P E D'
          STOP
          
        ENDIF
        
        
            
      ENDSUBROUTINE nbstConnects
  
      !----------------------------------------------------------------------
    
      SUBROUTINE nbstMidpoints(elms, nelem, lnods, midpoints)
        
        ! Determines array with MIDPOINTS coords. of the NELEM elements
        ! whose (external) labels and connectivities are listed in ELMS
        ! and LNODS, respectively
              
        IMPLICIT NONE
        
        ! dummy args
        INTEGER, DIMENSION(:)       ::  elms        ! (external) labels of elements whose midpoints are to be determined
        INTEGER                     ::  nelem       ! number of elements and midpoints
        INTEGER, DIMENSION(:,:)     ::  lnods       ! elements connectivities
        REAL(8), DIMENSION(:,:)     ::  midpoints   ! array of midpoints coordinates
        
        ! internal vars
        INTEGER                     ::  iEle        ! index
        INTEGER                     ::  nodes(3)    ! labels
        
        
        DO i=1,nelem ! element index in BUSE or BLSE array
          !iEle = elms(i) ! element external label
          nodes(1:3) = lnods(:,i)!Ele)  ! internal labels of element nodes
          midpoints(1:3,i) = 1.0D0/3.0D0 * ( coord(:,nodes(1)) + coord(:,nodes(2)) + coord(:,nodes(3)) )
        ENDDO
        
      ENDSUBROUTINE nbstMidpoints
  
      !----------------------------------------------------------------------
      
      SUBROUTINE nbstIntersect(ns,posic,xa,n1,n2,n3,x1,x2,x3,ip,exists)
        
        ! find element, its nodes' coords and normal vector
        ! and
        ! determine intersection point
        
        IMPLICIT NONE
        
        ! dummy args
        INTEGER,                INTENT(in)  ::  ns
        TYPE(pair), POINTER,    INTENT(in)  ::  posic
        REAL(8),                INTENT(in)  ::  xa(3)
        INTEGER,                INTENT(out) ::  n1, n2, n3
        REAL(8),                INTENT(out) ::  x1(3), x2(3), x3(3), ip(3)
        LOGICAL,				INTENT(out) ::  exists
        
        ! local vars
        INTEGER :: ei
        REAL(8) :: v1(3), v2(3), n(3), r(3), s
        LOGICAL :: parallel, coincide
        
        
        ! find element nodes, its coords and normal vector
        ei = posic%elmIndxs(ns)  ! element index
        !--
        n1 = posic%lnods(1,ei) ! 1st node
        n2 = posic%lnods(2,ei) ! 2nd node
        n3 = posic%lnods(3,ei) ! 3rd node
        !--
        x1 = coord(1:3,n1) ! 1st node coordinates
        x2 = coord(1:3,n2) ! 2nd node coordinates
        x3 = coord(1:3,n3) ! 3rd node coordinates
        !--
        v1 = x2-x1
        v2 = x3-x1
        CALL cross_product(n, v1, v2) ! cross_product(out, in1, in2)
        ! check parallelism
        parallel = DABS(DOT_PRODUCT(n,thick)) .LT. 1D-5 ! .TRUE. if n and thick are normal
        IF(parallel)THEN
          ! check coincidence
          coincide = DABS(DOT_PRODUCT(n,xa-x1)) .LT. 1D-5
          IF(coincide)THEN
            exists = .TRUE.
            ip = xa
          ELSE
            exists = .FALSE.
          ENDIF
        ELSE
          exists = .TRUE.
          ! determine intersection point
          r = xa-x1
          s = - DOT_PRODUCT(n,r)/DOT_PRODUCT(n,thick)
          ip = xa+s*thick ! intersection point
        ENDIF
        
      ENDSUBROUTINE nbstIntersect
  
      !----------------------------------------------------------------------
      
      SUBROUTINE nbstLCoords(x1,x2,x3,xyzp,lcoor,projects)
        
        ! determines Local COORDinateS (LCOOR) of point with coords. XYZP
        ! in NBST triangle element defined by nodes whose coordinates are
        ! x1, x2 and x3
        ! (node order called i,j,k)
        ! (point P must lay in the plane of the triangle)
        
        IMPLICIT NONE
        
        ! dummy args
        LOGICAL,    INTENT(out)     ::  projects
        REAL(8),    INTENT(in)      ::  x1(3),x2(3),x3(3),xyzp(3)
        REAL(8),    INTENT(out)     ::  lcoor(3)
        
        ! internal vars
        REAL(8) :: vj(3), vk(3), vn(3), vt(3), auxi(3), area2
        
        ! side vectors
        vj = x2 - x1 ! side i-->j
        vk = x3 - x1 ! side i-->k
        ! normal vector  Vn = Vj x Vk
        CALL cross_product(vn, vj, vk) ! cross_product(out, in1, in2)
        ! compute twice the area of the triangle and normalizes Vn
        area2 = dnrm2(3,vn,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
        vn = vn/area2 ! unit normal
        ! relative position of P wrt LNODS(1)
        vt = xyzp - x1
        !  3rd local coordinate
        CALL cross_product(auxi, vj, vt) ! cross_product(out, in1, in2) - auxi = Vj x Vt
        lcoor(3) = DOT_PRODUCT(auxi,vn)/area2 ! 3rd local coordinate (associated with node k)
        !  2nd local coordinate
        CALL cross_product(auxi, vt, vk) ! cross_product(out, in1, in2) - auxi = Vt x Vk
        lcoor(2) = DOT_PRODUCT(vn,auxi)/area2 ! 2nd local coordinate (associated with node j)
        !  1st local coordinate
        lcoor(1) = 1d0 - lcoor(2) - lcoor(3)  ! 1st local coordinate (associated with node i)
        
        ! determine projection
        projects = lcoor(1) >= 0d0 .AND. lcoor(2) >= 0d0 .AND. lcoor(3) >= 0d0
        
      ENDSUBROUTINE nbstLCoords
  
      !----------------------------------------------------------------------
      
      SUBROUTINE nbstIrrDep2(xyzp,distsUp,distsLo,eel,n1,n2,n3,lcoor)
        
        ! irregular dependency
        ! aero node out of structural mesh, depending on an element
        ! using triangular local coords. and versor normal to element plane
        
        IMPLICIT NONE
        
        ! dummy args
        REAL(8),                    INTENT(in)  ::  xyzp(3)
        REAL(8),    DIMENSION(:),   INTENT(in)  ::  distsUp, distsLo
        INTEGER,                    INTENT(out) ::  eel, n1, n2, n3
        REAL(8),    DIMENSION(4),   INTENT(out) ::  lcoor
        
        ! internal vars
        LOGICAL                     ::  projects
        INTEGER                     ::  ind, eci
        REAL(8)                     ::  dUp, dLo, area2
        REAL(8),    DIMENSION(3)    ::  x1, x2, x3, l1, l2, t, r, pp
        
        ! FIND NEAREST ELEMENT (MIDPOINT) AND ITS NODES' COORDS
        dUp = MINVAL( distsUp, 1 ) ! distance to nearest element midpoint
        dLo = MINVAL( distsLo, 1 ) ! distance to nearest element midpoint
        IF(dUp .LE. dLo)THEN
          ind = MINLOC( distsUp, 1 ) ! index of nearest element midpoint in posic%bed%buse array
          eel = posic%bed%buse(ind) ! Element External Label
        ELSE
          ind = MINLOC( distsLo, 1 ) ! index of nearest element midpoint in posic%bed%blse array
          eel = posic%bed%blse(ind) ! Element External Label
        ENDIF
        eci = posic%elmIndxs(eel) ! Elemen Connectivities array Index
        !--
        n1 = posic%lnods(1,eci) ! 1st node
        n2 = posic%lnods(2,eci) ! 2nd node
        n3 = posic%lnods(3,eci) ! 3rd node
        !--
        x1 = coord(:,n1) ! 1st node coordinates
        x2 = coord(:,n2) ! 2nd node coordinates
        x3 = coord(:,n3) ! 3rd node coordinates
        
        ! NORMAL VERSOR
        l1 = x3-x2
        l2 = x1-x3
        CALL cross_product(t, l1, l2) ! cross_product(out, in1, in2)
        ! compute twice the area of the triangle and normalize t
        area2 = dnrm2(3,t,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
        t = t/area2 ! unit normal
        
        ! DISTANCE BETWEEN POINT AND PLANE
        r = xyzp-x1
        lcoor(4) = DOT_PRODUCT(t,r) ! signed distance
        
        ! PROJECTION OF POINT ON PLANE
        pp = xyzp-lcoor(4)*t
        
        ! LOCAL COORDS. OF PROJECTION
        CALL nbstLCoords(x1,x2,x3,pp,lcoor(1:3),projects)
        
      ENDSUBROUTINE nbstIrrDep2
  
      ! ------------------------------- debug -------------------------------
      ! =====================================================================
      
      SUBROUTINE nbstMstrsTecplot( posic )
        
        ! prints TecPlot file with connections between aero nodes and its
        ! structural nodes masters
        
        ! prints a commented section (at the beginning) with the data in a table
        
        ! Mauro S. Maza - 23/10/2015
        
        IMPLICIT NONE
        
        ! dummy args
        TYPE(pair), POINTER :: posic
        
        ! internal vars
        LOGICAL                             ::  fExist
        INTEGER                             ::  npoin, nelem, nIrrD, i, j, n1, lowLim, upeLim
        INTEGER, ALLOCATABLE, DIMENSION(:)  ::  nodsIndxs
        
        
        INQUIRE(FILE='nbstMstrsTecplot.tec', exist=fExist)
        IF(fExist)THEN
          ! nothing
        ELSE
          OPEN(UNIT = 242, FILE='nbstMstrsTecplot.tec', STATUS='new', ACTION='write')
          
          ! DATA TABLE
          ! header
          WRITE(242,400)
          WRITE(242,500,advance='no')  '# aero    regular  masters                                         local coordinates                                                       distancies       '
          WRITE(242,500,advance='no')  '# node    depen_   upper surface           lower surface           upper surface                       lower surface                       upper       lower'
          WRITE(242,500,advance='no')  '# label   dency    node#1  node#2  node#3  node#1  node#2  node#3  node#1      node#2      node#3      node#1      node#2      node#3      surf.       surf.'
          WRITE(242,500,advance='no')  '# AERO GRID NODES ------------------------------------------------------------------------------------------------------------------------------------------'
          DO i=1,posic%nan
            WRITE(242,600,advance='no')  '#', i, posic%reguDep(i), posic%mnods(:,i), posic%xita(:,i), posic%dists(:,i)
          ENDDO
          WRITE(242,500,advance='no')  '# CONTROL POINTS -------------------------------------------------------------------------------------------------------------------------------------------'
          DO i=1,posic%np
            WRITE(242,600,advance='no')  '#', i, posic%cpReguDep(i), posic%cpmnods(:,i), posic%cpxita(:,i), posic%cpdists(:,i)
          ENDDO
          WRITE(242,400)
          
          
          ! ACTUAL TECPLOT DATA FILE -------------------
          !   GRID NODES *******************************
          ! determine some header data
          npoin = posic%nan+posic%numbnodes
          nelem = 0
          nIrrD = 0 ! number of aero nodes with irregular dependency
          lowLim = 1
          upeLim = posic%nan
          DO i=lowLim,upeLim
            IF(posic%reguDep(i))THEN
              nelem = nelem+6
            ELSE
              nelem = nelem+3
              nIrrD = nIrrD+1
            ENDIF
          ENDDO
          
          ! SLAVE-MASTER LINKS
          ! header
          WRITE( 242 , * ) 'ZONE T = "Slave-Master Links (grid)" , F = FEPOINT, N = ', npoin, ' , E = ', nelem, ' , ET = LINESEG'
          ! aero nodes
          DO i=1,posic%nan
            WRITE( 242 , 100 ) posic%nodsaer(i)%xyz(1), posic%nodsaer(i)%xyz(2), posic%nodsaer(i)%xyz(3)
          ENDDO
          ! structural nodes and its order
          ALLOCATE( nodsIndxs(MINVAL(posic%nodes):MAXVAL(posic%nodes)) )
          nodsIndxs = 0
          DO i=1,posic%numbnodes
            n1 = posic%nodes(i)
            nodsIndxs(n1) = i ! inverted relation between node labes and indexes of that expressed in posic%nodes
            !xyz_e1 = coord(:,n1)
            WRITE( 242 , 100 ) coord(:,n1) !xyz_e1(1) , xyz_e1(2) , xyz_e1(3)
          ENDDO
          WRITE( 242 , 400 )
          ! connectivities
          DO i=lowLim,upeLim
            IF(posic%reguDep(i))THEN
              DO j=1,6
                n1 = posic%mnods(j,i)
                WRITE( 242 , 200 ) i, nodsIndxs(n1)+posic%nan
              ENDDO
            ELSE
              DO j=1,3
                n1 = posic%mnods(j,i)
                WRITE( 242 , 200 ) i, nodsIndxs(n1)+posic%nan
              ENDDO
            ENDIF
          ENDDO
          WRITE( 242 , 400 )
        
          ! SLAVES WITH IRREGULAR DEPENDENCY
          ! header
          IF(nIrrD.NE.0)THEN
            WRITE( 242 , * ) 'ZONE T = "IrrDep Nodes (grid)" , I = ', nIrrD
            DO i=lowLim,upeLim
              IF(.NOT.posic%reguDep(i))THEN
                WRITE( 242 , 100 ) posic%nodsaer(i)%xyz(1), posic%nodsaer(i)%xyz(2), posic%nodsaer(i)%xyz(3)
              ENDIF
            ENDDO
            WRITE( 242 , 400 )
          ENDIF
          
          
          !   CONTROL POINTS ***************************
          ! determine some header data
          npoin = posic%np+posic%numbnodes
          nelem = 0
          nIrrD = 0 ! number of aero nodes with irregular dependency
          lowLim = 1
          upeLim = posic%np
          DO i=lowLim,upeLim
            IF(posic%cpReguDep(i))THEN
              nelem = nelem+6
            ELSE
              nelem = nelem+3
              nIrrD = nIrrD+1
            ENDIF
          ENDDO
          
          ! SLAVE-MASTER LINKS
          ! header
          WRITE( 242 , * ) 'ZONE T = "Slave-Master Links (CPs)" , F = FEPOINT, N = ', npoin, ' , E = ', nelem, ' , ET = LINESEG'
          ! aero nodes
          DO i=1,posic%np
            WRITE( 242 , 100 ) posic%panel(i)%xyzcp(1), posic%panel(i)%xyzcp(2), posic%panel(i)%xyzcp(3)
          ENDDO
          ! structural nodes
          DO i=1,posic%numbnodes
            n1 = posic%nodes(i)
            WRITE( 242 , 100 ) coord(:,n1)
          ENDDO
          WRITE( 242 , 400 )
          ! connectivities
          DO i=lowLim,upeLim
            IF(posic%cpReguDep(i))THEN
              DO j=1,6
                n1 = posic%cpmnods(j,i)
                WRITE( 242 , 200 ) i, nodsIndxs(n1)+posic%np
              ENDDO
            ELSE
              DO j=1,3
                n1 = posic%cpmnods(j,i)
                WRITE( 242 , 200 ) i, nodsIndxs(n1)+posic%np
              ENDDO
            ENDIF
          ENDDO
          WRITE( 242 , 400 )
        
          ! SLAVES WITH IRREGULAR DEPENDENCY
          ! header
          IF(nIrrD.NE.0)THEN
            WRITE( 242 , * ) 'ZONE T = "IrrDep Nodes (CPs)" , I = ', nIrrD
            DO i=lowLim,upeLim
              IF(.NOT.posic%cpReguDep(i))THEN
                WRITE( 242 , 100 ) posic%panel(i)%xyzcp(1), posic%panel(i)%xyzcp(2), posic%panel(i)%xyzcp(3)
              ENDIF
            ENDDO
            WRITE( 242 , 400 )
          ENDIF
          
        ENDIF
        
        CLOSE( UNIT = 242 )
        
        
        100 FORMAT(3(1X,F9.3))  ! coords.
        200 FORMAT(2(1X,I5))    ! 2 node connects.
        300 FORMAT(3(1X,I6))    ! 3 node connects.
        400 FORMAT(/)           ! new line
        500 FORMAT(A/)          ! text and new line
        600 FORMAT(A,x,I6,2x,L1,8x,6(I6,2x),8(E10.3E2,2x),/) ! table
        
      ENDSUBROUTINE nbstMstrsTecplot
      
  ENDSUBROUTINE nbstMstrs
  
  !--------------------------------------------------------------------------
  
  ! ----------------------- general  ----------------------
  ! =======================================================
  
  SUBROUTINE distancies(dists, pointSet, targt)
    
    IMPLICIT NONE
    
    
    ! dummy args
    REAL,   INTENT(out)  :: dists(:)
    REAL,   INTENT(in)   :: pointSet(:,:)
    REAL,   INTENT(in)   :: targt(3)
    
    ! internal vars
    INTEGER                             ::  i
    REAL, ALLOCATABLE, DIMENSION(:,:)   ::  relPoss
    
    
    ALLOCATE( relPoss(3,SIZE(pointSet,DIM=2)) )
    DO i=1,3
      relPoss(i,:) = pointSet(i,:) - targt(i)
    ENDDO
    DO i= LBOUND(relPoss,2), UBOUND(relPoss,2)    ! loop over relative positions
      dists(i) = dnrm2(3, relPoss(:,i), 1)  ! MKL library routine
                                            ! ?nrm2(n,x,s): euclidean norm of x(n) with stride s
    ENDDO
    DEALLOCATE( relPoss )
    
  ENDSUBROUTINE distancies
  
ENDMODULE mastrnods_sr

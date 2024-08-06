
MODULE inter_sr
  
  USE npo_db,         ONLY:   coora, euler, velnp
  USE inter_db
  USE Constants
  USE introut_db,     ONLY:   msg
  USE introut_sr,     ONLY:   introut_mng
  USE param_db,       ONLY:   mich
  
  IMPLICIT NONE
  
  INCLUDE 'mkl_blas.fi' ! INTERFACE blocks for BLAS/MKL routines
  
CONTAINS ! ==================================================================
  
  ! ----------- Subroutines for interaction ------------
  ! ====================================================
  
  SUBROUTINE inter_ini
    
    ! INTERaction_INItial computations
    
    ! Mauro S. Maza - 15/09/2015
    
    ! Initial calculations for interaction
    ! Determine:    global coordinates for aero. nodes and CPs
    !               projection element for each aerodynamical node and CP
    !               projection position in the element (local coordinates)
    !               initial distance in master nodes local reference frames
    
    USE StructuresA,    ONLY:   phi, beta, theta1, theta2, theta3, &
                                q1, q2, q3, q4, q5, &
                                nn, nnh, nnn, nnt, &
                                np, nph, npn, npt, &
                                dens
    USE mastrnods_sr,   ONLY:   mngMstrs ! this implies an interdependence between inter_sr and mastrnods_sr modules, which generates a compilation problem when solution rebuilded - it can be avoided commenting this line for a first compilation
    USE wind_sr,        ONLY:   DtCalc
    USE ele13_db,       ONLY:   nods13, srch_ele13, srch_ele13e,    & ! routines
                                ele13_set, ele13,                   & ! types
                                head                                  ! variables
    USE ctrl_db,        ONLY:   npoin
    USE npo_db,         ONLY:   coord
    
    IMPLICIT NONE
    
    LOGICAL                             ::  found
    INTEGER                             ::  labels(npoin)   ! in order to find INTERNAL labels with NODS13
    INTEGER,            POINTER         ::  nodes(:)        ! subroutine NODS13 dummy argument is pointer
    TYPE(pair),         POINTER         ::  posic
    TYPE(ele13_set),    POINTER         ::  anterset, set
    
    
    ! Step variables for data flow control - both first used in dynamic.f90
    nest = 0        ! current struct. step inside an aero. step - between tow aero. calculations
    numnest = 0     ! number of struct. steps for each aero. step
    
    ! Reference quantities
    CALL Geo_references        ! length and area
    CF2F = 0.5D+0 * dens       ! Force coefficient
    
    ! Local coordinates for blades 2 and 3
    CALL Copy_blades
    
    ! Initial global coordinates for aero. nodes
    phi    = phi    * deg2rad
    beta   = beta   * deg2rad
    theta1 = theta1 * deg2rad
    theta2 = theta2 * deg2rad
    theta3 = theta3 * deg2rad
    q1     = q1     * deg2rad
    q2     = q2     * deg2rad
    q3     = q3     * deg2rad
    q4     = q4     * deg2rad
    q5     = q5     * deg2rad
    CALL Blade1_Kinematics
    CALL Blade2_Kinematics
    CALL Blade3_Kinematics
    CALL Hub_Kinematics
    CALL Nacelle_Kinematics
    CALL Tower_Kinematics
    
    ! Global coordinates for Control Points
    CALL CP_coords
    
    ! Pairing
    IF (ASSOCIATED (headpair)) THEN        ! Check if a list is empty
      nIrrDCP = 0
      posic => headpair
      DO ! Loop over pairs
        
        SELECT CASE (TRIM(posic%m_aer))
        CASE ('blade1')
          posic%eType = TRIM(bEType)        ! element type for structural model
          posic%nan = nn                    ! Number of Aero. Nodes
          posic%nodsaer => blade1nodes      ! Array with nodal info of posic%m_aer
          posic%np = np                     ! Number of Control Points
          posic%panel   => blade1panels     ! Array with panels' info of posic%m_aer
          posic%section => blade1sections   ! Array with sections info of blade1 lifting surface
          
        CASE ('blade2')
          posic%eType = TRIM(bEType)        ! element type for structural model
          posic%nan = nn                    ! Number of Aero. Nodes
          posic%nodsaer => blade2nodes      ! Array with nodal info of posic%m_aer
          posic%np = np                     ! Number of Control Points
          posic%panel   => blade2panels     ! Array with panels' info of posic%m_aer
          posic%section => blade2sections   ! Array with sections info of blade2 lifting surface
          
        CASE ('blade3')
          posic%eType = TRIM(bEType)        ! element type for structural model
          posic%nan = nn                    ! Number of Aero. Nodes
          posic%nodsaer => blade3nodes      ! Array with nodal info of posic%m_aer
          posic%np = np                     ! Number of Control Points
          posic%panel   => blade3panels     ! Array with panels' info of posic%m_aer
          posic%section => blade3sections   ! Array with sections info of blade3 lifting surface
          
        CASE ('hub')
          posic%eType = 'RIGID'        ! element type for structural model
          posic%nan = nnh              ! Number of Aero. Nodes
          posic%nodsaer => hubnodes    ! Array with nodal info of posic%m_aer
          posic%np = nph               ! Number of Control Points
          posic%panel   => hubpanels   ! Array with panels' info of posic%m_aer
        CASE ('nacelle')
          posic%eType = 'BEAM '            ! element type for structural model - nacelle aero. grid is linked to union estr. mesh
          posic%nan = nnn                  ! Number of Aero. Nodes
          posic%nodsaer => nacellenodes    ! Array with nodal info of posic%m_aer
          posic%np = npn                   ! Number of Control Points
          posic%panel   => nacellepanels   ! Array with panels' info of posic%m_aer
        CASE ('tower')
          posic%eType = 'BEAM '          ! element type for structural model
          posic%nan = nnt                ! Number of Aero. Nodes
          posic%nodsaer => towernodes    ! Array with nodal info of posic%m_aer
          posic%np = npt                 ! Number of Control Points
          posic%panel   => towerpanels   ! Array with panels' info of posic%m_aer
        ENDSELECT
        
        
        IF(TRIM(posic%eType)=='NBST')THEN ! shell model for blades
          CALL findBed(posic)             ! Blade Elements Distribution
          CALL srch_ele13(head, anterSet, set, posic%m_est, found) ! outputs SET pointing to set named elmSetName=posic%m_est
                                                                   ! in the list of NBST elements (whose first link is HEAD)
          posic%nelem = set%nelem
          labels = [1:npoin] ! in order to find INTERNAL labels with NODS13
          CALL nods13(npoin, set%nelem, posic%numbnodes, nodes, set%head, labels) ! NODES is vector with LABELS of nodes defining the elements of ELSET
          ALLOCATE( posic%nodes(posic%numbnodes) )
          posic%nodes(:) = nodes(:) ! subroutine NODS13 dummy argument is pointer
        ENDIF
        
        
        CALL mngMstrs( posic )
        
        
        IF (ASSOCIATED (posic%next)) THEN    ! Check if it's not the end of the list
          posic => posic%next
        ELSE
          exit
        ENDIF
        
      ENDDO
    ENDIF
    
    IF(nIrrDCP.NE.0)THEN
      ! allocate array
      ALLOCATE(nbstH(27,nIrrDCP))
    ENDIF
    
    ! More reference quantities
    CALL DtCalc(.FALSE.)    ! aero. time step - needs hub_master to be determined
    
    
  ENDSUBROUTINE inter_ini
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_kin
    
    ! INTERaction_KINematic
  
    ! Mauro S. Maza - 11/10/2011
    !                 04/01/2016
    
    ! Determines position and velocity of aero. nodes
    ! and control points from structural nodes data
    
    IMPLICIT NONE
    
    ! internal vars
    TYPE(pair),         POINTER     ::  posic
   
    INTEGER (kind=4)                ::  n1, n2, nn1(3), nn2(3)
    INTEGER                         ::  i, j
    
    REAL(kind = 8)                  ::  d1(3), d2(3), d1l(3), d2l(3),       &
                                        o1(3), o2(3), o1l(3), o2l(3),       &
                                        v1(3), v2(3),                       &
                                        lin_vel1(3), lin_vel2(3), xyz1(3),  &
                                        lambda1(3,3), lambda2(3,3),         &
                                        xx1(3,3), xx2(3,3),                 &
                                        vv1(3,3), vv2(3,3)
    REAL(8)                         ::  l1(3), l2(3), l3(3), t(3), area2,    &
                                        tRow(1,3), h(3,9), no1(3), no2(3), no3(3)
    REAL                            ::  xita1
    
    
    IF( ASSOCIATED(headpair) )THEN      ! Check if a list is empty
      posic => headpair
      
      DO ! Loop over pairs
        SELECT CASE (TRIM(posic%eType))
        CASE('RIGID') ! ---------------------
          ! General data
          n1 = posic%mnods(1,1)       ! master node
          lambda1 = RESHAPE(euler(1:9,n1),(/3,3/))	! Rotation matrix for master node local system
          xyz1 = coora(:,n1)   ! master node position vector
          v1 = velnp(1:3,n1)   ! master node linear velocity vector
          o1l = velnp(4:6,n1)  ! master node angular velocity vector (rotation matrix axial vector) in local coordinates of master node
          ! Aero. nodes
          DO i=1,posic%nan
            ! Position
            d1l = posic%dists(1:3,i)    ! Initial distance vector in local coords. of master node n1
            d1 = MATMUL(lambda1,d1l)    ! Distance vector from master node to aero. node n1-->an
            posic%nodsaer(i)%xyz(1:3) = xyz1 + d1   ! aero node position
		  ENDDO
          ! Control Points
          DO i=1,posic%np
            ! Position
            d1l = posic%cpdists(1:3,i)  ! Initial distance vector in local coords. of master node n1
            d1 = MATMUL(lambda1,d1l)    ! Distance vector from master node to aero. node n1-->an
            posic%panel(i)%xyzcp = xyz1 + d1    ! aero node position
            ! Velocity
            o1 = MATMUL(lambda1,o1l)  ! master node angular velocity vector (rotation matrix axial vector) in global coordinates
            CALL cross_product( lin_vel1 , o1 , d1 )   ! aero. node linear velocity due to master node angular velocity: o1 x d1l = lin_vel
            posic%panel(i)%velcp = v1 + lin_vel1   ! aero. node linear velocity
		  ENDDO
          
        CASE('BEAM') ! ---------------------
          ! Aero. nodes
          DO i=1,posic%nan
            ! General data
            n1 = posic%mnods(1,i)       ! master node n1
            n2 = posic%mnods(2,i)       ! master node n2
            d1l = posic%dists(1:3,i)    ! Initial distance vector in local coords. of master node n1
            d2l = posic%dists(4:6,i)    ! Initial distance vector in local coords. of master node n2
            xita1 = posic%xita(1,i)       ! coordinate for projection (relative to n1)            
            lambda1 = RESHAPE(euler(1:9,n1),(/3,3/))	! Rotation matrix for n1 local system
            lambda2 = RESHAPE(euler(1:9,n2),(/3,3/))	! Rotation matrix for n2 local system
            ! Position
            d1 = MATMUL(lambda1,d1l)    ! Distance vector from master node to aero. node n1-->an
            d2 = MATMUL(lambda2,d2l)    ! Distance vector from master node to aero. node n2-->an
            posic%nodsaer(i)%xyz = (1-xita1) * ( coora(:,n1) + d1 ) +  & ! aero node position
                                      xita1  * ( coora(:,n2) + d2 )
          ENDDO
          ! Control Points
          DO i=1,posic%np
            ! General data
            n1 = posic%cpmnods(1,i)       ! master node n1
            n2 = posic%cpmnods(2,i)       ! master node n2
            d1l = posic%cpdists(1:3,i)    ! Initial distance vector in local coords. of master node n1
            d2l = posic%cpdists(4:6,i)    ! Initial distance vector in local coords. of master node n2
            xita1 = posic%cpxita(1,i)       ! coordinate for projection (relative to n1)
            lambda1 = RESHAPE(euler(1:9,n1),(/3,3/))	! Rotation matrix for n1 local system
            lambda2 = RESHAPE(euler(1:9,n2),(/3,3/))	! Rotation matrix for n2 local system
            ! Position
            d1 = MATMUL(lambda1,d1l)    ! Distance vector from master node to aero. node n1-->an
            d2 = MATMUL(lambda2,d2l)    ! Distance vector from master node to aero. node n2-->an
            posic%panel(i)%xyzcp = (1-xita1) * ( coora(:,n1) + d1 ) +  & ! aero node position
                                      xita1  * ( coora(:,n2) + d2 )
            ! Velocity
            v1 = velnp(1:3,n1)  ! master node n1 linear velocity vector
            v2 = velnp(1:3,n2)  ! master node n2 linear velocity vector
            o1l = velnp(4:6,n1)  ! master node n1 angular velocity vector (rotation matrix axial vector) in local coordinates of master node n1
            o2l = velnp(4:6,n2)  ! master node n2 angular velocity vector (rotation matrix axial vector) in local coordinates of master node n2
            o1 = MATMUL(lambda1,o1l)    ! master node n1 angular velocity vector (rotation matrix axial vector) in global coordinates
            o2 = MATMUL(lambda2,o2l)    ! master node n2 angular velocity vector (rotation matrix axial vector) in global coordinates
            CALL cross_product( lin_vel1 , o1 , d1 ) ! aero. node linear velocity due to master node angular velocity: o1 x d1l = lin_vel1
            CALL cross_product( lin_vel2 , o2 , d2 ) ! aero. node linear velocity due to master node angular velocity: o2 x d2l = lin_vel2
            posic%panel(i)%velcp = (1-xita1) * ( v1 + lin_vel1 ) +  & ! aero node linear velocity
                                      xita1  * ( v2 + lin_vel2 )
		  ENDDO
          
        CASE('NBST') ! ---------------------
          ! Aero. nodes
          DO i=1,posic%nan
            ! data
            nn1(1:3) = posic%mnods(1:3,i)
            xx1(1:3,1:3) = coora(:,nn1)
            
            IF(posic%reguDep(i))THEN
              ! more data
              nn2(1:3) = posic%mnods(4:6,i)
              xx2(1:3,1:3) = coora(:,nn2)
              
              ! position
              ! note that the order of dists(1:2,i) is inverted because they act as weights in the sum
              posic%nodsaer(i)%xyz = posic%dists(2,i)*MATMUL(xx1,posic%xita(1:3,i)) + posic%dists(1,i)*MATMUL(xx2,posic%xita(4:6,i))
            ELSE
              ! NORMAL VERSOR
              l1 = xx1(:,3) - xx1(:,2)
              l2 = xx1(:,1) - xx1(:,3)
              CALL cross_product(t, l1, l2) ! cross_product(out, in1, in2)
              ! compute twice the area of the triangle and normalize t
              area2 = dnrm2(3,t,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
              t = t/area2 ! unit normal
              
              ! POSITION
              !   SUM Li*xi + d*t
              !   Li= posic%xita(1:3,i)
              !   xi=xx1
              !   Li*xi= MATMUL(xx1,posic%xita(1:3,i))
              !   d=posic%xita(4,i)
              posic%nodsaer(i)%xyz = MATMUL(xx1,posic%xita(1:3,i)) + posic%xita(4,i)*t
            ENDIF
          ENDDO
          
          ! Control Points
          DO i=1,posic%np
            ! data
            nn1(1:3) = posic%cpmnods(1:3,i)
            xx1(1:3,1:3) = coora(:,nn1)
            vv1(1:3,1:3) = velnp(1:3,nn1)
            
            IF(posic%cpReguDep(i))THEN
              ! more data
              nn2(1:3) = posic%cpmnods(4:6,i)
              xx2(1:3,1:3) = coora(:,nn2)
              vv2(1:3,1:3) = velnp(1:3,nn2)
              
              ! position
              ! note that the order of cpdists(1:2,i) is inverted because they act as weights in the sum
              posic%panel(i)%xyzcp = posic%cpdists(2,i)*MATMUL(xx1,posic%cpxita(1:3,i)) + posic%cpdists(1,i)*MATMUL(xx2,posic%cpxita(4:6,i))
              
              ! velocity
              ! note that the order of cpdists(1:2,i) is inverted because they act as weights in the sum
              posic%panel(i)%velcp = posic%cpdists(2,i)*MATMUL(vv1,posic%cpxita(1:3,i)) + posic%cpdists(1,i)*MATMUL(vv2,posic%cpxita(4:6,i))
            ELSE
              ! NORMAL VERSOR
              l1 = xx1(:,3) - xx1(:,2)
              l2 = xx1(:,1) - xx1(:,3)
              CALL cross_product(t, l1, l2) ! cross_product(out, in1, in2)
              ! compute twice the area of the triangle and normalize t
              area2 = dnrm2(3,t,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
              t = t/area2 ! unit normal
              
              ! POSITION
              !   SUM Li*xi + d*t
              !   Li= posic%xita(1:3,i)
              !   xi=xx1
              !   Li*xi= MATMUL(xx1,posic%xita(1:3,i))
              !   d=posic%xita(4,i)
              posic%panel(i)%xyzcp = MATMUL(xx1,posic%cpxita(1:3,i)) + posic%cpxita(4,i)*t
              
              ! VELOCITY
              !   SUM Hi*vi
              !   Hi=Id3x3*Li+d/2A*kron(ni,t)
              !   ni=li x t
              l3 = xx1(:,2) - xx1(:,1)
              CALL cross_product(no1, l1, t) ! cross_product(out, in1, in2)
              CALL cross_product(no2, l2, t) ! cross_product(out, in1, in2)
              CALL cross_product(no3, l3, t) ! cross_product(out, in1, in2)
              tRow = RESHAPE(t,[1,3])
              h(1:3,1:3) = MATMUL(RESHAPE(no1,[3,1]),tRow)
              h(1:3,4:6) = MATMUL(RESHAPE(no2,[3,1]),tRow)
              h(1:3,7:9) = MATMUL(RESHAPE(no3,[3,1]),tRow)
              h = h * posic%cpxita(4,i) / area2
              FORALL(j=1:3) h(j,j  ) = h(j,j  )+posic%cpxita(1,i)
              FORALL(j=1:3) h(j,j+3) = h(j,j+3)+posic%cpxita(2,i)
              FORALL(j=1:3) h(j,j+6) = h(j,j+6)+posic%cpxita(3,i)
              
              posic%panel(i)%velcp = MATMUL(h,RESHAPE(vv1,[9]))
              
              nbstH(1:27, posic%cpmnods(4,i) ) = RESHAPE(h,[27])
            ENDIF
        ENDDO
        
        ENDSELECT
        
        IF (ASSOCIATED (posic%next)) THEN   ! Check if it's not the end of the list
		  posic => posic%next
		ELSE
		  exit
		ENDIF
        
      ENDDO
    ENDIF
    
    !CALL CPsVelsText ! for debugging purposes only - file is never deleted by the program; data is simply added to the end - you may delet it manually
    !CALL nbstOmTecplot ! for debugging purposes only - file is never deleted by the program; data is simply added to the end - you may delet it manually
  
  CONTAINS
      
      ! ------------------------------- debug -------------------------------
      ! =====================================================================
      
      SUBROUTINE CPsVelsText
        
        ! prints txt file with translational velocities of control points
        
        ! Mauro S. Maza - 26/07/2016
        
        USE ctrl_db,        ONLY:   ttime
        
        IMPLICIT NONE
        
        ! internal vars
        LOGICAL                             ::  fExist
        INTEGER                             ::  i
        
        
        INQUIRE(FILE='CPsVelsText.txt', exist=fExist)
        IF(fExist)THEN
          OPEN(UNIT = 248, FILE='CPsVelsText.txt', STATUS='old', POSITION='append', ACTION='write')
        ELSE
          OPEN(UNIT = 248, FILE='CPsVelsText.txt', STATUS='new', ACTION='write')
          ! Heading
          WRITE( 248 , 400 )
          WRITE( 248 , * ) 'translational velocities of control points'
          WRITE( 248 , * ) '******************************************'
          WRITE( 248 , * ) 'ONLY BLADE2 NOW'
          WRITE( 249 , 400 )
          WRITE( 249 , * ) 'for debugging purposes only'
          WRITE( 249 , * ) 'file is never deleted by the program; data is simply added to the end'
          WRITE( 249 , * ) 'you may delet it manually'
          WRITE( 249 , 400 )
        ENDIF
        
        ! step title
        WRITE( 248 , * ) 'aeroStep            time'
        WRITE( 248 , * ) n, ttime
        IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	      posic => headpair
          DO ! Loop over pairs
            
            IF( TRIM(posic%m_aer)=='blade2' )THEN
              DO i=1,posic%np
                WRITE(248,600) posic%panel(i)%velcp(1:3)
              ENDDO
              WRITE( 244 , 400 )
            ENDIF
            
            IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
              posic => posic%next
            ELSE
              exit
            ENDIF
	    	
          ENDDO
        ENDIF
        
        CLOSE(UNIT = 248)
        
        400 FORMAT(/) ! new line
        600 FORMAT(3(1X,E12.5E2)) ! 3 times vel. comps.
        
      ENDSUBROUTINE CPsVelsText
      
      ! ---------------------------------------------------------------------
      
      SUBROUTINE nbstOmTecplot
        
        ! prints Tecplot file with translational velocities of structural
        ! nodes and control points
        
        ! Mauro S. Maza - 12/11/2015
        
        USE ctrl_db,        ONLY:   ttime
        
        IMPLICIT NONE
        
        ! internal vars
        LOGICAL                             ::  fExist
        INTEGER                             ::  strandID, i, nRegD, nIrrD
        
        
        INQUIRE(FILE='nbstOmTecplot.tec', exist=fExist)
        IF(fExist)THEN
          OPEN(UNIT = 244, FILE='nbstOmTecplot.tec', STATUS='old', POSITION='append', ACTION='write')
        ELSE
          OPEN(UNIT = 244, FILE='nbstOmTecplot.tec', STATUS='new', ACTION='write')
          ! Heading
          WRITE( 244 , * ) 'TITLE = "Structural nodes and control points translational velocities"'
          WRITE( 244 , * ) 'VARIABLES = "X" "Y" "Z" "v_X" "v_Y" "v_Z"'
          WRITE( 244 , 400 )
        ENDIF
        
        strandID = 0
        
        IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	      posic => headpair
          DO ! Loop over pairs
            IF( TRIM(posic%eType)=='NBST' )THEN
              
              ! STRUCTURAL NODES
              strandID = strandID + 1
              WRITE( 244 , * ) 'ZONE T = "StruNodes - ', TRIM(posic%m_est),' n=', n, '" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , ZONETYPE = ORDERED, I = ', posic%numbnodes, ' , C = BLUE , DATAPACKING = POINT'
              DO i=1,posic%numbnodes ! loop over structural nodes
                WRITE(244,600) coora(1:3,posic%nodes(i)), velnp(1:3,posic%nodes(i)) 
              ENDDO
              
              ! CONTROL POINTS
              ! determine some header data
              nRegD = 0 ! number of aero nodes with regular dependency
              nIrrD = 0 ! number of aero nodes with irregular dependency
              DO i=1,posic%np
                IF(posic%cpReguDep(i))THEN
                  nRegD = nRegD+1
                ELSE
                  nIrrD = nIrrD+1
                ENDIF
              ENDDO
              
              ! regular dependency
              IF(nRegD.NE.0)THEN
                strandID = strandID + 1
                WRITE( 244 , * ) 'ZONE T = "CPsRegDep - ', TRIM(posic%m_est),' n=', n, '" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , ZONETYPE = ORDERED, I = ', nRegD, ' , C = GREEN , DATAPACKING = POINT'
                DO i=1,posic%np
                  IF( posic%cpReguDep(i) )THEN
                    WRITE(244,600) posic%panel(i)%xyzcp(1:3), posic%panel(i)%velcp(1:3)
                  ENDIF
                ENDDO
              ENDIF
              
              ! irregular dependency
              IF(nIrrD.NE.0)THEN
                strandID = strandID + 1
                WRITE( 244 , * ) 'ZONE T = "CPsIrrDep - ', TRIM(posic%m_est),' n=', n, '" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , ZONETYPE = ORDERED, I = ', nIrrD, ' , C = RED , DATAPACKING = POINT'
                DO i=1,posic%np
                  IF( .NOT.posic%cpReguDep(i) )THEN
                    WRITE(244,600) posic%panel(i)%xyzcp(1:3), posic%panel(i)%velcp(1:3)
                  ENDIF
                ENDDO
              ENDIF
              
              WRITE( 244 , 400 )
              
            ENDIF
            
            IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
              posic => posic%next
            ELSE
              exit
            ENDIF
	    	
          ENDDO
        ENDIF
        
        CLOSE(UNIT = 244)
        
        400 FORMAT(/) ! new line
        600 FORMAT(3(1X,E12.5E2),3(1X,E12.5E2)) ! 3 times coords. plus 3 times vel. comps.
        
      ENDSUBROUTINE nbstOmTecplot
    
  ENDSUBROUTINE inter_kin
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_loa
    
    ! INTERaction_LOAds
  
    ! Mauro S. Maza - 11/10/2011
    !                 04/01/2016
    
    ! Determines forces acting on structural nodes
    ! from those acting on control points of lifting surfaces
    
    USE StructuresA,    ONLY: nsec, npan
    
    IMPLICIT NONE
    
    TYPE(pair), POINTER     :: posic
    
    INTEGER (kind=4)        :: n1, n2, nn6(6), nn3(3)
    INTEGER                 :: i, j, m, fp
    
    REAL(kind = 8)          :: F(3), F1(3), F2(3),           &
                               d1l(3), d2l(3), d1(3), d2(3), &
                               m1l(3), m2l(3), m1(3), m2(3), &
                               lambda1(3,3), lambda2(3,3),   &
                               Fa(3,1), FF6(3,6), FF3(3,3),  &
                               fces(9), alphas(1,6)
    REAL                    :: xita1
    LOGICAL, PARAMETER :: printMore01=.FALSE.
    
    loa(:,:) = 0
    
    IF( ASSOCIATED(headpair) )THEN    ! Check if a list is empty
    posic => headpair
      
    DO ! Loop over blade pairs
        
        IF( TRIM(posic%m_aer)=='blade1' .OR. &      ! no forces for non-lifting surfaces
            TRIM(posic%m_aer)=='blade2' .OR. &
            TRIM(posic%m_aer)=='blade3'        )THEN
          
          SELECT CASE (TRIM(posic%eType))
          CASE('BEAM') ! ---------------------
            DO i=1,nsec   ! loop over sections in lifting surface
              fp = posic%section(i)%panels(1)  ! section first panel
              DO j=1,npan ! loop over panels in present section
                ! General data
                m = fp - 1 + j    ! panel (control point) index
                n1 = posic%cpmnods(1,m)     ! master node n1
                n2 = posic%cpmnods(2,m)     ! master node n2
                d1l = posic%cpdists(1:3,m)  ! Initial distance vector in local coords. of master node n1
                d2l = posic%cpdists(4:6,m)  ! Initial distance vector in local coords. of master node n2
                lambda1 = RESHAPE(euler(1:9,n1),(/3,3/))	! Rotation matrix for n1 local system
                lambda2 = RESHAPE(euler(1:9,n2),(/3,3/))	! Rotation matrix for n2 local system
                d1 = MATMUL(lambda1,d1l)  ! Distance vector from master node to aero. node n1-->an
                d2 = MATMUL(lambda2,d2l)  ! Distance vector from master node to aero. node n2-->an
                xita1 = posic%cpxita(1,m)         ! coordinate for projection (relative to n1)
                F = posic%panel(m)%CF * CF2F  ! Force, normal to panel m, acting on its control point
                F1 = F * (1-xita1)
                F2 = F *    xita1
                ! Forces
                loa(1:3,n1) = loa(1:3,n1) + F1   ! contribution of aero. node m force to struc. node n1 force
                loa(1:3,n2) = loa(1:3,n2) + F2   ! contribution of aero. node m force to struc. node n2 force
                ! Moments
                CALL cross_product( m1 , d1 , F1 )  ! Moment applied on master node due to translation of force, global coordinates
                CALL cross_product( m2 , d2 , F2 )  ! Moment applied on master node due to translation of force, global coordinates
                m1l = MATMUL(m1,lambda1) ! Moment in local coordinates of master node
                m2l = MATMUL(m2,lambda2) ! Moment in local coordinates of master node
                loa(4:6,n1) = loa(4:6,n1) + m1l
                loa(4:6,n2) = loa(4:6,n2) + m2l
                
                if(printMore01) print '(I4,3(xF15.2))', m, F
              ENDDO
            ENDDO
            
          CASE('NBST') ! ---------------------
            DO i=1,nsec*npan
              F = posic%panel(i)%CF*CF2F
              IF( posic%cpReguDep(i) )THEN
                Fa = RESHAPE( F, [3,1] )
                nn6 = posic%cpmnods(1:6,i)
                ! note that the order of cpdists(1:2,i) is inverted because they act as weights in the sum
                alphas = RESHAPE( [ posic%cpdists(2,i)*posic%cpxita(1:3,i) , posic%cpdists(1,i)*posic%cpxita(4:6,i) ] , [1,6] )
                FF6 = MATMUL( Fa(:,:) , alphas(:,:) )
                loa(1:3,nn6) = loa(1:3,nn6) + FF6
              ELSE
                fces(1:9) = MATMUL(TRANSPOSE(RESHAPE(nbstH(1:27,posic%cpmnods(4,i)),[3,9])),F)
                nn3 = posic%cpmnods(1:3,i)
                FF3 = RESHAPE(fces,[3,3])
                loa(1:3,nn3) = loa(1:3,nn3) + FF3
              ENDIF
            ENDDO
            
          ENDSELECT
            
        ENDIF
        
    IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
      posic => posic%next
    ELSE
      exit
    ENDIF
    
      ENDDO
    ENDIF
    
    !CALL nbstFcsTecplot   ! for debugging purposes only
    !CALL FcsTecplotCurves ! for debugging purposes only
    CALL AeroFcsOnStruc   ! for debugging purposes only
    
  CONTAINS
  
      ! ------------------------------- debug -------------------------------
      ! =====================================================================
      
      SUBROUTINE nbstFcsTecplot
        
        ! prints Tecplot file with forces on control points and struct. nodes
        
        USE ctrl_db,        ONLY:   ttime
        
        IMPLICIT NONE
        
        ! internal vars
        LOGICAL                             ::  fExist
        INTEGER                             ::  strandID, i, nRegD, nIrrD
        
        
        INQUIRE(FILE='nbstFcsTecplot.tec', exist=fExist)
        IF(fExist)THEN
          OPEN(UNIT = 245, FILE='nbstFcsTecplot.tec', STATUS='old', POSITION='append', ACTION='write')
        ELSE
          OPEN(UNIT = 245, FILE='nbstFcsTecplot.tec', STATUS='new', ACTION='write')
          ! Heading
          WRITE( 245 , * ) 'TITLE = "Forces on structural nodes and control points"'
          WRITE( 245 , * ) 'VARIABLES = "X" "Y" "Z" "f_X" "f_Y" "f_Z"'
          WRITE( 245 , 400 )
        ENDIF
        
        strandID = 0
        
        IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	      posic => headpair
          DO ! Loop over pairs
            IF( TRIM(posic%eType)=='NBST' )THEN
              
              ! STRUCTURAL NODES
              strandID = strandID + 1
              WRITE( 245 , * ) 'ZONE T = "StruNodes - ', TRIM(posic%m_est),' n=', n-1, '" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , ZONETYPE = ORDERED, I = ', posic%numbnodes, ' , C = BLUE , DATAPACKING = POINT'
              DO i=1,posic%numbnodes ! loop over structural nodes
                WRITE(245,600) coora(:,posic%nodes(i)), loa(1:3,posic%nodes(i))
              ENDDO
              WRITE( 245 , 400 )
              
              ! CONTROL POINTS
              ! determine some header data
              nRegD = 0 ! number of aero nodes with regular dependency
              nIrrD = 0 ! number of aero nodes with irregular dependency
              DO i=1,posic%np
                IF(posic%cpReguDep(i))THEN
                  nRegD = nRegD+1
                ELSE
                  nIrrD = nIrrD+1
                ENDIF
              ENDDO
              
              ! regular dependency
              IF(nRegD.NE.0)THEN
                strandID = strandID + 1
                WRITE( 245 , * ) 'ZONE T = "CPs_RegularDep - ', TRIM(posic%m_est),' n=', n-1, '" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , ZONETYPE = ORDERED, I = ', nRegD, ' , C = GREEN , DATAPACKING = POINT'
                DO i=1,posic%np
                  IF( posic%cpReguDep(i) )THEN
                    WRITE(245,600) posic%panel(i)%xyzcp, posic%panel(i)%CF*CF2F
                  ENDIF
                ENDDO
              ENDIF
              WRITE( 245 , 400 )
                
              ! irregular dependency
              IF(nIrrD.NE.0)THEN
                strandID = strandID + 1
                WRITE( 245 , * ) 'ZONE T = "CPs_IrregularDep - ', TRIM(posic%m_est),' n=', n-1, '" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , ZONETYPE = ORDERED, I = ', nIrrD, ' , C = RED , DATAPACKING = POINT'
                DO i=1,posic%np
                  IF( .NOT.posic%cpReguDep(i) )THEN
                    WRITE(245,600) posic%panel(i)%xyzcp, posic%panel(i)%CF*CF2F
                  ENDIF
                ENDDO
              ENDIF
              WRITE( 245 , 400 )
              
            ENDIF
            
            IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
              posic => posic%next
            ELSE
              exit
            ENDIF
	    	
          ENDDO
        ENDIF
        
        CLOSE(UNIT = 245)
        
        400 FORMAT(/) ! new line
        600 FORMAT(3(1X,E12.5E2),3(1X,E12.5E2)) ! 3 times coords. plus 3 times vel. comps.
        
      ENDSUBROUTINE nbstFcsTecplot
      
      ! ---------------------------------------------------------------------
      
      SUBROUTINE FcsTecplotCurves
        
        ! prints file for XY-line comparison with Tecplot between total force
        ! and total moment in aero grid and struc mesh
        
        USE ctrl_db,        ONLY:   ttime
        
        IMPLICIT NONE
        
        ! internal vars
        LOGICAL, PARAMETER      ::  printMore1=.FALSE.
        LOGICAL                 ::  fExist
        INTEGER                 ::  i, ns
        REAL(8)                 ::  write03, write04, write05, write06
        REAL(8), DIMENSION(3)   ::  Fs, Ms, Fa, Ma, f, r, m
        REAL(8), ALLOCATABLE    ::  fArray(:,:), mArray(:,:), auxArr(:)
        
        
        INQUIRE(FILE='FcsTecplotCurves.tec', exist=fExist)
        IF(fExist)THEN
          OPEN(UNIT = 247, FILE='FcsTecplotCurves.tec', STATUS='old', POSITION='append', ACTION='write')
        ELSE
          OPEN(UNIT = 247, FILE='FcsTecplotCurves.tec', STATUS='new', ACTION='write')
          ! Heading
          WRITE( 247 , * ) '## for debugging purposes only'
          WRITE( 247 , * ) '## file is never deleted by the program'
          WRITE( 247 , * ) '## data is simply added to the end'
          WRITE( 247 , 400 )
          WRITE( 247 , * ) '## Designed for NBST element type'
          WRITE( 247 , * ) '## for wich there are no moments applied on structural nodes'
          WRITE( 247 , 400 )
          WRITE( 247 , * ) 'TITLE = "Total forces and moments comparison"'
          WRITE( 247 , * ) 'VARIABLES = "n" "t" "|Fs|/|Fa|" "cos(angF)" "|Ms|/|Ma|" "cos(angM)"'
          WRITE( 247 , 400 )
        ENDIF
        
        IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
          
          Fs = 0
          Ms = 0
          Fa = 0
          Ma = 0
	      
          posic => headpair
          DO ! Loop over pairs
            !IF( TRIM(posic%eType)=='NBST' )THEN
            !IF( TRIM(posic%m_aer)=='blade1' .OR. &      ! no forces for non-lifting surfaces
            !    TRIM(posic%m_aer)=='blade2' .OR. &
            !    TRIM(posic%m_aer)=='blade3'        )THEN
            IF( TRIM(posic%m_aer)=='blade2' )THEN
              
              ! STRUCTURAL NODES
              ALLOCATE(fArray(3,posic%numbnodes), &
                       mArray(3,posic%numbnodes), &
                       auxArr(  posic%numbnodes))
              
              
              if(printMore1) print*, 'loads on structural nodes'
              DO i=1,posic%numbnodes ! loop over structural nodes
                ns = posic%nodes(i)
                f = loa(1:3,ns)
                !Fs = Fs + f
                fArray(:,i) = f
                r = coora(1:3,ns)
                CALL cross_product(m,r,f) ! moment wrt global origin - cross_product(out, in1, in2)
                !Ms = Ms + m
                mArray(:,i) = m
                if(printMore1) print '(I4,6(xF15.2))', i, f, m
              ENDDO
              
              DO i=1,3
                auxArr = fArray(i,:)
                CALL sort(auxArr)
                fArray(i,:) = auxArr
                
                auxArr = mArray(i,:)
                CALL sort(auxArr)
                mArray(i,:) = auxArr
              ENDDO
              
              IF(posic%numbnodes.GT.1)THEN
                DO i=2,posic%numbnodes ! loop over structural nodes
                  fArray(:,1) = fArray(:,1) + fArray(:,i)
                  mArray(:,1) = mArray(:,1) + mArray(:,i)
                ENDDO
              ENDIF
              Fs = Fs + fArray(:,1)
              Ms = Ms + mArray(:,1)
              
              DEALLOCATE(fArray,mArray,auxArr)
              
              
              
              ! CONTROL POINTS
              ALLOCATE(fArray(3,nsec*npan), &
                       mArray(3,nsec*npan), &
                       auxArr(  nsec*npan))
              
              
              if(printMore1) print*, 'loads on control points'
              DO i=1,nsec*npan ! lifting surface's CPs
                f = posic%panel(i)%CF*CF2F
                !Fa = Fa + f
                fArray(:,i) = f
                r = posic%panel(i)%xyzcp
                CALL cross_product(m,r,f) ! moment wrt global origin - cross_product(out, in1, in2)
                !Ma = Ma + m
                mArray(:,i) = m
                if(printMore1) print '(I4,6(xF10.2))', i, f, m
              ENDDO
              
              DO i=1,3
                auxArr = fArray(i,:)
                CALL sort(auxArr)
                fArray(i,:) = auxArr
                
                auxArr = mArray(i,:)
                CALL sort(auxArr)
                mArray(i,:) = auxArr
              ENDDO
              
              IF(nsec*npan.GT.1)THEN
                DO i=2,nsec*npan ! lifting surface's CPs
                  fArray(:,1) = fArray(:,1) + fArray(:,i)
                  mArray(:,1) = mArray(:,1) + mArray(:,i)
                ENDDO
              ENDIF
              Fa = Fa + fArray(:,1)
              Ma = Ma + mArray(:,1)
              
              DEALLOCATE(fArray,mArray,auxArr)
              
            ENDIF
            
            IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
              posic => posic%next
            ELSE
              exit
            ENDIF
	    	
          ENDDO
          
          ! write data into file
          write03 = dnrm2(3,Fs,1)/dnrm2(3,Fa,1)
          write04 = DOT_PRODUCT(Fa,Fs) /dnrm2(3,Fa,1) /dnrm2(3,Fs,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
          write05 = dnrm2(3,Ms,1)/dnrm2(3,Ma,1)
          write06 = DOT_PRODUCT(Ma,Ms) /dnrm2(3,Ma,1) /dnrm2(3,Ms,1) ! MKL library routine - ?nrm2(n,x,s): euclidean norm of x(n) with stride s
          WRITE(247,800) n, ttime, write03, write04, write05, write06
          
        ENDIF
        
        CLOSE(UNIT = 247)
        
        400 FORMAT(/) ! new line
        800 FORMAT(x,I6, x,E11.5E2, 4(x,F11.8)) ! 1 int + 1 gen real + 4 reals around 1
        
      ENDSUBROUTINE FcsTecplotCurves
      
      ! ---------------------------------------------------------------------
      
      SUBROUTINE AeroFcsOnStruc
        
        ! prints file with aerodynamic loads on estructural (master) nodes
        ! of BLADE2
        
        USE ctrl_db,        ONLY:   ttime
        
        IMPLICIT NONE
        
        ! internal vars
        LOGICAL                 ::  fExist
        INTEGER                 ::  i, ns
        REAL(8)                 ::  write03, write04, write05, write06
        REAL(8), DIMENSION(3)   ::  Fs, Ms, Fa, Ma, f, r, m
        REAL(8), ALLOCATABLE    ::  fArray(:,:), mArray(:,:), auxArr(:)
        
        
        INQUIRE(FILE='AeroFcsOnStruc.dat', exist=fExist)
        IF(fExist)THEN
          OPEN(UNIT = 252, FILE='AeroFcsOnStruc.dat', STATUS='old', POSITION='append', ACTION='write')
        ELSE
          OPEN(UNIT = 252, FILE='AeroFcsOnStruc.dat', STATUS='new', ACTION='write')
          ! Heading
          WRITE( 252 , 400 )
          WRITE( 252 , * ) 'Aerodynamic loads on estructural (master) nodes'
          WRITE( 252 , 400 )
          WRITE( 252 , * ) 'FORMAT'
          WRITE( 252 , * ) 'iter time'
          WRITE( 252 , * ) 'numbNodes'
          WRITE( 252 , * ) 'nodeId fx fy fz mx my mz'
          WRITE( 252 , * ) 'nodeId fx fy fz mx my mz'
          WRITE( 252 , * ) '...'
          WRITE( 252 , * ) 'nodeId fx fy fz mx my mz'
          WRITE( 252 , 400 )
          WRITE( 252 , * ) 'Note: forces are in global coords., but'
          WRITE( 252 , * ) '      moments are in local nodes coords.'
          WRITE( 252 , 400 )
        ENDIF
        
        IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
          
          ! write data into file
          WRITE(252,600) n, ttime
          
          posic => headpair
          DO ! Loop over pairs
            IF( TRIM(posic%m_aer)=='blade1' .OR. &      ! no forces for non-lifting surfaces
                TRIM(posic%m_aer)=='blade2' .OR. &
                TRIM(posic%m_aer)=='blade3'        )THEN
            ! IF( TRIM(posic%m_aer)=='blade2' )THEN
              
              ! write data into file
              WRITE(252,500) posic%numbnodes
              DO i=1,posic%numbnodes ! loop over structural nodes
                ns = posic%nodes(i)
                f = loa(1:3,ns)
                m = loa(4:6,ns)
                ! write data into file
                WRITE(252,800) ns, f(1), f(2), f(3), m(1), m(2), m(3)
              ENDDO
              
            ENDIF
            
            IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
              posic => posic%next
            ELSE
              exit
            ENDIF
	    	
          ENDDO
          
          ! write data into file
          WRITE( 252 , 400 )
          
        ENDIF
        
        CLOSE(UNIT = 252)
        
        400 FORMAT(/) ! new line
        500 FORMAT(x,I6) ! 1 int
        600 FORMAT(x,I6, x,1(  E11.5E2)) ! 1 int + 1 gen real
        800 FORMAT(x,I6,   6(x,E12.5E2)) ! 1 int + 6 gen real
        
      ENDSUBROUTINE AeroFcsOnStruc
      
      ! ---------------------------------------------------------------------
      
      SUBROUTINE sort(arr)
        
        ! sorts in ascending order real rank-1 arrays
        
        ! Mauro S. Maza - 06/01/2016
        
        IMPLICIT NONE
        
        INTEGER                                             ::  sizeArr, i, ind
        REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(inout)   ::  arr
        
        
        sizeArr = SIZE(arr)
        DO i=1,sizeArr
          ind = MINLOC( arr(i:sizeArr), 1 ) + i - 1
          CALL swap(arr(i),arr(ind))
        ENDDO
        
      ENDSUBROUTINE sort
      
      ! ---------------------------------------------------------------------
      
      SUBROUTINE swap(a,b)
        
        IMPLICIT NONE
        
        REAL(8), INTENT(INOUT)  ::  a, b
        REAL(8)                 ::  Temp
        
        
        Temp = a
        a    = b
        b    = Temp
        
      ENDSUBROUTINE swap
      
  ENDSUBROUTINE inter_loa
  
ENDMODULE inter_sr

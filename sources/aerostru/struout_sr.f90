
MODULE struout_sr
  
  USE struout_db
  USE introut_db,       ONLY:   newTPFileS
  USE inter_db,         ONLY:   pair,       &   ! types
                                n, headpair     !variables
  
  IMPLICIT NONE
  
CONTAINS ! =============================================
    
    
    ! Subroutines for output
    !   inter_tp6           Structural mesh in TecPlot output
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_tp6
    
    ! INTERaction_TecPlot6
  
    ! Mauro S. Maza - 18/12/2012
    !                 02/09/2015
    
    ! CUIDADO - Para los cuerpos rígidos produce un elemento de viga de dos
    ! nodos con las mismas coordenadas, por lo que sólo se ve UN nudo maestro
    
    ! Prints structural mesh in TecPlot output
    
    USE ctrl_db,        ONLY:   ttime, npoin
    USE npo_db,         ONLY:   coora, label
    USE ele13_db,       ONLY:   srch_ele13,			 & ! routines
                                ele13_set, ele13,    & ! types
                                head                   ! variables
    USE aeroout_db,     ONLY:   aerGriOut
    
    IMPLICIT NONE
    
    ! internal vars
    LOGICAL                             ::  found
    INTEGER                             ::  i, n1, n2, n3, strandID, countb, lBou, uBou
    INTEGER, ALLOCATABLE, DIMENSION(:)  ::  nodsIndxs
    !INTEGER, POINTER                    ::  nodes(:) ! subroutine NODS13 dummy argument is pointer
    REAL                                ::  xyz_e1(3), xyz_e2(3)
    TYPE(pair),         POINTER         ::  posic
    TYPE(ele13_set),    POINTER         ::  anterset, set
    TYPE(ele13),        POINTER         ::  cEle
    
    
    ! Struc. Meshes always after Aero. Grids and before Panel Forces
    IF( aerGriOut )THEN
      strandID = 9 ! 9 aero. grid zones per iteration
    ELSE
      strandID = 0
    ENDIF
    
    IF( newTPFileS )THEN ! first step for current file/strategy
      
      newTPFileS = .FALSE. ! do not write connectivities again for this file
      
      IF( ASSOCIATED(headpair) )THEN  ! Check if a list is empty
        posic => headpair
        DO ! Loop over pairs
          strandID = strandID + 1
          
          SELECT CASE (TRIM(posic%eType))
          CASE('RIGID')
            WRITE( 237 , * ) 'ZONE T = "STRC - ', TRIM(posic%m_est),' n=',n,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ', 2, ' , E = ', 1, ' , ET = LINESEG , C = BLACK , PASSIVEVARLIST=[4-7]'
            n1 = posic%mnods(1,1)
            xyz_e1 = coora(:,n1)
            WRITE( 237 , 100 ) xyz_e1(1), xyz_e1(2), xyz_e1(3)
            WRITE( 237 , 100 ) xyz_e1(1), xyz_e1(2), xyz_e1(3)
            WRITE( 237 , 400 )
            WRITE( 237 , * ) 1 , 2
            WRITE( 237 , 400 )
            
          CASE('BEAM')
            WRITE( 237 , * ) 'ZONE T = "STRC - ', TRIM(posic%m_est),' n=',n,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ', 2*posic%nelem, ' , E = ', posic%nelem, ' , ET = LINESEG , C = BLACK , PASSIVEVARLIST=[4-7]'
            DO i=1,posic%nelem
              n1 = posic%lnods(1,i)
              n2 = posic%lnods(2,i)
              xyz_e1 = coora(:,n1)
              xyz_e2 = coora(:,n2)
              WRITE( 237 , 100 ) xyz_e1(1) , xyz_e1(2) , xyz_e1(3)
              WRITE( 237 , 100 ) xyz_e2(1) , xyz_e2(2) , xyz_e2(3)
            ENDDO
            WRITE( 237 , 400 )
            DO i=1,posic%nelem
              WRITE( 237 , 200 ) 2*i-1 , 2*i
            ENDDO
            WRITE( 237 , 400 )
            
          CASE('NBST')
            CALL srch_ele13(head, anterSet, set, posic%m_est, found) ! outputs SET pointing to set named elmSetName=posic%m_est
                                                                     ! in the list of NBST elements (whose first link is HEAD)
            lBou = MINVAL(posic%nodes,1)
            uBou = MAXVAL(posic%nodes,1)
            ALLOCATE( nodsIndxs(lBou:uBou) )
            nodsIndxs = 0
            ! zone header
            WRITE( 237 , * ) 'ZONE T = "STRC - ', TRIM(posic%m_est),' n=',n,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ', posic%numbnodes, ' , E = ', posic%nelem, ' , ET = TRIANGLE , C = BLACK , PASSIVEVARLIST=[4-7]'
            ! coords... and more
            DO i=1,posic%numbnodes
              n1 = posic%nodes(i)
              nodsIndxs(n1) = i ! inverted relation between node labes and indexes of that expressed in posic%nodes
              xyz_e1 = coora(:,n1)
              WRITE( 237 , 100 ) xyz_e1(1) , xyz_e1(2) , xyz_e1(3)
            ENDDO
            WRITE( 237 , 400 )
            ! connectivities
            cEle => set%head ! Current ELEment initially points to head of list
            DO i=1,posic%nelem
              ! indexes of connectivity nodes in TecPlot zone
              WRITE( 237 , 300 ) nodsIndxs(cEle%lnods(1:3))
              cEle => cEle%next
            ENDDO
            DEALLOCATE( nodsIndxs )
            WRITE( 237 , 400 )
            
          ENDSELECT
          
          IF (ASSOCIATED (posic%next)) THEN    ! Check if it's not the end of the list
            posic => posic%next
          ELSE
            exit
          ENDIF
        
        ENDDO
      ENDIF
      
    ELSE ! not the first step - use CONNECTIVITYSHAREZONE
      countb = 0     ! counter for pairs in a step - for CONNECTIVITYSHAREZONE
      IF( ASSOCIATED(headpair) )THEN       ! Check if a list is empty
        posic => headpair
        DO ! Loop over pairs
          countb = countb + 1
          strandID = strandID + 1
          
          SELECT CASE (TRIM(posic%eType))
          CASE('RIGID')
            WRITE( 237 , * ) 'ZONE T = "STRC - ', TRIM(posic%m_est),' n=',n,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ', 2, ' , E = ', 1, ' , ET = LINESEG , C = BLACK , CONNECTIVITYSHAREZONE = ', countb, ' , PASSIVEVARLIST=[4-7]'
            n1 = posic%mnods(1,1)
            xyz_e1 = coora(:,n1)
            WRITE( 237 , 100 ) xyz_e1(1), xyz_e1(2), xyz_e1(3)
            WRITE( 237 , 100 ) xyz_e1(1), xyz_e1(2), xyz_e1(3)
            WRITE( 237 , 400 )
              
          CASE('BEAM')
            WRITE( 237 , * ) 'ZONE T = "STRC - ', TRIM(posic%m_est),' n=',n,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ', 2*posic%nelem, ' , E = ', posic%nelem, ' , ET = LINESEG , C = BLACK , CONNECTIVITYSHAREZONE = ', countb, ' , PASSIVEVARLIST=[4-7]'
            DO i=1,posic%nelem
              n1 = posic%lnods(1,i)
              n2 = posic%lnods(2,i)
              xyz_e1 = coora(:,n1)
              xyz_e2 = coora(:,n2)
              WRITE( 237 , 100 ) xyz_e1(1) , xyz_e1(2) , xyz_e1(3)
              WRITE( 237 , 100 ) xyz_e2(1) , xyz_e2(2) , xyz_e2(3)
            ENDDO
            WRITE( 237 , 400 )
            
          CASE('NBST')
            WRITE( 237 , * ) 'ZONE T = "STRC - ', TRIM(posic%m_est),' n=',n,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ', posic%numbnodes, ' , E = ', posic%nelem, ' , ET = TRIANGLE , C = BLACK , CONNECTIVITYSHAREZONE = ', countb, ' , PASSIVEVARLIST=[4-7]'
            ! coords
            DO i=1,posic%numbnodes
              n1 = posic%nodes(i)
              xyz_e1 = coora(:,n1)
              WRITE( 237 , 100 ) xyz_e1(1) , xyz_e1(2) , xyz_e1(3)
            ENDDO
            WRITE( 237 , 400 )
            
          ENDSELECT
          
          IF (ASSOCIATED (posic%next)) THEN    ! Check if it's not the end of the list
            posic => posic%next
          ELSE
            exit
          ENDIF
        
        ENDDO
      ENDIF
      
    ENDIF
    
    100 FORMAT(3(1X,F9.3)) ! coords.
    200 FORMAT(2(1X,I5))   ! 2 node connects.
    300 FORMAT(3(1X,I6))   ! 3 node connects.
    400 FORMAT(/)          ! new line
  
  ENDSUBROUTINE inter_tp6
  
ENDMODULE struout_sr

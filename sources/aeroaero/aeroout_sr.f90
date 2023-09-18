
MODULE aeroout_sr
  
  USE aeroout_db
  USE ctrl_db,      ONLY:   ttime ! variables
  USE inter_db,     ONLY:   pair,             & ! types
                            headpair, n, CF2F   ! variables
  
  IMPLICIT NONE
  
  TYPE(pair),       POINTER	    :: posic
  
CONTAINS ! =============================================

    ! Subroutines for output
    !   inter_tp2           Aerodynamic grid in TecPlot output
    !   inter_tp4           Panels' forces data in TecPlot output
    !
    !   inter_aet           "Time for aerodynamic calculations" in root.AT.dat
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_tp2
    
    ! INTERaction_TecPlot2
  
    ! Mauro S. Maza - 18/12/2012
    !                 02/09/2015
        
    ! Prints aerodynamic grid in TecPlot output
    
    USE StructuresA, ONLY: twknode,                         & ! types
                           wk1nodes, wk2nodes, wk3nodes,    & ! variables
                           NSA,nnwke
    USE ArraysA,     ONLY: connectwk1, connectwk2, connectwk3 ! variables
    USE introut_db,  ONLY: newTPFileA
    
    IMPLICIT NONE
    
    TYPE(twknode), POINTER :: wkinodes(:)
    INTEGER                :: strandID, i, n_out, s, j, &
                              countb ! counter for pairs in a step
    INTEGER,       POINTER :: connectwki(:,:)
    
    
    strandID = 0 ! if aero. grids printed, always first zones
    
    ! Bound sheets
    IF( newTPFileA )THEN ! first step for current file/strategy
      
      newTPFileA = .FALSE. ! do not write connectivities again for this file
      
	  IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	    posic => headpair
        DO ! Loop over pairs
          
          strandID = strandID + 1
          WRITE( 232 , * ) 'ZONE T = "AERO - ', TRIM(posic%m_aer),' n=', n, '" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ', posic%nan, ' , E = ', posic%np, ' , ET = QUADRILATERAL , C = BLUE , PASSIVEVARLIST=[4-7]'
          DO i=1,posic%nan
            WRITE( 232 , 100 ) posic%nodsaer(i)%xyz(1), posic%nodsaer(i)%xyz(2), posic%nodsaer(i)%xyz(3)
          ENDDO
          WRITE( 232 , 400 )
          DO i=1,posic%np
            WRITE( 232 , 200 ) posic%panel(i)%nodes(1), posic%panel(i)%nodes(2), posic%panel(i)%nodes(3), posic%panel(i)%nodes(4)
          ENDDO
          WRITE( 232 , 400 )
          
          IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
            posic => posic%next
          ELSE
            exit
          ENDIF
	  	
        ENDDO
      ENDIF
      
    ELSE ! not the first step - use CONNECTIVITYSHAREZONE
      countb = 0     ! counter for pairs in a step - for CONNECTIVITYSHAREZONE
	  IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	    posic => headpair
        DO ! Loop over pairs
          
          countb = countb + 1
          strandID = strandID + 1
          WRITE( 232 , * ) 'ZONE T = "AERO - ', TRIM(posic%m_aer),' n=',n,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ', posic%nan, ' , E = ', posic%np, ' , ET = QUADRILATERAL , C = BLUE , CONNECTIVITYSHAREZONE = ', countb, ' , PASSIVEVARLIST=[4-7]'
          DO i=1,posic%nan
            WRITE( 232 , 100 ) posic%nodsaer(i)%xyz(1), posic%nodsaer(i)%xyz(2), posic%nodsaer(i)%xyz(3)
          ENDDO
          WRITE( 232 , 400 )
          
          IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
            posic => posic%next
          ELSE
            exit
          ENDIF
	  	
        ENDDO
      ENDIF
      
    ENDIF
    
    ! Wakes
    
    IF( n<NSA )THEN
      n_out = n
    ELSE
      n_out = NSA
    ENDIF
    
    IF( n/=0 )THEN
      DO i=1,3
        strandID = strandID + 1
        SELECT CASE (i) ! determines variables to use and prints zones' titles
        CASE (1)
          connectwki => connectwk1
          wkinodes => wk1nodes
          s = 0
          DO j=1,( nnwke-1) * (n_out)
            IF( connectwki(j,5)==1 )THEN
              s = s + 1
            ENDIF
          ENDDO
          WRITE( 232 , * ) 'ZONE T = " WAKE 1 n=', n ,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ' , (nnwke)*(n_out+1) , ', E = ' , s, ', ET = QUADRILATERAL , C = RED , PASSIVEVARLIST=[4-7]'
          WRITE( 232 , 400 )
        CASE (2)
          connectwki => connectwk2
          wkinodes => wk2nodes
          s = 0
          DO j=1,( nnwke-1) * (n_out)
            IF( connectwki(j,5)==1 )THEN
              s = s + 1
            ENDIF
          ENDDO
          WRITE( 232 , * ) 'ZONE T = " WAKE 2 n=', n ,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ' , (nnwke)*(n_out+1) , ', E = ' , s, ', ET = QUADRILATERAL , C = YELLOW , PASSIVEVARLIST=[4-7]'
          WRITE( 232 , 400 )
        CASE (3)
          connectwki => connectwk3
          wkinodes => wk3nodes
          s = 0
          DO j=1,( nnwke-1) * (n_out)
            IF( connectwki(j,5)==1 )THEN
              s = s + 1
            ENDIF
          ENDDO
          WRITE( 232 , * ) 'ZONE T = " WAKE 3 n=', n ,'" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , F = FEPOINT, N = ' , (nnwke)*(n_out+1) , ', E = ' , s, ', ET = QUADRILATERAL , C = GREEN , PASSIVEVARLIST=[4-7]'
          WRITE( 232 , 400 )
        ENDSELECT
        
        ! print coordinates
        DO j=1,(nnwke) * (n_out+1)
          WRITE( 232 , 100 ) wkinodes(j)%xyz(1), wkinodes(j)%xyz(2), wkinodes(j)%xyz(3)
        ENDDO
        WRITE( 232 , 400 )
        ! print connectivities
        DO j=1,(nnwke-1) * (n_out)
          IF( connectwki(j,5)==1 )THEN
            WRITE( 232 , 200 ) connectwki(j,1), connectwki(j,2), connectwki(j,3), connectwki(j,4)
          ENDIF
        ENDDO
        WRITE( 232 , 400 )
      ENDDO
    ELSE ! no wakes in first iteration
      strandID = strandID + 1
      WRITE( 232 , * ) 'ZONE T = " WAKE 1 n=           0 " , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , C = RED , PASSIVEVARLIST=[1-7]'
      strandID = strandID + 1
      WRITE( 232 , * ) 'ZONE T = " WAKE 2 n=           0 " , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , C = YELLOW , PASSIVEVARLIST=[1-7]'
      strandID = strandID + 1
      WRITE( 232 , * ) 'ZONE T = " WAKE 3 n=           0 " , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , C = GREEN , PASSIVEVARLIST=[1-7]'
    ENDIF
    
    100 FORMAT(3(1X,F9.3)) ! coords.
    200 FORMAT(4(1X,I5))   ! connects.
    400 FORMAT(/)          ! new line
  
  ENDSUBROUTINE inter_tp2
    
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_tp4
    
    ! INTERaction_TecPlot4
  
    ! Mauro S. Maza - 26/09/2012
    !                 02/09/2015
    
    ! Prints panels' forces data in TecPlot output
    
    USE StructuresA, ONLY: nsec, npan ! variables
    USE struout_db,  ONLY: strMesOut
    
    IMPLICIT NONE
    
    INTEGER                :: strandID,i, j, m, fp
    
    
    ! Panel Forces always after Aero. Grids and Struc. Meshes
    strandID = 0
    IF( aerGriOut )THEN
      strandID = strandID + 9 ! 9 aero. grid zones per iteration
    ENDIF
    IF( strMesOut )THEN
      strandID = strandID + 6 ! 6 struct. mesh zones per iteration
    ENDIF
    
    IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	  posic => headpair
      DO ! Loop over pairs
        IF( TRIM(posic%m_aer)=='blade1' .OR. &      ! no forces for non-lifting surfaces
            TRIM(posic%m_aer)=='blade2' .OR. &
            TRIM(posic%m_aer)=='blade3'        )THEN
          
          strandID = strandID + 1
          WRITE( 235 , * ) 'ZONE T = "PanelFRCES - ', TRIM(posic%m_aer),' n=', n-1, '" , STRANDID = ', strandID, ' , SOLUTIONTIME = ', ttime, ' , ZONETYPE = ORDERED, I = ', nsec*npan, ' , C = PURPLE , PASSIVEVARLIST=[7] , DATAPACKING = POINT'
          DO i=1,nsec   ! loop over sections in lifting surface
            fp = posic%section(i)%panels(1)  ! section first panel
            DO j=1,npan ! loop over panels in present section
              m = fp - 1 + j
              
              WRITE( 235 , 600 ) posic%panel(m)%xyzcp(1),   posic%panel(m)%xyzcp(2),   posic%panel(m)%xyzcp(3), & ! control point coordinates
                                 posic%panel(m)%CF(1)*CF2F, posic%panel(m)%CF(2)*CF2F, posic%panel(m)%CF(3)*CF2F   ! cartesian panel force components
              
            ENDDO
          ENDDO
          WRITE( 235 , 400 )
          
        ENDIF
        
        IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
          posic => posic%next
        ELSE
          exit
        ENDIF
		
      ENDDO
    ENDIF
    
    400 FORMAT(/) ! new line
    600 FORMAT(3(1X,F9.3),3(1X,E10.3)) ! 3 times coords. plus 3 times force comps.
  
  ENDSUBROUTINE inter_tp4
    
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_aet
    
    ! INTERaction_AEroTime
    
    ! Mauro S. Maza - 18/12/2012
    
    ! Prints "Time for aerodynamic calculations" in root.AT.dat
    
    USE outp_db, ONLY: time
    
    IMPLICIT NONE
    
    
    WRITE( 234 , * ) n, time(20)
    
  
  ENDSUBROUTINE inter_aet
  
ENDMODULE aeroout_sr

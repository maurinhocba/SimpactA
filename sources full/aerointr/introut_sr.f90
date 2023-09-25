
MODULE introut_sr
  
  USE introut_db
  USE aeroout_db,   ONLY:   aerGriOut, forcesOut
  USE struout_db,   ONLY:   strMesOut
  USE param_db,     ONLY:   milb, mlen, mlin, mich ! Global parameters
  USE name_db,      ONLY:   output
  USE inter_db,     ONLY:   pair,             & ! types
                            headpair            ! variables
  
  IMPLICIT NONE
  
    
CONTAINS ! =============================================
  
  ! ---------------- Subroutines for output -----------------
  ! =========================================================
    
    ! Subroutines for output
    !   introut_mng         MANAGE
    !
    !   projProbl           print messages
    !
    !   inter_tp1           OPEN
    !   inter_tp3           CLOSE some
    !
    !   inter_tp5 (Tecplot) Layout files        - opened, written and closed in one shot
    !
    !   inter_fam           "Blades Forces and Effective Moment" in root.FM.dat
    !                       (and "Generator Registry" in root.GE.dat)-NOT NOW
    !   inter_reg           "Aeroelastic Registry" in root.AR.dat
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE introut_mng(doWhat,actio)
    
    ! INTERaction OUTput_MaNaGe
    
    ! Mauro S. Maza - 20/03/2015
    !                 18/09/2015
    
    ! Manages output
    
    USE aeroout_sr,     ONLY:   inter_tp2, inter_tp4, inter_aet
    USE struout_sr,     ONLY:   inter_tp6
    
    IMPLICIT NONE
    
    INTEGER                         :: j
    CHARACTER(len=*), INTENT(IN)    :: doWhat
    CHARACTER(len=milb), OPTIONAL   :: actio
    CHARACTER(len=1)                :: ch       !character
    
    
    IF( PRESENT(actio) )THEN
      IF( TRIM(actio)=='NEW' .OR. TRIM(actio)=='RESTAR' )THEN
        DO j=1,LEN_TRIM(output)
          ch = output(j:j)
          IF (ch == '.') THEN 
            EXIT
          ENDIF
          outputRoot(j:j) = ch
          outputRootLen = j
        ENDDO
      ENDIF
    ENDIF
    
    
    SELECT CASE (doWhat)
    CASE('repProbl')
      CALL projProbl
    
    CASE('initial')
      ! Max duration of strategy not to exceed max number of zones Tecplot can manage
      zonesPerIter = 0
      IF( aerGriOut )THEN
        zonesPerIter = zonesPerIter + 9
      ENDIF
      IF( strMesOut )THEN
        zonesPerIter = zonesPerIter + 6
      ENDIF
      IF( forcesOut )THEN
        zonesPerIter = zonesPerIter + 3
      ENDIF
      ! Start some as 0 - done also in introut_db
      dtimeOld = 0
      DeltatOld = 0
      
    CASE('openFiles')
      IF( TRIM(actio)/='NEW' .AND. TRIM(actio)/='RESTAR' )THEN ! neither new problem nor restart, then new strategy
        CALL inter_tp3      ! close output files
      ENDIF
      CALL inter_tp1(actio) ! Creates output files and writes headers
    
    CASE('gridAndMesh')
      IF( aerGriOut )THEN
        CALL inter_tp2
      ENDIF
      IF( strMesOut )THEN
        CALL inter_tp6
      ENDIF
    
    CASE('panelForces')
      IF( forcesOut )THEN
        CALL inter_tp4  ! Print panel forces
      ENDIF
      
    CASE('frcesMomets')
      CALL inter_fam ! "Blades Forces and Effective Moment" in root.FM.dat
                     !NOT NOW - and "Generator Registry" in root.GE.dat
      
    CASE('aerodynTime')
      CALL inter_aet ! Print "Time for aerodynamic calculations" in root.AT.dat
      
    CASE('dtime')
      CALL inter_reg('dtime')
    
    CASE('end')
      CALL inter_tp3      ! Closes "last strategy's TecPlot output file"
      CLOSE( UNIT = 233 ) !    and "total forces and moments aplied over each blade"
      CLOSE( UNIT = 234 ) !    and "time for aerodynamic calculations"
      CLOSE( UNIT = 238 ) !    and "aeroelastic registry"
      CLOSE( UNIT = 239 ) !    and "generator registry"
      CALL inter_tp5('tec') ! Generates Tecplot .lay file
      CALL inter_tp5('plt') ! Generates Tecplot .lay file
      
    ENDSELECT
    
  ENDSUBROUTINE introut_mng
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE projProbl
    
    ! Mauro S. Maza - 18/09/2015
    
    ! Report problem in aeroprobls.dat file
    
    IMPLICIT NONE
    
    ! dummy args
    
    
    ! internal vars
    LOGICAL            :: fExist
    CHARACTER(len=16)  :: date
    
    
    INQUIRE(FILE='aeroprobls.dat', exist=fExist)
    IF(fExist)THEN
      OPEN(UNIT = 241, FILE='aeroprobls.dat', STATUS='old', POSITION='append', ACTION='write')
    ELSE
      OPEN(UNIT = 241, FILE='aeroprobls.dat', STATUS='new', ACTION='write')
    ENDIF
    CALL get_DDMonYY(date)
    WRITE(241, '(A)', advance='no' ) TRIM(date//' - '//msg)
    CLOSE(UNIT = 241)
    
  CONTAINS
  
      SUBROUTINE get_DDMonYY(date)
        CHARACTER(len=16), INTENT(out) :: date
      
        CHARACTER(len=2) :: dd
        CHARACTER(len=3) :: mons(12)
        CHARACTER(len=4) :: yyyy
        CHARACTER(len=2) :: hh
        CHARACTER(len=2) :: mm
        CHARACTER(len=2) :: ss
        INTEGER :: values(8)
      
        mons = ['Jan','Feb','Mar','Apr','May','Jun',&
          'Jul','Aug','Sep','Oct','Nov','Dec']
      
        CALL DATE_AND_TIME(VALUES=values)
      
        WRITE(  dd,'(i2)') values(3)
        WRITE(yyyy,'(i4)') values(1)
        WRITE(  hh,'(i2)') values(5)
        WRITE(  mm,'(i2)') values(6)
        WRITE(  ss,'(i2)') values(7)
      
        date = dd//mons(values(2))//yyyy(3:4)//'-'//hh//':'//mm//':'//ss
      ENDSUBROUTINE get_DDMonYY
    
  ENDSUBROUTINE projProbl
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_tp1(actio)
    
    ! INTERaction_TecPlot1 ! historical reasons
    
    ! Mauro S. Maza - 18/12/2012
    !                 02/09/2015
    
    ! Creates output files and writes headers
    ! .OR.
    ! Opens old output files for restart
    
    IMPLICIT NONE
    
    INTEGER             :: countBlades
    CHARACTER(len=mlin) :: fileName
    CHARACTER(len=milb) :: actio
    TYPE(pair), POINTER :: posic
    
    
    ! AERODYNAMIC GRID
    IF( aerGriOut )THEN
      fileName = TRIM(output)//TRIM('.Aero.tec')
      IF( TRIM(actio)=='RESTAR' )THEN  ! append data to end of file
        ! just open existing file
        OPEN( UNIT = 232 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'append' )
        WRITE( 232 , * ) '## FILE RESTARTED FROM HERE'
      ELSE ! new strategy file
        ! Create file
        OPEN( UNIT = 232 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'sequential' )
        ! Heading
        WRITE( 232 , * ) 'TITLE = "Aerodynamic Grid"'
        WRITE( 232 , * ) 'VARIABLES = "X" "Y" "Z" "n(X)" "n(Y)" "n(Z)" "Cp"'
        WRITE( 232 , 400 )
        ! Re-write connectivities
        newTPFileA = .TRUE.
      ENDIF
    ENDIF
    
    ! STRUCTURAL MESH
    IF( strMesOut )THEN
      fileName = TRIM(output)//TRIM('.Stru.tec')
      IF( TRIM(actio)=='RESTAR' )THEN  ! append data to end of file
        ! just open existing file
        OPEN( UNIT = 237 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'append' )
        WRITE( 237 , * ) '## FILE RESTARTED FROM HERE'
      ELSE ! new strategy file
        ! Create file
        OPEN( UNIT = 237 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'sequential' )
        ! Heading
        WRITE( 237 , * ) 'TITLE = "Structural Mesh"'
        WRITE( 237 , * ) 'VARIABLES = "X" "Y" "Z" "n(X)" "n(Y)" "n(Z)" "Cp"'
        WRITE( 237 , 400 )
        ! Re-write connectivities
        newTPFileS = .TRUE.
      ENDIF
    ENDIF
    
    ! FORCES IN PANELS
    IF( forcesOut )THEN
      fileName = TRIM(output)//TRIM('.PFcs.tec')
      IF( TRIM(actio)=='RESTAR' )THEN  ! append data to end of file
        ! just open existing file
        OPEN( UNIT = 235 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'append' )
        WRITE( 235 , * ) '# FILE RESTARTED FROM HERE'
      ELSE ! new strategy file
        ! Create file
        OPEN( UNIT = 235 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'sequential' )
        ! Heading
        WRITE( 235 , * ) 'TITLE = "Panels'' Forces"'
        WRITE( 235 , * ) 'VARIABLES = "X" "Y" "Z" "n(X)" "n(Y)" "n(Z)" "Cp"'
        WRITE( 235 , 400 )
      ENDIF
    ENDIF
    
    ! FORCES AND MOMENTS     and
    ! AERODYNAMIC TIME       and
    ! AERODYNAMIC REGISTRY   and
    ! GENERATOR REGISTRY
    IF( TRIM(actio)=='NEW' )THEN ! Create files and write headers
      
      ! FORCES AND MOMENTS --------------------------------------------------
      fileName = outputRoot(1:outputRootLen)//'.FM.dat'
      OPEN( UNIT = 233 , FILE = TRIM(fileName), STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'sequential' )
      WRITE( 233 , 400 )
      WRITE( 233 , * ) '  Aerodynamic forces applied on each blade'
      WRITE( 233 , * ) '  ----------------------------------------'
      WRITE( 233 , * ) '  (Blade order is the one used in the definition of pairs)'
      WRITE( 233 , 525 , advance='no' ) '                            '
      
      IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	    posic => headpair
      
        countBlades = 0
        DO ! Loop over pairs
          
          IF( TRIM(posic%m_aer)=='blade1' .OR. &      ! no forces for non-lifting surfaces
              TRIM(posic%m_aer)=='blade2' .OR. &
              TRIM(posic%m_aer)=='blade3'        )THEN
            
            countBlades = countBlades + 1
            IF( countBlades==3 )THEN
              WRITE( 233 , 510 , advance='no' ) TRIM(posic%m_aer), '                                                                     '
            ELSE
              WRITE( 233 , 525 , advance='no' ) TRIM(posic%m_aer), '                                                                     '
            ENDIF
            
          ENDIF
          
          IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
		    posic => posic%next
		  ELSE
		    exit
		  ENDIF
        
        ENDDO
      ENDIF
      WRITE( 233 , 550 , advance='no' ) '   t                        f_x                      f_y                      f_z                      f_x                      f_y                      f_z                      f_x                      f_y                      f_z'
      
      ! AERODYNAMIC TIME ----------------------------------------------------
      fileName = outputRoot(1:outputRootLen)//'.AT.dat'
      OPEN( UNIT = 234 , FILE = TRIM(fileName), STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'sequential' )
      WRITE( 234 , 400 )
      WRITE( 234 , * ) 'Time invested in aerodynamic calculations'
      WRITE( 234 , * ) '(time after iteration end)'
      WRITE( 234 , 400 )
      WRITE( 234 , * ) '     iter #   t [s]'
      
      ! AEROELASTIC REGISTRY ------------------------------------------------
      fileName = outputRoot(1:outputRootLen)//'.AR.dat'
      OPEN( UNIT = 238 , FILE = TRIM(fileName), STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'sequential' )
      WRITE( 238 , 400 )
      WRITE( 238 , * ) 'AEROELASTIC REGISTRY'
      WRITE( 238 , 400 )
      WRITE( 238 , * ) 'PROBLEM: ', outputRoot(1:outputRootLen)
      WRITE( 238 , 400 )
      CALL inter_reg('init')
      
      ! GENERATOR REGISTRY --------------------------------------------------
      fileName = outputRoot(1:outputRootLen)//'.GE.dat'
      OPEN( UNIT = 239 , FILE = TRIM(fileName), STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'sequential' )
      WRITE( 239 , 400 )
      WRITE( 239 , * ) 'GENERATOR REGISTRY'
      WRITE( 239 , 400 )
      WRITE( 239 , * ) 'PROBLEM: ', outputRoot(1:outputRootLen)
      WRITE( 239 , 400 )
      WRITE( 239 , * ) 'Data in columns:'
      WRITE( 239 , * ) '  [01] iteration number'
      WRITE( 239 , * ) '  [02] time'
      WRITE( 239 , * ) '  [03] rotor angular position'
      WRITE( 239 , * ) '  [04] rotor angular velocity'
      WRITE( 239 , * ) '  [05] rotor angular acceleration'
      WRITE( 239 , * ) '  [06] damping coefficient'
      WRITE( 239 , * ) '  [07] equivalent viscous moment'
      WRITE( 239 , * ) '  [08] power absorbed by de generator'
      WRITE( 239 , * ) '  [09] applied brake moment'
      WRITE( 239 , * ) '  [10] target advance coefficient (lambda)'
      WRITE( 239 , * ) '  [11] target power coefficient (c_P-tar)'
      WRITE( 239 , * ) '  [12] target blades'' pitch angle (q_tar)'
      WRITE( 239 , * ) '  [13] q1 - pith angle for balde 1'
      WRITE( 239 , * ) '  [14] q2 - pith angle for balde 2'
      WRITE( 239 , * ) '  [15] q3 - pith angle for balde 3'
      WRITE( 239 , * ) '  [16] total aerodynamic effective moment'
      WRITE( 239 , * ) '  [17] mechanical power extracted from the air flow'
      WRITE( 239 , * ) '  [18] moments ratio [07]/[16]'
      WRITE( 239 , * ) '  [19] wind speed'
      WRITE( 239 , 400 )
      WRITE( 239 , 550 , advance='no'  ) ' [01]           [02]           [03]           [04]           [05]           [06]           [07]           [08]           [09]    [10]    [11]     [12]     [13]     [14]     [15]           [16]           [17]           [18]         [19]'
      WRITE( 239 , 550 , advance='no'  ) ' ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
      
    ELSEIF( TRIM(actio)=='RESTAR' )THEN ! just open existing file and append data to end
      
      ! FORCES AND MOMENTS --------------------------------------------------
      fileName = outputRoot(1:outputRootLen)//'.FM.dat'
      OPEN( UNIT = 233 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'append' )
      WRITE( 233 , * ) '## FILE RESTARTED FROM HERE'
      
      ! AERODYNAMIC TIME ----------------------------------------------------
      fileName = outputRoot(1:outputRootLen)//'.AT.dat'
      OPEN( UNIT = 234 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'append' )
      WRITE( 234 , * ) '## FILE RESTARTED FROM HERE'
      
      ! AEROELASTIC REGISTRY ------------------------------------------------
      fileName = outputRoot(1:outputRootLen)//'.AR.dat'
      OPEN( UNIT = 238 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'append' )
      WRITE( 238 , * ) '## FILE RESTARTED FROM HERE'
      
      ! GENERATOR REGISTRY --------------------------------------------------
      fileName = outputRoot(1:outputRootLen)//'.GE.dat'
      OPEN( UNIT = 239 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'append' )
      WRITE( 239 , * ) '## FILE RESTARTED FROM HERE'
      
    ENDIF
    
    
    400 FORMAT(/)
    510 FORMAT(A,A,/)
    525 FORMAT(A,A)
    550 FORMAT(A,/)
    
  ENDSUBROUTINE inter_tp1
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_tp3
    
    ! INTERaction_TecPlot3
    
    ! Mauro S. Maza - 18/12/2012
    
    ! Closes output files
    
    IMPLICIT NONE
    
    CLOSE( UNIT = 232 ) ! root_Aero.@i.tec
    CLOSE( UNIT = 235 ) ! root_PF.@i.tec
    CLOSE( UNIT = 237 ) ! root_Stru.@i.tec
    
  ENDSUBROUTINE inter_tp3
    
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_tp5(ext)
    
    ! INTERaction_TecPlot5
  
    ! Mauro S. Maza - 18/12/2012
    
    ! Prints Tecplot's .lay files
    
    USE ctrl_db,        ONLY:   nstra  ! global control parameters
    USE inter_db,       ONLY:   n
    USE StructuresA,    ONLY:   Deltat
    
    IMPLICIT NONE
    
    CHARACTER(len=3), INTENT(IN) :: ext ! files' extension
    CHARACTER(len=mich):: inttoch ! FUNCTION: Puts an integer into a string
    CHARACTER(len=mlin):: fileName
    INTEGER :: i, fieldMapNum
    
    
    ! OPEN FILE
    fileName = outputRoot(1:outputRootLen)//'.TP.'//TRIM(ext)//TRIM('.lay')
    OPEN( UNIT = 236 , FILE = TRIM(fileName) , STATUS = 'unknown' , FORM = 'formatted' , ACCESS = 'sequential' )
    
    ! BODY
    ! Header
    WRITE( 236 , '(a)' ) '#!MC 1100'
    ! Files involved
    WRITE( 236 , 100 , advance='no'  ) '$!VarSet |LFDSFN1| = '''
    IF( aerGriOut )THEN
      DO i=1,nstra
        fileName = outputRoot(1:outputRootLen)//'.@'//TRIM(inttoch(i,0))//'.Aero.'//TRIM(ext)
        WRITE( 236 , 200 , advance='no'  ) '"', TRIM(fileName), '" '
      ENDDO
    ENDIF
    IF( strMesOut )THEN
      DO i=1,nstra
        fileName = outputRoot(1:outputRootLen)//'.@'//TRIM(inttoch(i,0))//'.Stru.'//TRIM(ext)
        WRITE( 236 , 200 , advance='no'  ) '"', TRIM(fileName), '" '
      ENDDO
    ENDIF
    IF( forcesOut )THEN
      DO i=1,nstra
        fileName = outputRoot(1:outputRootLen)//'.@'//TRIM(inttoch(i,0))//'.PFcs.'//TRIM(ext)
        WRITE( 236 , 200 , advance='no'  ) '"', TRIM(fileName), '" '
      ENDDO
    ENDIF
    ! some info ¿?
    WRITE( 236 , 100 ) ''''
    WRITE( 236 , '(a)' ) '$!VarSet |LFDSVL1| = ''"X" "Y" "Z" "n(X)" "n(Y)" "n(Z)" "Cp"'''
    WRITE( 236 , '(a)' ) '$!SETSTYLEBASE FACTORY'
    WRITE( 236 , '(a)' ) '### Frame Number 1 ###'
    WRITE( 236 , '(a)' ) '$!READDATASET  ''|LFDSFN1|'' '
    WRITE( 236 , '(a)' ) '  INITIALPLOTTYPE = CARTESIAN3D'
    WRITE( 236 , '(a)' ) '  INCLUDETEXT = NO'
    WRITE( 236 , '(a)' ) '  INCLUDEGEOM = NO'
    WRITE( 236 , '(a)' ) '  ASSIGNSTRANDIDS = YES'
    WRITE( 236 , '(a)' ) '  VARLOADMODE = BYNAME'
    WRITE( 236 , '(a)' ) '  VARNAMELIST = ''|LFDSVL1|'' '
    WRITE( 236 , '(a)' ) '$!REMOVEVAR |LFDSVL1|'
    WRITE( 236 , '(a)' ) '$!REMOVEVAR |LFDSFN1|'
    WRITE( 236 , '(a)' ) '$!FRAMELAYOUT '
    WRITE( 236 , '(a)' ) '  SHOWHEADER = NO'
    WRITE( 236 , '(a)' ) '$!PLOTTYPE  = CARTESIAN3D'
    WRITE( 236 , '(a,i,a)' ) '$!ACTIVEFIELDMAPS  =  [1-', NINT(zonesPerIter), ']'
    WRITE( 236 , '(a)' ) '$!GLOBALTIME '
    WRITE( 236 , '(a,e)' ) '  SOLUTIONTIME = ', (n-1)*Deltat
    WRITE( 236 , '(a)' ) '$!GLOBALTHREEDVECTOR '
    WRITE( 236 , '(a)' ) '  UVAR = 4'
    WRITE( 236 , '(a)' ) '  VVAR = 5'
    WRITE( 236 , '(a)' ) '  WVAR = 6'
    WRITE( 236 , '(a)' ) '  RELATIVELENGTH = 0.005'   
    ! zones' styles
    fieldMapNum = 0
    IF( aerGriOut )THEN
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLUE}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLUE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLUE}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLUE}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLUE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLUE}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLUE}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLUE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLUE}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLUE}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLUE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLUE}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLUE}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLUE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLUE}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLUE}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLUE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLUE}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = RED}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = RED USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = RED}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = YELLOW}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = YELLOW USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = YELLOW}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = GREEN}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = GREEN USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = GREEN}'
    ENDIF
    IF( strMesOut )THEN
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLACK}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLACK USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLACK}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLACK}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLACK USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLACK}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLACK}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLACK USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLACK}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLACK}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLACK USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLACK}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLACK}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLACK USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLACK}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {COLOR = BLACK}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = BLACK USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  SHADE {COLOR = BLACK}'
    ENDIF
    IF( forcesOut )THEN
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = PURPLE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {COLOR = PURPLE LINETHICKNESS = 0.2}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = PURPLE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {COLOR = PURPLE LINETHICKNESS = 0.2}'
      fieldMapNum = fieldMapNum + 1
      WRITE( 236 , '(a,i,a)' ) '$!FIELDMAP  [', fieldMapNum, ']'
      WRITE( 236 , '(a)' ) '  MESH {SHOW = NO}'
      WRITE( 236 , '(a)' ) '  CONTOUR {CONTOURTYPE = FLOOD COLOR = PURPLE USELIGHTINGEFFECT = YES}'
      WRITE( 236 , '(a)' ) '  VECTOR {COLOR = PURPLE LINETHICKNESS = 0.2}'
    ENDIF
    WRITE( 236 , '(a)' ) '$!VIEW FIT'
    WRITE( 236 , '(a)' ) '$!THREEDAXIS '
    WRITE( 236 , '(a)' ) '  AXISMODE = XYZDEPENDENT'
    WRITE( 236 , '(a)' ) '  XYDEPXTOYRATIO = 1'
    WRITE( 236 , '(a)' ) '  DEPXTOYRATIO = 1'
    WRITE( 236 , '(a)' ) '  DEPXTOZRATIO = 1'
    WRITE( 236 , '(a)' ) '$!THREEDAXIS '
    WRITE( 236 , '(a)' ) '  XDETAIL {SHOWAXIS = NO}'
    WRITE( 236 , '(a)' ) '$!THREEDAXIS '
    WRITE( 236 , '(a)' ) '  YDETAIL {SHOWAXIS = NO}'
    WRITE( 236 , '(a)' ) '$!THREEDAXIS '
    WRITE( 236 , '(a)' ) '  ZDETAIL {SHOWAXIS = NO}'
    WRITE( 236 , '(a)' ) '$!GLOBALTHREED '
    WRITE( 236 , '(a)' ) '  AXISSCALEFACT {X = 1 Y = 1 Z = 1}'
    WRITE( 236 , '(a)' ) '$!FIELDLAYERS '
    WRITE( 236 , '(a)' ) '  SHOWVECTOR = YES'
    WRITE( 236 , '(a)' ) '  SHOWSHADE = YES'
    WRITE( 236 , '(a)' ) '  SHOWEDGE = NO'
    WRITE( 236 , '(a)' ) '$!STREAMTRACELAYERS '
    WRITE( 236 , '(a)' ) '  SHOW = NO'
    WRITE( 236 , '(a)' ) '$!SETSTYLEBASE CONFIG'
    
    ! CLOSE FILE
    CLOSE( UNIT = 236 )
    
    100 FORMAT(A)
    200 FORMAT(3A)
  
  ENDSUBROUTINE inter_tp5
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_fam
    
    ! INTERaction_ForcesAndMoments
  
    ! Mauro S. Maza - 14/03/2013
    !                 29/09/2015
    
    ! Calculates total forces applied over each blade in aero. grid.
    ! Calculates total effective moment (over rotor) of those forces.
    
    USE StructuresA,    ONLY: nsec, npan!, Deltat
    USE ctrl_db,        ONLY: ttime
    USE inter_db,       ONLY: hub_master, CF2F!,                &
    !                          rotorEqu, rotorDam,        &
    !                          theta, omOld, omNew, omDot
    !USE StructuresA,    ONLY: q1, q2, q3
    USE npo_db,         ONLY: coora, euler, velnp
    !USE control_db,     ONLY: ctrlData, ctrlAct, brkMom, ltar, cptar, qtar, wSpeed
    !USE control_sr,     ONLY: ssiWSpeed
    
    IMPLICIT NONE
    
    INTEGER                         :: i, j, m, fp, countBlades !, n1
    REAL                            :: Force(3), Faero(3),          &
                                       Momen(3), MoEff, MoEfT,      &
                                       xyz_a(3), xyz_r(3), r(3),    &
                                       t1(3) !, & !, viscPow
                                       !Mfrel(3), lambda1(3,3), Mfree(3), Mstru, Fstru(3)
    !REAL                            :: ltarW, cptarW, qtarW ! control variables for writting
    TYPE(pair), 		POINTER		:: posic
    
    ! FORCES AND MOMENTS ----------------------------------------------------
    ! First column
    WRITE( 233 , 500 , advance='no' ) ttime
    
    xyz_r = coora(:,hub_master) ! Hub Master location
    t1 = euler(1:3,hub_master)  ! Hub axis of rotation
    t1 = (1/SQRT(DOT_PRODUCT(t1,t1)))*t1 ! asure it is versor
    
    MoEfT = 0 ! one per iteration
    IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	  posic => headpair
      
      countBlades = 0
      DO ! Loop over pairs
        
        IF( TRIM(posic%m_aer)=='blade1' .OR. &      ! no forces for non-lifting surfaces
            TRIM(posic%m_aer)=='blade2' .OR. &
            TRIM(posic%m_aer)=='blade3'        )THEN
          
          countBlades = countBlades + 1
          
          Faero(:) = 0 ! one for each blade per iteration
          ! Aero. forces and moment
          DO i=1,nsec   ! loop over sections in lifting surface
            fp = posic%section(i)%panels(1)  ! section first panel
            DO j=1,npan ! loop over panels in present section
              m = fp - 1 + j
              
              Force = posic%panel(m)%CF * CF2F
              xyz_a = posic%panel(m)%xyzcp ! Control Point location
              r = xyz_a - xyz_r ! relative position
              CALL cross_product( Momen , r , Force ) ! Moment with respect to the Hub Master Node
              MoEff = DOT_PRODUCT(Momen,t1) ! effective moment
              
              Faero = Faero + Force
              MoEfT = MoEfT + MoEff ! Total Effective Moment
              
            ENDDO
          ENDDO
          
          ! -----------------------------------------------------------------
          ! -----------------------------------------------------------------
          !Fstru(:) = 0 ! one for each blade per iteration
          !! Struct. forces and moments
          !DO i=1,posic%numbnodes
          !  
          !  n1 = posic%nodes(i)
          !  Force = loa(1:3,n1)
          !  
          !  xyz_a = coora(:,n1)
          !  r = xyz_a - xyz_r ! relative position
          !  CALL cross_product( Momen , r , Force ) ! Moment with respect to the origin
          !  
          !  Mfrel = loa(4:6,n1) ! Free moment applied on struc. node n1, in local coordinates of node n1
          !  lambda1 = RESHAPE(euler(1:9,n1),(/3,3/))	! Rotation matrix for n1 local system
          !  Mfree = MATMUL(lambda1,Mfrel)   ! Free moment applied on struc. node n1, in global coordinates
          !  
          !  MoEff = DOT_PRODUCT(Momen+Mfree,t1) ! effective moment
          !  
          !  Fstru = Fstru + Force
          !  Mstru = Mstru + MoEff ! Total Effective Moment
          !  
          !ENDDO
          ! -----------------------------------------------------------------
          ! -----------------------------------------------------------------
          
          ! Print
          IF( countBlades==3 )THEN
            WRITE( 233 , 600 ) Faero(1), Faero(2), Faero(3)
          ELSE
            WRITE( 233 , 600 , advance='no' ) Faero(1), Faero(2), Faero(3)
          ENDIF
          
        ENDIF
        
        IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
          posic => posic%next
        ELSE
          exit
        ENDIF
        
      ENDDO
    ENDIF
    
    !! GENERATOR REGISTRY ----------------------------------------------------
    !omNew = velnp(4,hub_master) ! rotor angular velocity (new)
    !IF( n==1 )THEN
    !  theta = 0 ! rotor angular position
    !  omDot = 0 ! rotor angular acceleration
    !ELSE
    !  theta = theta + omNew*Deltat ! rotor angular position
    !  omDot = (omNew-omOld)/Deltat ! rotor angular acceleration
    !ENDIF
    !omOld = omNew
    !
    !IF( (.NOT.ctrlData) .OR. (.NOT.ctrlAct) )THEN ! no control system OR not active
    !  ltarW = 0.0
    !  cptarW = 0.0
    !  qtarW = 0.0
    !ELSE
    !  ltarW = ltar
    !  cptarW = cptar
    !  qtarW = qtar
    !ENDIF
    !
    !CALL ssiWSpeed
    !
    !! Data in columns: lambda, cptar, qtar y ponerlos a cero si no hay control
    !!   [01] iteration number
    !!   [02] time
    !!   [03] rotor angular position
    !!   [04] rotor angular velocity
    !!   [05] rotor angular acceleration
    !!   [06] damping coefficient
    !!   [07] equivalent viscous moment
    !!   [08] power absorbed by de generator
    !!   [09] applied brake moment
    !!   [10] advance coefficient (lambda)
    !!   [11] target power coefficient (c_P-tar)
    !!   [12] target blades' pitch angle (q_tar)
    !!   [13] q1 - pith angle for balde 1
    !!   [14] q2 - pith angle for balde 2
    !!   [15] q3 - pith angle for balde 3
    !!   [16] total aerodynamic effective moment
    !!   [17] mechanical power extracted from the air flow
    !!   [18] moments ratio [07]/[16]
    !!   [19] wind speed
    !!   [20] wind direction
    !IF( rotorEqu .GT. 0 )THEN !       [01]   [02]   [03]   [04]   [05]      [06]            [07]               [08]    [09]   [10]    [11]           [12]        [13]        [14]        [15]   [16]         [17]                  [18]   [19]
    !  WRITE( 239 , 700 , advance='no' ) n, ttime, theta, omNew, omDot, rotorDam, omNew*rotorDam, omNew**2*rotorDam, brkMom, ltarW, cptarW, qtarW*rad2deg, q1*rad2deg, q2*rad2deg, q3*rad2deg, MoEfT, omNew*MoEfT, omNew*rotorDam/MoEfT, wSpeed
    !ELSE
    !  WRITE( 239 , 800 , advance='no' ) n, ttime, theta, omNew, omDot,  '    -',        '    -',           '    -', brkMom, ltarW, cptarW, qtarW*rad2deg, q1*rad2deg, q2*rad2deg, q3*rad2deg, MoEfT, omNew*MoEfT,              '    -', wSpeed
    !ENDIF
    
    500 FORMAT(E)
    600 FORMAT(E,E,E)
    700 FORMAT( I5, 4E15.6, 3E15.6, E15.6, F8.2, F8.3, 4F9.2, 2E15.6, E15.6, E13.4, /)
    800 FORMAT( I5, 4E15.6, 3A15  , E15.6, F8.2, F8.3, 4F9.2, 2E15.6, A15  , E13.4, /)
  
  ENDSUBROUTINE inter_fam
  
  SUBROUTINE inter_fam_debug
    


    ! WARNING!!!!!!!!!!!
    ! THIS SUBROUTINE WAS WRITTEN FOR DEBUGING PURPOSSES ONLY
    ! RESTORE IT FOR NORMAL EXCECUTION
    
    
    
    ! INTERaction_ForcesAndMoments
  
    ! Mauro S. Maza - 14/03/2013
    !                 29/09/2015
    
    ! Calculates total force and total moment w.r.t. hub node
    ! applied over Blade2 aero. grid. and structural mesh
    ! FOR DEBUGGING PURPOSES
    
    USE StructuresA,    ONLY: nsec, npan!, Deltat
    USE ctrl_db,        ONLY: ttime
    USE inter_db,       ONLY: hub_master, CF2F, loa !,                &
    !                          rotorEqu, rotorDam,        &
    !                          theta, omOld, omNew, omDot
    !USE StructuresA,    ONLY: q1, q2, q3
    USE npo_db,         ONLY: coora, euler, velnp
    !USE control_db,     ONLY: ctrlData, ctrlAct, brkMom, ltar, cptar, qtar, wSpeed
    !USE control_sr,     ONLY: ssiWSpeed
    
    IMPLICIT NONE
    
    INTEGER                         :: i, j, m, fp, countBlades, n1
    REAL                            :: Force(3), Faero(3),          &
                                       Momen(3), Maero(3), MoEfT,      &
                                       xyz_a(3), xyz_r(3), r(3),    &
                                       t1(3), & !, viscPow
                                       Mfrel(3), lambda1(3,3), Mfree(3), Mstru(3), Fstru(3)
    !REAL                            :: ltarW, cptarW, qtarW ! control variables for writting
    TYPE(pair), 		POINTER		:: posic
    
    ! FORCES AND MOMENTS ----------------------------------------------------
    ! First column
    WRITE( 233 , 500 , advance='no' ) ttime
    
    xyz_r = coora(:,hub_master) ! Hub Master location
    t1 = euler(1:3,hub_master)  ! Hub axis of rotation
    t1 = (1/SQRT(DOT_PRODUCT(t1,t1)))*t1 ! asure it is versor
    
    IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
	  posic => headpair
      
      countBlades = 0
      DO ! Loop over pairs
        
        IF( TRIM(posic%m_aer)=='blade2' )THEN
          
          Faero(:) = 0 ! one per iteration
          Maero(:) = 0 ! one per iteration
          ! Aero. forces and moments
          DO i=1,nsec   ! loop over sections in lifting surface
            fp = posic%section(i)%panels(1)  ! section first panel
            DO j=1,npan ! loop over panels in present section
              m = fp - 1 + j
              
              Force = posic%panel(m)%CF * CF2F
              xyz_a = posic%panel(m)%xyzcp ! Control Point location
              r = xyz_a - xyz_r ! relative position
              CALL cross_product( Momen , r , Force ) ! Moment with respect to the Hub Master Node
              
              Faero = Faero + Force
              Maero = Maero + Momen
              
            ENDDO
          ENDDO
          
          ! -----------------------------------------------------------------
          ! -----------------------------------------------------------------
          Fstru(:) = 0 ! one per iteration
          Mstru(:) = 0 ! one per iteration
          ! Struct. forces and moments
          DO i=1,posic%numbnodes
           
            n1 = posic%nodes(i)
            Force = loa(1:3,n1)
            
            xyz_a = coora(:,n1)
            r = xyz_a - xyz_r ! relative position
            CALL cross_product( Momen , r , Force ) ! Moment with respect to the origin
            
            Mfrel = loa(4:6,n1) ! Free moment applied on struc. node n1, in local coordinates of node n1
            lambda1 = RESHAPE(euler(1:9,n1),(/3,3/))	! Rotation matrix for n1 local system
            Mfree = MATMUL(lambda1,Mfrel)   ! Free moment applied on struc. node n1, in global coordinates
            
            Fstru = Fstru + Force
            Mstru = Mstru + Momen+Mfree ! Total Effective Moment
           
          ENDDO
          ! -----------------------------------------------------------------
          ! -----------------------------------------------------------------
          
          ! Print
          WRITE( 233 , 600 , advance='no' ) Faero(1), Faero(2), Faero(3), Maero(1), Maero(2), Maero(3), Fstru(1), Fstru(2), Fstru(3), Mstru(1), Mstru(2), Mstru(3)
          
        ENDIF
        
        IF (ASSOCIATED (posic%next)) THEN	! Check if it's not the end of the list
          posic => posic%next
        ELSE
          exit
        ENDIF
        
      ENDDO
    ENDIF
    
    500 FORMAT(E)
    600 FORMAT(12E,/)
    700 FORMAT( I5, 4E15.6, 3E15.6, E15.6, F8.2, F8.3, 4F9.2, 2E15.6, E15.6, E13.4, /)
    800 FORMAT( I5, 4E15.6, 3A15  , E15.6, F8.2, F8.3, 4F9.2, 2E15.6, A15  , E13.4, /)
  
  ENDSUBROUTINE inter_fam_debug
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE inter_reg(stage)
    
    ! INTERaction_REGistry
    
    ! Mauro S. Maza - 18/01/2013
    !                 02/09/2015
    
    ! Prints read and calculated variables
    ! "Aeroelastic Registry" in root.AR.dat
    
    USE StructuresA,    ONLY:   NS, NSA, cutoff, &
                                lref, vref, Deltat, DtVel, DtFact, vfs, Alphadeg, expo, HVfs, &
                                nn, np, nseg, nsec, npan, nnwke, npu1, npu2, &
                                nnh, nnn, nnt, nph, npn, npt, nsegh, nsegn, nsegt, &
                                phi, beta, q1, q2, q3, q4, q5, &
                                Areat
    USE Constants               ! pi, deg2rad, rad2deg
    USE inter_db,       ONLY:   numnest
    USE genr_db,        ONLY:   rotorDam
    USE wind_db,        ONLY:   DtMax, rVel, rOm, nSmoothOm
    USE ctrl_db,        ONLY:   ttime, dtime
    
    IMPLICIT NONE
    
    CHARACTER(len=*), INTENT(IN)    :: stage
    
    IF( TRIM(stage) == 'init' )THEN
      WRITE( 238 , 300 ) ' INPUT DATA READ ============================================'
      WRITE( 238 , 300 ) ' SIMULATION -------------------------------------------------'
      WRITE( 238 , 500 ) ' Number of time steps (NS):                    ', NS
      WRITE( 238 , 500 ) ' Number of time steps for wakes (NSA):         ', NSA
      WRITE( 238 , 600 ) ' Cutoff (cutoff):                              ', cutoff
      WRITE( 238 , 600 ) ' Velocity for calculating time step (DtVel):   ', DtVel
      WRITE( 238 , 600 ) ' Factor for changing time step (DtFact):       ', DtFact
      WRITE( 238 , 600 ) ' Max allowable Aero Time Step (DtMax):         ', DtMax
      WRITE( 238 , 600 ) ' Free Stream Velocity Norm "with sign" (vfs):  ', vfs
      WRITE( 238 , 600 ) ' Free Stream Velocity Angle (Alphadeg):        ', Alphadeg
      WRITE( 238 , 600 ) ' Exponent for earth boundary layer (expo):     ', expo
      WRITE( 238 , 600 ) ' Reference velocity height (HVfs):             ', HVfs
      WRITE( 238 , 450 ) ' Use rotor speed for Deltat (rVel)?:           ', rVel
      WRITE( 238 , 600 ) ' Reference length if rVel=.TRUE. (rOm):        ', rOm
      WRITE( 238 , 500 ) ' Smoothening factor if rVel=.TRUE. (rOm):      ', nSmoothOm
      WRITE( 238 , 400 )
      WRITE( 238 , 300 ) ' BLADES -----------------------------------------------------'
      WRITE( 238 , 500 ) ' Number of nodes (nn):                         ', nn
      WRITE( 238 , 500 ) ' Number of panels (np):                        ', np
      WRITE( 238 , 500 ) ' Number of segments (nseg):                    ', nseg
      WRITE( 238 , 500 ) ' Number of sections (nsec):                    ', nsec
      WRITE( 238 , 500 ) ' Number of panels per section (npan):          ', npan
      WRITE( 238 , 500 ) ' Number of wake nodes (nnwke):                 ', nnwke
      WRITE( 238 , 500 ) ' Number of union panel in upper side (npu1):   ', npu1
      WRITE( 238 , 500 ) ' Number of union panel in lower side (npu2):   ', npu2
      WRITE( 238 , 400 )
      WRITE( 238 , 300 ) ' HUB --------------------------------------------------------'
      WRITE( 238 , 500 ) ' Number of nodes (nnh):                        ', nnh
      WRITE( 238 , 500 ) ' Number of panels (nph):                       ', nph
      WRITE( 238 , 500 ) ' Number of segments (nsegh):                   ', nsegh
      WRITE( 238 , 400 )
      WRITE( 238 , 300 ) ' NACELLE ----------------------------------------------------'
      WRITE( 238 , 500 ) ' Number of nodes (nnn):                        ', nnn
      WRITE( 238 , 500 ) ' Number of panels (npn):                       ', npn
      WRITE( 238 , 500 ) ' Number of segments (nsegn):                   ', nsegn
      WRITE( 238 , 400 )
      WRITE( 238 , 300 ) ' TOWER ------------------------------------------------------'
      WRITE( 238 , 500 ) ' Number of nodes (nnt):                        ', nnt
      WRITE( 238 , 500 ) ' Number of panels (npt):                       ', npt
      WRITE( 238 , 500 ) ' Number of segments (nsegt):                   ', nsegt
      WRITE( 238 , 400 )
      WRITE( 238 , 300 ) ' ENTIRE TURBINE GEOMETRY (some) -----------------------------'
      WRITE( 238 , 700 ) ' Pre-Tilt angle (phi):                         ', phi*rad2deg
      WRITE( 238 , 700 ) ' Pre-Cone angle (beta):                        ', beta*rad2deg
      WRITE( 238 , 700 ) ' Blade 1 initial pitch angle (q1):             ', q1*rad2deg
      WRITE( 238 , 700 ) ' Blade 2 initial pitch angle (q2):             ', q2*rad2deg
      WRITE( 238 , 700 ) ' Blade 3 initial pitch angle (q3):             ', q3*rad2deg
      WRITE( 238 , 700 ) ' Rotor initial rotation angle (q4):            ', q4*rad2deg
      WRITE( 238 , 700 ) ' Nacelle-Rotor initial yaw angle (q5):         ', q5*rad2deg
      WRITE( 238 , 400 )
      WRITE( 238 , 300 ) ' TECPLOT OUTPUT ---------------------------------------------'
      IF( aerGriOut )THEN
        WRITE( 238 , 300 ) ' Aerodynamic Grid:           YES'
      ELSE
        WRITE( 238 , 300 ) ' Aerodynamic Grid:            NO'
      ENDIF
      IF( strMesOut )THEN
        WRITE( 238 , 300 ) ' Structural Mesh:            YES'
      ELSE
        WRITE( 238 , 300 ) ' Structural Mesh:             NO'
      ENDIF
      IF( forcesOut )THEN
        WRITE( 238 , 300 ) ' Blades'' Panels'' Forces:     YES'
      ELSE
        WRITE( 238 , 300 ) ' Blades'' Panels'' Forces:      NO'
      ENDIF
      WRITE( 238 , 400 )
      WRITE( 238 , 300 ) ' GENERATOR --------------------------------------------------'
      WRITE( 238 , 600 ) ' Rotor damp coefficient (rotorDam)             ', rotorDam
      WRITE( 238 , 400 )
      WRITE( 238 , 300 ) ' CALCUALTED VARIABLES ======================================='
      WRITE( 238 , 600 ) ' Total area of one lifting surface (Areat):    ', Areat
      WRITE( 238 , 600 ) ' Reference length (lref):                      ', lref
      WRITE( 238 , 600 ) ' Reference velocity (vref):                    ', vref
      WRITE( 238 , 600 ) ' Initial Aero Time Step (Deltat):              ', Deltat
      WRITE( 238 , 300 ) ' Structural and Aerodynamic Time Steps, their Ratio and Max duration'
      WRITE( 238 , 300 ) ' of strategy (not to exceed max number of zones Tecplot can manage):'
      WRITE( 238 , 300 ) '  t             STS           ATS           ATS/STS       tpTS'
      WRITE( 238 , 300 ) ' ---------------------------------------------------------------------'
    ELSEIF( TRIM(stage) == 'dtime' )THEN
      ! IF( ( ABS( ( dtimeOld- dtime)/ dtimeOld ) .GE. 0.05 )  .OR. &
          ! ( ABS( (DeltatOld-Deltat)/DeltatOld ) .GE. 0.05 )  )THEN
      IF( ( dtimeOld  .NE. dtime )  .OR. &
          ( DeltatOld .NE. Deltat)  )THEN
        tMaxTecPlot = 32000 / zonesPerIter * Deltat
        WRITE( 238 , 650 ) ttime, dtime, Deltat, REAL(numnest), tMaxTecPlot
        dtimeOld = dtime
        DeltatOld = Deltat
      ENDIF
    ENDIF
    
    300 FORMAT(A)       ! text
    400 FORMAT(/)       ! newline
    450 FORMAT(A,L14)   ! text and logical
    500 FORMAT(A,I14)   ! text and integer
    600 FORMAT(A,E14.6) ! text and general numbers
    650 FORMAT(5E14.6)  ! general numbers in 5-column table
    700 FORMAT(A,F14.2) ! text and angles
  
  ENDSUBROUTINE inter_reg
  
ENDMODULE introut_sr

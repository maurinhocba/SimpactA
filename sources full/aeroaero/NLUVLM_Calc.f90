SUBROUTINE NLUVLM_Calc
  
  ! Mauro S. Maza - 13/10/2011
  ! From Cristian Gebhardt's NLUVLM_Prog
  
  ! Non-Linear Unsteady Vortex Latice Method
  ! Calculates, among others, aerodynamical loads
  ! on blades' lifting surfaces.
  
  USE StructuresA,  ONLY: Alpha, Alphadeg, Vinf, vfs,                   &
                          np, nph, npn, npt,                            &
                          blade1panels,  blade2panels, blade3panels,    &
                             hubpanels, nacellepanels,  towerpanels
  USE ArraysA,      ONLY: RHS, RHS1, RHS2, A, IPIV, gamma
  USE ctrl_db,      ONLY: ndofn, npoin
  USE inter_db,     ONLY: n, loa
  USE introut_db,   ONLY: msg
  USE introut_sr,   ONLY: introut_mng
  USE wind_sr,      ONLY: DtCalc
  USE Constants
    
  INTEGER                         :: INFO
  
  
  IF( n==0 )THEN
    ALLOCATE( loa(ndofn,npoin) )
  ELSE  ! Old (previous) rings' vorticity
    blade1panels( : )%gammaold  = blade1panels( : )%gamma
    blade2panels( : )%gammaold  = blade2panels( : )%gamma
    blade3panels( : )%gammaold  = blade3panels( : )%gamma
    hubpanels( : )%gammaold     = hubpanels( : )%gamma
    nacellepanels( : )%gammaold = nacellepanels( : )%gamma
    towerpanels( : )%gammaold   = towerpanels( : )%gamma
  ENDIF
  
  ! Normal versors at control points
  CALL CP_normal
  
  ! Right Hand Side for non-penetration condition
  CALL RHS1_vector    ! Free stream and CP's translational velocities RHS
  IF( n/=0 )THEN
    CALL FixConvect( n )  ! Updates coordinates of wake nodes that are on the trailing edge
    CALL RHS2_vector( n ) ! Wake RHS
    RHS = RHS1 + RHS2
  ELSE
    RHS = RHS1
  ENDIF
  
  ! Print aero. grid and str. mesh
  CALL introut_mng('gridAndMesh')
  
  ! Coeficient matrix
  CALL A_matrix
  
  ! Solution - Current rings' vorticity
  CALL DGESV(3*np+nph+npn+npt,1,A,3*np+nph+npn+npt,IPIV,RHS,3*np+nph+npn+npt,INFO)
  IF(INFO/=0)THEN
    WRITE(msg,'(A,I,A)') 'there was a problem solving A*G=RHS - DGESV output variable INFO=',INFO,' - matrix A may be bad conditioned'
    CALL introut_mng('repProbl')
    STOP
  ENDIF
  gamma = RHS
  blade1panels( : )%gamma  = gamma( 1 : np )
  blade2panels( : )%gamma  = gamma( np + 1 : 2 * np )
  blade3panels( : )%gamma  = gamma( 2 * np + 1 : 3 * np )
  hubpanels( : )%gamma     = gamma( 3 * np + 1 : 3 * np + nph )
  nacellepanels( : )%gamma = gamma( 3 * np + nph + 1 : 3 * np + nph + npn )
  towerpanels( : )%gamma   = gamma( 3 * np + nph + npn + 1 : 3 * np + nph + npn + npt )
     
  ! Segments' voticity
  CALL Segment_properties( gamma )
  
  ! Loads
  CALL Loads( n )
  
  ! Update aero. time step
  CALL DtCalc(.TRUE.)
  
  ! Wake convection
  CALL Convect( n )
  
  ! Next step counter
  n = n + 1
  
  
  !CALL ringsGsText ! for debugging purposes only - file is never deleted by the program; data is simply added to the end - you may delet it manually
  
CONTAINS
    
  ! ------------------------------- debug -------------------------------
  ! =====================================================================
  
  SUBROUTINE ringsGsText
    
    ! prints txt file with vorticity of vortex rings, G, in blade 2
    
    ! Mauro S. Maza - 29/07/2016
    
    USE ctrl_db,        ONLY:   ttime
    USE inter_db,       ONLY:   headpair, pair
    
    IMPLICIT NONE
    
    ! internal vars
    TYPE(pair),         POINTER     ::  posic
    LOGICAL                             ::  fExist
    INTEGER                             ::  i
    
    
    !IF( n==1 )THEN ! commenting this IF, a tridimensional matrix will be printed, wich allows to analyse evolution of A with n
      INQUIRE(FILE='ringsGsText.txt', exist=fExist)
      IF(fExist)THEN
        OPEN(UNIT = 249, FILE='ringsGsText.txt', STATUS='old', POSITION='append', ACTION='write')
      ELSE
        OPEN(UNIT = 249, FILE='ringsGsText.txt', STATUS='new', ACTION='write')
        ! Heading
        WRITE( 249 , 400 )
        WRITE( 249 , * ) 'vorticity of vortex rings, G'
        WRITE( 249 , * ) '****************************'
        WRITE( 249 , * ) 'ONLY BLADE 2 NOW'
        WRITE( 249 , 400 )
        WRITE( 249 , * ) 'for debugging purposes only'
        WRITE( 249 , * ) 'file is never deleted by the program; data is simply added to the end'
        WRITE( 249 , * ) 'you may delet it manually'
        WRITE( 249 , 400 )
      ENDIF
      
      ! step title
      WRITE( 249 , * ) 'aeroStep            time'
      WRITE( 249 , * ) n, ttime
      IF( ASSOCIATED(headpair) )THEN		! Check if a list is empty
        posic => headpair
        DO ! Loop over pairs
          
          IF( TRIM(posic%m_aer)=='blade2' )THEN
            DO i=1,posic%np
              WRITE(249,600) blade2panels(i)%gamma
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
      
      CLOSE(UNIT = 249)
    !ENDIF
    
    400 FORMAT(/) ! new line
    600 FORMAT(1X,E12.5E2) ! 1 time vorticity
    
  ENDSUBROUTINE ringsGsText
      
END SUBROUTINE NLUVLM_Calc
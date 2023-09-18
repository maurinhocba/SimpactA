
MODULE struout_db
  
  IMPLICIT NONE
  
  LOGICAL                       :: strMesOut ! generate output files with Structural Mesh
  
CONTAINS ! =============================================
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_struout
    
    ! DUMPINg_STRUctural OUTput data
    
    ! Mauro S. Maza - 18/01/2013
    
    ! Dumps output data in restart file
    
    IMPLICIT NONE
    
    WRITE(50,ERR=9999) strMesOut
                       !outputRootLen, outputRoot: determined in restart too
    
    RETURN
    
    9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_struout
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_struout
    
    ! RESTARt_STRUctural OUTput data
    
    ! Mauro S. Maza - 18/01/2013
    
    ! Reads output data from restart file
    
    IMPLICIT NONE
    
    READ (51)          strMesOut
    
  ENDSUBROUTINE restar_struout
  
ENDMODULE struout_db

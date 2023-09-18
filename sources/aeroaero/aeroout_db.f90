
MODULE aeroout_db
  
  !USE param_db,         ONLY:   mlen ! Global parameters
  
  IMPLICIT NONE
  
  LOGICAL                       :: aerGriOut,           & ! generate output files with Aerodynamic Grid
                                   forcesOut              ! generate output files with forces in panels
  
CONTAINS ! =============================================
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_aeroout
    
    ! DUMPINg_AEROdynamic OUTput data
    
    ! Mauro S. Maza - 20/03/2015
    
    ! Dumps output data in restart file
    
    IMPLICIT NONE
    
    WRITE(50,ERR=9999) aerGriOut, forcesOut
    
    RETURN
    
    9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_aeroout
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_aeroout
    
    ! RESTARt_AEROdynamic OUTput data
    
    ! Mauro S. Maza - 20/03/2015
    
    ! Reads output data from restart file
    
    IMPLICIT NONE
    
    READ (51)          aerGriOut, forcesOut
    
  ENDSUBROUTINE restar_aeroout
  
ENDMODULE aeroout_db

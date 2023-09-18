
MODULE mastrnods_db
  
  ! Mauro S. Maza - 15/09/2015
  
  IMPLICIT NONE
  
CONTAINS ! =============================================
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_mastrnods
    
    ! DUMPINg_MASTEr NODeS data
    
    ! Mauro S. Maza - 07/01/2016
    
    ! Dumps output data in restart file
    
    IMPLICIT NONE
    
    
    !WRITE(50,ERR=9999) 
    !
    !RETURN
    !
    !9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_mastrnods
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_mastrnods
    
    ! RESTARt_MASTEr NODeS data
    
    ! Mauro S. Maza - 07/01/2016
    
    ! Reads output data from restart file
    
    IMPLICIT NONE
    
    
    READ (51)          
    
  ENDSUBROUTINE restar_mastrnods
  
ENDMODULE mastrnods_db

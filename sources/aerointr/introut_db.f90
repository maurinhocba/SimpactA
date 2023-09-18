
MODULE introut_db
  
  USE param_db,     ONLY:   mlen ! Global parameters
  
  IMPLICIT NONE
  
  LOGICAL                       :: newTPFileA,          & ! Aero Grid - .TRUE. if a new TecPlot file is to be started (connectivities should be printed again)
                                   newTPFileS             ! Structural Mesh
  
  REAL( kind = 8 )              :: tMaxTecPlot            ! Max duration of strategy not to exceed max number of zones Tecplot can manage
  REAL( kind = 8 )              :: zonesPerIter           ! Number of zones per iteration in TecPlot output
  REAL( kind = 8 )              :: dtimeOld  = 0.0        ! print dtime  in reg only if it changes
  REAL( kind = 8 )              :: DeltatOld = 0.0        ! print Deltat in reg only if it changes
  
  INTEGER                       :: outputRootLen
  CHARACTER(len=mlen)           :: outputRoot
  
  CHARACTER(len=200)            :: msg                    ! problem message to be reported
  
CONTAINS ! =============================================
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_introut
    
    ! DUMPINg_INTeRaction OUTput data
    
    ! Mauro S. Maza - 20/03/2015
    
    ! Dumps output data in restart file
    
    IMPLICIT NONE
    
    
    WRITE(50,ERR=9999) newTPFileA, newTPFileS, &
                       tMaxTecPlot, zonesPerIter, dtimeOld, DeltatOld
                       !outputRootLen, outputRoot: determined in restart too
                       !msg: not nec
    
    RETURN
    
    9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_introut
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_introut
    
    ! RESTARt_INTeRaction OUTput data
    
    ! Mauro S. Maza - 20/03/2015
    
    ! Reads output data from restart file
    
    IMPLICIT NONE
    
    
    READ (51)          newTPFileA, newTPFileS, &
                       tMaxTecPlot, zonesPerIter, dtimeOld, DeltatOld
    
  ENDSUBROUTINE restar_introut
  
ENDMODULE introut_db

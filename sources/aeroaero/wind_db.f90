MODULE wind_db
  
  ! Mauro S. Maza - 15/04/2014
  
  USE nrtype,       ONLY:   SP
  
  IMPLICIT NONE
  
  LOGICAL                           :: windFlag                 ! true if wind velocity is defined as a function of time
  LOGICAL                           :: rVel                     ! true if Deltat must be calculated with relative velocity
  INTEGER                           :: pointNum,              & ! number of data points in table defining wind velocity as a function of time
                                       nSmoothOm
  REAL(SP), ALLOCATABLE             :: windTab(:,:)             ! Wind module and angle table (for time dependent wind velocity)
  REAL( kind = 8 )                  :: DtMax,                 & ! Limiting value for Deltat (in case of small relative velocity)
                                       rOm                      ! Radius for translation velocity due to rotation of blades

CONTAINS !===================================================================
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_win
    
    ! DUMPINg_WINd data
    
    ! Mauro S. Maza - 21/04/2014
    
    ! Dumps wind data in restart file
    
    IMPLICIT NONE
    
    INTEGER::   i, j
    
    
    WRITE(50,ERR=9999) windFlag, rVel
    
    IF(windFlag)THEN
      WRITE(50,ERR=9999) pointNum, nSmoothOm
      DO i=1,pointNum
        WRITE(50,ERR=9999) (windTab(i,:j),j=1,6)
      ENDDO
    ENDIF
    
    WRITE(50,ERR=9999) DtMax, rOm
    
    
    RETURN
    
    9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_win
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_win
    
    ! RESTARt_WINd data
    
    ! Mauro S. Maza - 21/04/2014
    
    ! Reads wind data from restart file
    
    IMPLICIT NONE
    
    INTEGER::   i,j
    
    
    READ (51)          windFlag, rVel
    
    IF(windFlag)THEN
      READ (51)          pointNum, nSmoothOm
      ALLOCATE( windTab(pointNum,6) )
      DO i=1,pointNum
        READ (51)          (windTab(i,:j),j=1,6)
      ENDDO
    ENDIF
    
    READ (51)          DtMax, rOm
    
    
  ENDSUBROUTINE restar_win
  
ENDMODULE wind_db
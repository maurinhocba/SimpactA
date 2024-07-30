
MODULE wind_sr
  
  USE wind_db
  
  IMPLICIT NONE
  
  ! VARIABLES: INPUT, DETERMINATION AND USE
  !
  ! > Vinf(3)/StructuresA/aero_mods.f90
  !   - (unperturbed) velocity (vector) field, not considering earth's
  !     boundary layer
  !   - determined in instantWV/wind_sr.f90 from vfs and Alphadeg, wich can
  !     be set from
  !     + simulation.dat (IF .NOT.windFlag), or
  !     + table in wind.dat (IF windFlag)
  !   - used in Convect/Convect.f90 and RHS1_vector/NLUVLM.f90
  ! > vref/StructuresA/aero_mods.f90
  !   - reference wind speed (absolute value)
  !   - determined in dtClac/wind_sr.f90 from
  !     + DtVel in extras.dat (IF .NOT.windFlag), or
  !     + vfs from table in wind.dat (IF windFlag)
  !     may be modified in any case with rotor angular speed (IF rVel)
  !   - used in DtCalc/wind_sr.f90 to determine Deltat/StructuresA/aero_mods.f90

CONTAINS !===================================================================
  
  SUBROUTINE rdwind
    
    ! Mauro S. Maza - 15/04/2014
    ! Reads table defining wind velocity norm and angle as a function of time
    ! and determines cubic splines for later interpolation
    
    ! Time values must be ascending, not necessarily starting at 0 nor having
    ! a maximum value equal or greater than that of the maximum simulated t.
    ! Each time value mus be distinct from the others.
    
    ! Angles are in [deg]
    
    USE GenUseNR,       ONLY:   spline
    
    IMPLICIT NONE
    
    INTEGER :: i, j
    REAL(4), PARAMETER :: SP_zero=0
    
    
    OPEN( UNIT = 218, FILE = 'wind.dat', STATUS = 'old', FORM = 'formatted', ACCESS = 'sequential', IOSTAT = j )
    
    IF( j .GT. 0 )THEN ! wind velocity is constant
      windFlag = .FALSE.
    ELSE ! wind velocity is function of time
      windFlag = .TRUE.
      
      READ( 218, '(/i)' ) pointNum
      READ( 218, '(2(/))' )
      ALLOCATE( windTab(pointNum,6) )
      DO i=1,pointNum
        READ( 218,  * ) windTab(i,1:4)
      ENDDO
      CLOSE( UNIT = 218 )
      
      ! splines - 0.0 are for natural boundary condition (zero second derivative)
      CALL spline(windTab(:,2), windTab(:,3), SP_zero, SP_zero, windTab(:,5)) ! vel. norm
      CALL spline(windTab(:,2), windTab(:,4), SP_zero, SP_zero, windTab(:,6)) ! vel. angle
    ENDIF
    
  ENDSUBROUTINE rdwind

  ! Ejemplo de archivo de entrada - este renglón no va
  ! Número de puntos para la definición
  ! 2
  ! Tabla:
  ! #punto   tiempo     norma       ángulo
  ! #        [t]        [long/t]    [deg]
  !  1       0.000      180.0       180.0
  !  2       100.0      180.0       180.0
  
  !--------------------------------------------------------------------------
  
  SUBROUTINE instantWV
    
    ! Mauro S. Maza - 16/04/2014
    ! Determines Instantaneous Wind Velocity (vfs, Alphadeg, Alpha and Vinf)
    
    USE ctrl_db,        ONLY:   ttime
    USE StructuresA,    ONLY:   vfs, Vinf, Alphadeg, Alpha
    USE Constants
    USE GenUseNR,       ONLY:   splintMM
    
    IMPLICIT NONE
    
    
    ! interpolation
    IF(windFlag)THEN
      IF( ttime.LE.windTab(1,2) )THEN
        vfs      = windTab(1,3)
        Alphadeg = windTab(1,4)
      ELSEIF( ttime.GE.windTab(pointNum,2) )THEN
        vfs      = windTab(pointNum,3)
        Alphadeg = windTab(pointNum,4)
      ELSE
        vfs      = splintMM(windTab(:,2), windTab(:,3), windTab(:,5), REAL(ttime,SP)) ! vel. norm
        Alphadeg = splintMM(windTab(:,2), windTab(:,4), windTab(:,6), REAL(ttime,SP)) ! vel. angle
      ENDIF
    ENDIF
    
    Alpha = Alphadeg * deg2rad
    Vinf( 1 ) = vfs * DCOS( Alpha )
    Vinf( 2 ) = vfs * DSIN( Alpha )
    Vinf( 3 ) = 0.0D+0
    
  ENDSUBROUTINE instantWV
  
  !--------------------------------------------------------------------------
  
  SUBROUTINE DtCalc(flag)
    
    ! Mauro S. Maza - 16/04/2014
    
    ! Determines:
    !   - reference velocity (vref),
    !   - aerodynamic time interval (Deltat).
    ! Updates:
    !   - numnest: number of str. steps for each aero step
    !   - dtime: so Deltat=dtime*numnest
    
    USE inter_db,       ONLY:   hub_master, CF2F, numnest
    USE introut_sr,     ONLY:   introut_mng
    USE StructuresA,    ONLY:   vfs, vref, lref, Deltat, DtFact, DtVel, dens
    USE npo_db,         ONLY:   velnp
    USE ctrl_db,		ONLY:   dtime
    
    IMPLICIT NONE
    
    LOGICAL             :: flag
    REAL(kind=8)        :: omeg
    REAL(kind=8),SAVE   :: omegOld=0
    
    
    CALL instantWV ! determine vfs and Vinf
    
    IF(windFlag)THEN ! wind velocity is function of time
      vref = DABS(vfs)
    ELSE ! wind velocity is constant
      vref = DABS(DtVel)
    ENDIF
    
    IF(rVel)THEN ! reference velocity is the one relative to the blades
      omeg = velnp(4,hub_master) ! rotor speed
      omeg = (omeg  + omegOld*nSmoothOm)/(1+nSmoothOm) ! smoothening
      vref = SQRT( (vref)**2.0 + (rOm*omeg)**2.0 )
      omegOld = omeg
    ENDIF
    
	Deltat = lref / vref * DtFact   ! aero time step
    IF(Deltat.GT.DtMax)THEN
      Deltat = DtMax
    ENDIF
    
    ! modify dtime to arrive exactly to next aerodyn. calculation
    IF(flag)THEN
      !         Calculation of critical time step
      CALL timing(17,1)
      CALL incdlt ( )   ! determine dtime, no longer calculated until a new aerodyn. step is done
      CALL timing(17,2)
      numnest = MAX(CEILING(deltat/dtime),1)  ! determines num of estruct. steps for each aero. step - assumes that numnest > 0
      dtime = deltat/numnest                  ! so dtime*numnest=Deltat
      CALL introut_mng('dtime')               ! print new numnest in aero. reg.
    ENDIF
    
  ENDSUBROUTINE DtCalc
  
  !--------------------------------------------------------------------------
  
  SUBROUTINE WSpeed( z , vcoeff )
    
    ! Mauro S. Maza - 24/09/2013
    ! Determines variation of wind velocity norm due to earth boundary layer
    ! as for NREL FAST's HubHeight Wind Data File
    
    USE StructuresA, ONLY: HVfs, expo
    
    IMPLICIT NONE
    
  
    REAL(kind=8),INTENT(in):: z
    REAL(kind=8),INTENT(out):: vcoeff
  
    vcoeff = (ABS( z / HVfs ))**expo
  
    RETURN
  
  ENDSUBROUTINE WSpeed

ENDMODULE wind_sr
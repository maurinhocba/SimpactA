
MODULE genr_db
  
  !USE param_db,     ONLY:   mnam
  !USE ele08_db,     ONLY:   ele08_set, ele08
  !USE npo_db,       ONLY:   coord, coora, eule0, euler, velnp
  !USE StructuresA,  ONLY:   tnode, tpanel, tsection,                        & ! types
  !                          nn, nnh, nnn, nnt, np, nph, npn, npt,           & ! variables
  !                          blade1nodes, blade2nodes,  blade3nodes,         &
  !                          hubnodes,    nacellenodes, towernodes,          &
  !                          blade1panels, blade2panels,  blade3panels,      &
  !                          hubpanels,    nacellepanels, towerpanels,       &
  !                          blade1sections, blade2sections, blade3sections, &
  !                          vref, lref, Deltat, vfs, DtVel, DtFact,         &
  !                          phi, beta, theta1, theta2, theta3,              &
  !                          q1, q2, q3, q4, q5
  !USE Constants             ! pi, deg2rad, rad2deg
  
  IMPLICIT NONE
  
  INTEGER   ::  rotorEqu    ! Equation associated with hub axis rotation
  REAL      ::  rotorDam, & ! Rotor damp coefficient c, so F_visc = c*v
                rotorMas, & ! Mass corresponding to rotorEqu equation
                rotorAlp, & ! Rotor damp coefficient alpha (to add to damp), so F_visc = (2*alpha*m)*v
                theta,    & ! rotor angular position
                omOld,    & ! rotor angular velocity (old)
                omNew,    & ! rotor angular velocity (new)
                omDot       ! rotor angular acceleration
  
CONTAINS ! =============================================
  
  ! ---------- Subroutines for restart files managing -----------
  ! =============================================================
  
  SUBROUTINE dumpin_gen
    
    ! DUMPINg_GENerator data
    
    ! Mauro S. Maza - 18/01/2013
    
    ! Dumps generator data in restart file
    
    IMPLICIT NONE
    
    WRITE(50,ERR=9999) rotorEqu, rotorDam, rotorMas, rotorAlp, &
                       theta, omOld, omNew, omDot
    
    RETURN
    
    9999 CALL runen2('')
    
  ENDSUBROUTINE dumpin_gen
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE restar_gen
    
    ! RESTARt_GENerator data
    
    ! Mauro S. Maza - 18/01/2013
    
    ! Reads generator data from restart file
    
    IMPLICIT NONE
    
    READ (51)          rotorEqu, rotorDam, rotorMas, rotorAlp, &
                       theta, omOld, omNew, omDot
    
  ENDSUBROUTINE restar_gen
  
  
ENDMODULE genr_db


MODULE genr_sr
  
  USE genr_db
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
  
  
CONTAINS ! =============================================
  
  SUBROUTINE gener_inp
    
    ! Mauro S. Maza - 13/03/2013
    
    ! Reads generator data file
    
    IMPLICIT NONE
    
    OPEN( UNIT = 217, FILE = 'generator.dat', STATUS = 'old', FORM = 'formatted', ACCESS = 'sequential' )
    
    READ( 217 , '(3(/))' )
	READ( 217 , * ) rotorDam ! damp coefficient c, so F_visc = c*v
    
    CLOSE( UNIT = 217 )
    
  ENDSUBROUTINE gener_inp
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE include_damp
    
    ! Mauro S. Maza - 13/03/2013
    
    ! Calculates damp coefficient "alpha" for the rotor and adds it
    
    ! rotorDam: damp coefficient c,     so F_visc =      c     *v
    ! rotorAlp: damp coefficient alpha, so F_visc = (2*alpha*m)*v, needed to
    !           be used in dampin (..\mainp\dampin.f90)
    
    USE npo_db,     ONLY:   ifpre, ymass
    USE damp_db,    ONLY:   damp
    USE control_db, ONLY:   ctrlData
    !USE inter_db,   ONLY:   hub_master
    
    IMPLICIT NONE
    
    
    ! Find associated equation
    rotorEqu = ifpre(4,hub_master)
    
    ! Calculates alpha coefficient and adds it up
    IF( (rotorEqu .GT. 0)                       .AND. &
        ( (rotorDam .NE. 0.0) .OR. (ctrlData) ) )THEN
      rotorMas = 1/ymass(rotorEqu) ! ymass has inverses of DoFs' inertias
      rotorAlp = rotorDam / (2*rotorMas)
      damp(rotorEqu) = damp(rotorEqu) + rotorAlp
    ENDIF
    
    
  ENDSUBROUTINE include_damp
  
! ---------------------------------------------------------------------------
  
  SUBROUTINE update_damp
    
    ! Mauro S. Maza - 14/03/2013
    
    ! ReCalculates damp coeff. "alpha" for the rotor in case inertia changes
    
    ! rotorDam: damp coefficient c,     so F_visc =      c     *v
    ! rotorAlp: damp coefficient alpha, so F_visc = (2*alpha*m)*v, needed to
    !           be used in dampin (..\mainp\dampin.f90)
    
    USE npo_db,     ONLY:   ymass
    USE damp_db,    ONLY:   damp
    !USE inter_db,   ONLY:   hub_master
    
    IMPLICIT NONE
    
    
    IF( rotorEqu .GT. 0 )THEN
      IF( (rotorMas .LT. 0.99/ymass(rotorEqu)) .OR. &
          (rotorMas .GT. 1.01/ymass(rotorEqu))      )THEN ! rotorMas not in 1% of 1/ymass(rotorEqu)
        ! extract old damp coeff.
        damp(rotorEqu) = damp(rotorEqu) - rotorAlp
        ! recalculate and update
        rotorMas = 1/ymass(rotorEqu) ! ymass has inverses of DoFs inertias
        rotorAlp = rotorDam / (2*rotorMas)
        damp(rotorEqu) = damp(rotorEqu) + rotorAlp
      ENDIF
    ENDIF
    
    
  ENDSUBROUTINE update_damp
  
  
ENDMODULE genr_sr

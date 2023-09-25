
        SUBROUTINE Input_extras

! Cristian G. Gebhardt, Mauro S. Maza
! 13/03/2015

        USE StructuresA,    ONLY:   NSA, NSAL, cutoff, &
                                    DtVel, DtFact, &
                                    RBT, RTT, LT, &
                                    phi, beta, theta1, theta2, theta3, &
                                    lnb, lba, &
                                    lac1, lac2, lac3, &
                                    lad1, lad2, lad3, &
                                    lae1, lae2, lae3, &
                                    q1, q2, q3, q4, q5, &
                                    HVfs, expo, dens
        USE inter_db,       ONLY:   ini_pairs, add_pair,    &   ! Subroutines
                                    pair,                   &   ! Types
                                    headpair, tailpair,     &   ! Variables
                                    bEType
        USE aeroout_db,     ONLY:   aerGriOut, forcesOut
        USE struout_db,     ONLY:   strMesOut
        USE wind_db,        ONLY:   DtMax, rVel, rOm, nSmoothOm
        
        IMPLICIT NONE
        
        INTEGER                     :: i, & ! Indice de conteo.
                                       aeroGrid, struMesh, forces!, cps
        !LOGICAL                     :: restart=.FALSE. ! input for Alloc
        TYPE(pair), 		POINTER :: newpair
  
  
! Lectura de datos extras de la simulación y otros

		OPEN( UNIT = 216 , FILE = 'extras.dat' , STATUS = 'old' , FORM = 'formatted' , ACCESS = 'sequential' )

        READ( 216 , '(12(/))' ) 
		
        READ( 216 , * ) NSA
        READ( 216 , '(4(/))' )

        READ( 216 ,  * ) NSAL
        READ( 216 , '(4(/))' )

        READ( 216 ,  * ) cutoff
        READ( 216 , '(4(/))' )

        READ( 216 ,  * ) DtVel, DtFact, DtMax
        READ( 216 , '(4(/))' )

        READ( 216 ,  * ) rVel, rOm, nSmoothOm
        READ( 216 , '(9(/))' )

		READ( 216 ,  * ) RBT, RTT, LT
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) phi, beta
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) theta1, theta2, theta3
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) lnb, lba
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) lac1, lac2, lac3
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) lad1, lad2, lad3
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) lae1, lae2, lae3
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) q1, q2, q3
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) q4, q5
        READ( 216 , '(9(/))' )
        
        READ( 216 ,  * ) HVfs, expo
        READ( 216 , '(4(/))' )

		READ( 216 ,  * ) dens
        READ( 216 , '(9(/))' )
		
		CALL ini_pairs (headpair, tailpair)
		DO i=1,6
		  ALLOCATE(newpair)
		  CALL add_pair (newpair, headpair, tailpair)
		  READ( 216 ,  * ) newpair%m_aer, newpair%m_est
		ENDDO
        READ( 216 , '(/)' )
        READ( 216 ,  * ) bEType ! blade element type
        READ( 216 , '(9(/))' )

		READ( 216 ,  * ) aeroGrid, struMesh, forces
        IF( aeroGrid == 0 )THEN !---------
          aerGriOut = .FALSE.
        ELSE
          aerGriOut = .TRUE.
        ENDIF
        IF( struMesh == 0 )THEN !---------
          strMesOut = .FALSE.
        ELSE
          strMesOut = .TRUE.
        ENDIF
        IF( forces == 0 )THEN !-----------
          forcesOut = .FALSE.
        ELSE
          forcesOut = .TRUE.
        ENDIF
        
        CLOSE( UNIT = 216 )

        
        ENDSUBROUTINE Input_extras

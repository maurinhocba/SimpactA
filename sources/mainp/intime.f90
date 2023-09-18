      SUBROUTINE intime( )
!***********************************************************************
!
!     read initial conditions
!
!***********************************************************************
      USE ctrl_db, ONLY: ndime,nrotd,ndofn,neulr,npoin,initial_displacements,ndoft,itemp
      USE npo_db, ONLY: label,euler,coorc,velnp,tempe,dtemp,psia
      USE outp_db, ONLY :  iwrit
      USE c_input
      USE nsets_db

      IMPLICIT NONE
      REAL(KIND=8) :: facto = 57.295779513082320876798154814105      ! facto= 180d0/pi
      REAL    (kind=8) :: xg(9)
      INTEGER (kind=4) :: i,n, neu
      CHARACTER(len=mnam) :: sname
      LOGICAL found
      TYPE (nset), POINTER :: ns
      INTEGER (kind=4) chnode

      INTERFACE
        INCLUDE 'elemnt.h'
        INCLUDE 'coract.h'
      END INTERFACE

!----------------------------------- reading initial conditions

      CALL listen('INTIME')

      ! INITIAL dislacements are computed as INCREMENTAL over previous strategy.
      ! coora = coora + displ
      ! This assumes that the data is consistent with boundary conditions and
      ! constraints.
      ! Also some factors and data associated to slave (dependant)
      ! nodes should be recomputed (FF)

      IF (exists('INITIA')) THEN       !key word INITIAL CONDITIONS

        CALL listen('INTIME')      !read a card

        DO                             !loop for each type of variables
          IF(exists('ENDINI'))EXIT     !key word END_INITIAL_CONDITIONS found
                                                 !=> exit loop
          !     *** initial displacements

          IF (exists('DISPLA')) THEN   !key word DISPLACEMENTS
            ! initial displacements are applied at INITIAL_COMP/WRTPOS/ACTCD
            IF( .NOT.initial_displacements )THEN
              initial_displacements = .TRUE.
              coorc = 0d0
              IF(ndofn == 8) psia = 0d0
            END IF
                          !print header
            IF(iwrit == 1) WRITE(lures,"(//,6x,'Node   Initial X-disp.', &
                           &            '   Initial Y-disp.   Initial Z-disp.'/)",ERR=9999)

            CALL listen('INTIME')      !read a card
            DO                               !loop for nodal displacements
              IF(exists('ENDDIS'))EXIT       !key-word END_DISPLACEMENTS found
              xg(1:ndime) = param(2:ndime+1) !get displacements
              IF (exists('SET   ',i))THEN    !Associated to a set?
                sname = get_name(posin=i,stype='NSET')    !set name
                CALL nsdb_search (sname, found, ns)   !search if node set exists
                IF (found) THEN              !set found, go ahead
                  CALL ns_head (ns)          !point to first value in the set
                  DO                         !loop over the nodes in the set
                    IF ( end_ns(ns) ) EXIT   !last node processed, EXIT
                    n = get_label(ns)        !get nodal label
                    i = chnode(n)            !get internal numeration
                    coorc(1:ndime,i) =  xg(1:ndime) !store to be summed later
                    IF(iwrit == 1) WRITE(lures,"(i10,3e15.5)",ERR=9999) label(i), &
                                   xg(1:ndime)
                    CALL ns_next (ns)        !point to next node in the list
                  END DO
                ELSE          !set not found => ERROR
                  WRITE (lures, "(' Set ',A,'  not found')",ERR=9999) TRIM(sname)
                  CALL runend('Intime: Set not found ')
                END IF

              ELSE       ! displacements of a single node
                n  = INT(param(1))   !nodal label
                i = chnode(n)        !internal node
                coorc(1:ndime,i) =  xg(1:ndime) !store to be summed later
                IF(iwrit == 1) WRITE(lures,"(i10,3e15.5)",ERR=9999) label(i), &
                               xg(1:ndime)
              END IF
              CALL listen('INTIME')
            END DO
          END IF

          !     *** initial Euler angles

          IF (exists('LOCALS')) THEN    !key-word LOCAL_SYSTEMS
            ! initial Euler Angles are applied at INITIAL_COMP/EQNUMS
            neu = 2*ndime - 3
            IF(neulr == 9) THEN         ! 3-D problems
              IF(iwrit == 1)WRITE(lures,"(//,'   Node       Initial ', &
            & 'alpha      beta         gamma.'/)",ERR=9999)
             ELSE IF(neulr == 1)THEN     ! 2-D problems
               IF(iwrit == 1)WRITE(lures,"(//,'   Node       ANGLE'/)",ERR=9999)
             ELSE                        ! error
              CALL runend('INTIME: Check NEULR vs LOCAL_System')
            END IF

            CALL listen('INTIME')       !read a card
            DO                          !loop for nodal local systems
              IF(exists('ENDLOC'))EXIT  !key-word END_LOCALS => exit loop
              xg(1:neulr) = param(2:neulr+1)    !get local system
              IF (exists('SET   ',i))THEN       !data for a SET
                sname = get_name(posin=i,stype='NSET')       !set name
                CALL nsdb_search (sname, found, ns)  !search for the set of nodes
                IF (found) THEN                   !IF set found
                  CALL ns_head (ns)               !point to first node in the set
                  DO                              !loop over the node
                    IF ( end_ns(ns) ) EXIT        !if last node, exit loop
                    n = get_label(ns)             !nodal label
                    i = chnode(n)                 !internal number
                    euler(1:neu,i) = xg(1:neu)/facto  !associate data
                    IF( neu == 3 )euler(4:9,i) = 0d0  !associate data
                    IF(iwrit == 1) WRITE(lures,"(i7,3x,(3f20.16))",ERR=9999) & !ECHO
                                               label(i),xg(1:neu)
                    CALL ns_next (ns)             !point to next node
                  END DO
                ELSE                              !set not found => error
                  WRITE(lures,"(' Set ',A,'  not found')",ERR=9999) TRIM(sname)
                  CALL runend('Intime: Set not found ')
                END IF

              ELSE          !data for individual nodes
                n  = INT(param(1))   !node label
                i = chnode(n)        !internal number
                euler(1:neu,i) = xg(1:neu)/facto  !associate data
                IF( neu == 3 )euler(4:9,i) = 0d0  !associate data
                IF(iwrit == 1) WRITE(lures,"(i7,3x,3f20.16)",ERR=9999)  & !ECHO
                                                label(i),xg(1:neu)
              END IF
              CALL listen('INTIME')
            END DO  !END_LOCALS key-word found
          END IF                     !

          IF (exists('ADDISP')) THEN    !key-word AD_DISP
            ! initial additional displacements
            IF(ndofn == 8) THEN         !
              IF(iwrit == 1)WRITE(lures,"(//,'   Node       Add-Dx  Add-Dy'/)",ERR=9999)
            ELSE                        ! error
              CALL runend('INTIME: Check NDOFN vs Ad-Displacem')
            END IF
            CALL listen('INTIME')       !read a card
            DO                          !loop for nodal local systems
              IF(exists('ENDADD'))EXIT  !key-word END_AD_DIS => exit loop
              xg(1:2) = param(2:3)    !get local system
              IF (exists('SET   ',i))THEN       !data for a SET
                sname = get_name(posin=i,stype='NSET')       !set name
                CALL nsdb_search (sname, found, ns)  !search for the set of nodes
                IF (found) THEN                   !IF set found
                  CALL ns_head (ns)               !point to first node in the set
                  DO                              !loop over the node
                    IF ( end_ns(ns) ) EXIT        !if last node, exit loop
                    n = get_label(ns)             !nodal label
                    i = chnode(n)                 !internal number
                    psia(1:2,i) = xg(1:2)         !associate data
                    IF(iwrit == 1) WRITE(lures,"(i7,3x,(2f20.16))",ERR=9999) & !ECHO
                                               label(i),xg(1:2)
                    CALL ns_next (ns)             !point to next node
                  END DO
                ELSE                              !set not found => error
                  WRITE(lures,"(' Set ',A,'  not found')",ERR=9999) TRIM(sname)
                  CALL runend('Intime: Set not found ')
                END IF

              ELSE          !data for individual nodes
                n  = INT(param(1))   !node label
                i = chnode(n)        !internal number
                psia(1:2,i) = xg(1:2)        !associate data
                IF(iwrit == 1) WRITE(lures,"(i7,3x,3f20.16)",ERR=9999)  & !ECHO
                                                label(i),xg(1:2)
              END IF
              CALL listen('INTIME')
            END DO  !END_AD_DIS key-word found
          END IF                     !

          !     *** initial velocities

          IF (exists('VELOCI')) THEN   !key-word VELOCITIES (all DOFs)
                              ! print header (OMEGA missing ) (FF)
            IF(iwrit == 1)WRITE(lures,"(//,'   Node   Init X-velo. ', &
                         &      '  Init Y-velo.   Init   Z-velo.'/)",ERR=9999)

            CALL listen('INTIME')          !read a card
            DO                             !loop for data
              IF(exists('ENDVEL'))EXIT     !key-word END_VELOCITIES =>exit loop
              xg = param(2:nrotd+1)        !get data
              IF (exists('SET   ',i))THEN  !If data is for a SET
                sname = get_name(posin=i,stype='NSET')  !set name
                CALL nsdb_search (sname, found, ns)   !search if set exist
                IF (found) THEN            !if set found
                  CALL ns_head (ns)        !point to first node in the set
                  DO                       !loop over the nodes
                    IF ( end_ns(ns) ) EXIT !if last node, exit loop
                    n = get_label(ns)      !node label
                    i = chnode(n)          !internal number
                    velnp(1:nrotd,i) = xg(1:nrotd)  !associate data
                    IF(iwrit == 1) WRITE(lures,"(i8,6e12.4)",ERR=9999) label(i),  & !ECHO
                                            xg(1:nrotd)
                    CALL ns_next (ns)
                  END DO
                ELSE              !set not found =>error
                  WRITE (lures, "(' Set ',A,'  not found')",ERR=9999) TRIM(sname)
                  CALL runend('Intime: Set not found ')
                END IF

              ELSE                !data for individual node
                n  = INT(param(1))  !node label
                i = chnode(n)       !internal number
                velnp(1:nrotd,i) = xg(1:nrotd)    !associate data
                IF(iwrit == 1) WRITE(lures,"(i8,6e12.4)",ERR=9999) label(i),    & !ECHO
                                            xg(1:nrotd)
              END IF
              CALL listen('INTIME')    !read a card
            END DO        !END_VELOCITIES  card found
          END IF       !key-word VELOCITIES exist

          !     *** initial temperatures
          !    although not checked, this has only sense for Coupled Thermo-Mechanical analysis

          IF (exists('TEMPER')) THEN   !key word TEMPERATURE
                          !print header
            IF(iwrit == 1)THEN
              IF (ndoft == 1)WRITE(lures,"(//,6x,'Node  Initial Temper.',/)",ERR=9999)
              IF (ndoft == 2)WRITE(lures,"(//,6x,'Node  Initial TI  Initial TS ',/)",ERR=9999)
              IF (ndoft == 3)WRITE(lures,"(//,6x,'Node  Initial TN  Initial TI  Initial TS',/)",ERR=9999)
            END IF

            itemp = .TRUE.
            IF( .NOT.ASSOCIATED (tempe) )THEN
              ALLOCATE (tempe(ndoft,npoin), dtemp(ndoft,npoin))
              tempe = 0d0
            END IF
            CALL listen('INTIME')      !read a card
            DO                               !loop for nodal temperatures
              IF(exists('ENDTEM'))EXIT       !key-word END_TEMPERATURES found
              xg(1:ndoft) = param(2:ndoft+1) !get temperatures
              IF (exists('SET   ',i))THEN    !Associated to a set?
                sname = get_name(posin=i,stype='NSET')    !set name
                CALL nsdb_search (sname, found, ns)   !search if node set exists
                IF (found) THEN              !set found, go ahead
                  CALL ns_head (ns)          !point to first value in the set
                  DO                         !loop over the nodes in the set
                    IF ( end_ns(ns) ) EXIT   !last node processed, EXIT
                    n = get_label(ns)        !get nodal label
                    i = chnode(n)            !get internal numeration
                    tempe(1:ndoft,i) = xg(1:ndoft) !update
                    IF(iwrit == 1) WRITE(lures,"(i10,3e12.4)",ERR=9999) label(i),xg(1:ndoft)
                    CALL ns_next (ns)        !point to next node in the list
                  END DO
                ELSE              !set not found =>error
                  WRITE (lures, "(' Set ',a,'  not found')",ERR=9999)TRIM(sname)
                  CALL runend('Intime: Set not found ')
                END IF

              ELSE       ! temperatures of a single node
                n  = INT(param(1))            !nodal label
                i = chnode(n)  !internal node
                tempe(1:ndoft,i) = xg(1:ndoft)             !nodal
                IF(iwrit == 1) WRITE(lures,"(i10,3e12.4)",ERR=9999) label(i),xg(1:ndoft)
              END IF
              CALL listen('INTIME')
            END DO
            dtemp = 0d0 !tempe !... in time == 0d0 must be dtemp == 0d0
          END IF       !key-word TEMPERATURES exist

          !     *** initial Gaussian variables
          !     (presently not for all elements)

          IF (exists('GAUSSI',i)) THEN

            sname = get_name(posin=i,stype='ESET')       !set name

            IF(iwrit == 1) WRITE(lures,"(// &
                   &     '  Initial Gaussian variables '/ &
                   &     '    for element set: ',A /)",ERR=9999) TRIM(sname)   !echo

            !found = .FALSE.  !unnecessary
            CALL elemnt ('INIGAU', name=sname, flag2=found)
            IF(.NOT.found) WRITE (lures,"(' Set ',A, &
                 &    ' with Initial Internal Variables NOT FOUND')",ERR=9999) &
                 &      TRIM(sname)

          END IF

          CALL listen('INTIME')      !read a card

        END DO      !EXIT on END_INITIAL

      ELSE         !no initial data

        backs = .TRUE.       !one line back

      END IF

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE intime

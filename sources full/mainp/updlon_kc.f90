      SUBROUTINE updlon_kc (naris,nvelr,ndime)

      !   updates kinematic condition databases

      USE param_db,ONLY: mnam
      USE ifx_db
      USE nes_db
      USE nar_db
      USE ndp_db
      USE rve_db
      USE rfdp_db
      USE curv_db, ONLY : del_cur
      IMPLICIT NONE

      INTEGER (kind=4) :: naris,nvelr,ndime

      INTEGER (kind=4) :: i,j,k,lab,chnode,nn
      TYPE (ifx_nod), POINTER :: fix, afix
      TYPE (rve_set), POINTER :: rves,anter
      TYPE (rve_nod), POINTER :: rven,aven
      TYPE (rfdp_nod), POINTER :: rfn,arfn
      TYPE (nes_nod), POINTER :: mpcs,heads,tails
      LOGICAL :: keep,now
      CHARACTER (len=mnam) :: lbl

      !---------------FIXITIES--------------
      IF (nifx > 0) THEN  !fixed node exist
        nn = 0            !counter
        fix => ihead
        NULLIFY(afix)       !nullify auxiliar value
        NULLIFY(itail)      !nullify last node position
        DO i=1,nifx         !for each value
          lab = chnode(fix%ifix(1))   !associated label
          IF( lab > 0) THEN             !if node exist
            nn = nn + 1                 !updates number of fixed nodes
            afix => fix                 !point to last known value
            itail => afix               !update last node
            fix => fix%next             !go to next
            CYCLE
          ELSE IF( nn > 0 )THEN      !there are nodes in the list
            afix%next => fix%next       !jump node to next
            DEALLOCATE(fix)             !release memory
            fix => afix%next            !point to next node in the list
          ELSE                       !no nodes in the lista
            ihead => ihead%next         !update header
            DEALLOCATE(fix)             !release memory
            fix => ihead                !point no next node in the list
          END IF
        END DO
        nifx = nn                              !update dependencies
      END IF

      !---------------SLAVE DOFS -----------------------------------
      IF (nnes > 0) THEN  !slave node exist
        CALL ini_nes(heads,tails)
        nn = 0            !counter of slave DOFs
        DO i=1,nnes         !for each value
          IF( nesdat(i)%factor == 0d0 )THEN   !IF slave dof
            lab = chnode(nesdat(i)%node)      !check slave node exist
            keep =  lab > 0
            j = i                             !initializes number of dependencies
            CYCLE                             !go for master DOFs
          END IF
          lab = chnode(nesdat(i)%node)        !master node label
          IF( lab == 0 ) keep = .FALSE.       !master does not exist
          now =  i == nnes                     !end of list
          IF( .NOT.now ) now = nesdat(i+1)%factor == 0d0  !end for this node
          IF( now )THEN                        !end of dependencies
            IF( keep )THEN
              DO k=j,i
                ALLOCATE (mpcs)                   !get memory
                mpcs%node = nesdat(k)%node        !master (slave) node
                mpcs%ndof = nesdat(k)%ndof        !master (slave) DOF
                mpcs%factor = nesdat(k)%factor    !factor (0 indicates a slave)
                CALL add_nes (mpcs, heads, tails)   !add to list
                nn = nn + 1
              END DO
            END IF
          END IF
        END DO
        IF( nnes /= nn )THEN     !if mpcs have changed
          nnes = nn                !update size
          DEALLOCATE (nesdat)      !release memory for the array
          IF( nnes > 0 )THEN          !if mpcs remains
            ALLOCATE (nesdat(nnes))     !get memory for the array
            CALL store_nes(heads)        !transfer data from list to array
                                        ! and release memory
            CALL ini_nes(heads,tails)     !nullify pointers
          END IF
        END IF
      END IF

      !---------------NODES ON an ARISTA ---------------------------
      IF (naris > 0) THEN
        nn = 0
        DO i=1,naris
          keep = .TRUE.
          DO j=1,3
            lab = chnode(nardat(j,i))
            IF( lab == 0 )keep = .FALSE.
          END DO
          IF( keep )THEN
            nn = nn+1
            IF( nn == i )CYCLE
            nardat(:,nn) = nardat(:,i)
          END IF
        END DO
        naris = nn
      END IF

      !---------------DEPENDENT NODES ------------------------------
      IF (nndp > 0) THEN
        nn = 0
        DO i=1,nndp
          keep = .TRUE.
          DO j=1,2
            lab = chnode(ndpdat(j,i))
            IF( lab == 0 )keep = .FALSE.
          END DO
          IF( keep )THEN
            nn = nn+1
            IF( nn == i )CYCLE
            ndpdat(:,nn) = ndpdat(:,i)
          END IF
        END DO
        nndp = nn
      END IF

      !---------------PRESCRIBED VELOCITIES ------------------------
      IF (nvelr > 0) THEN

        rves => headv          !set of nodes
        NULLIFY(anter)
        DO                     !for each set
          IF( .NOT.ASSOCIATED(rves) )EXIT
          rven => rves%head      !point to first node of the set
          NULLIFY(aven)
          nn = 0                   !initializes
          DO j=1,rves%nrv           !loop over each node
            lab = chnode(rven%node)   !present label
            IF( lab == 0 ) THEN         !does not exist
              IF (.NOT.ASSOCIATED (aven)) THEN  !if first node
                rves%head => rven%next              !new head
                DEALLOCATE  (rven)           !release memory
                rven => rves%head
              ELSE
                aven%next => rven%next         !skip posic
                DEALLOCATE  (rven)           !release memory
                rven = aven%next
              END IF
              CYCLE
            END IF
            nn = nn+1
            aven => rven
            rven => rven%next
          END DO
          IF( nn /= 0 ) THEN
            rves%nrv = nn     !number of nodes kept
            tailv => rves
            anter => rves
            rves => rves%next
          ELSE                !no nodes exist, delete set
            lbl = rves%lc               !set label
            CALL del_rves (headv, tailv, anter, rves)   !delete set assoc to LC
            CALL del_cur (lbl)          !delete curve data associated LC
            nvelr = nvelr - 1           !correct counter
            IF( ASSOCIATED(anter) )THEN   !point to next set
              rves => anter%next
            ELSE
              rves => headv
            END IF
          END IF
        END DO
      END IF

      !---------------ROTATION FREE DEPENDENCIES -------------------
      IF (nrfdp > 0) THEN
        IF( planar )THEN
          k = ndime
        ELSE
          IF( ndime == 2 ) k = 3
          IF( ndime == 3 ) k = 6
        END IF
        rfn => rfd_head
        NULLIFY(arfn)
        nn = 0
        DO i=1,nrfdp
          lab = chnode(rfn%slave)
          IF( lab > 0 )THEN
            keep = .TRUE.
            DO j=1,k
              lab = rfn%lnods(j)
              IF( lab == 0 )EXIT
              lab = chnode(lab)
              IF( lab == 0 )THEN
                keep = .FALSE.
                EXIT
              END IF
            END DO
          END IF
          IF( keep )THEN
            arfn = rfn
            rfn => rfn%next
            nn = nn+1
          ELSE
            IF( ASSOCIATED(arfn) )THEN
              arfn%next = rfn%next
              DEALLOCATE(rfn%lnods)
              DEALLOCATE(rfn)
              rfn => arfn%next
            ELSE
              rfd_head => rfn%next
              DEALLOCATE(rfn%lnods)
              DEALLOCATE(rfn)
              rfn => rfd_head
            END IF
          END IF
        END DO
        nrfdp = nn
      END IF

      RETURN

      END SUBROUTINE updlon_kc

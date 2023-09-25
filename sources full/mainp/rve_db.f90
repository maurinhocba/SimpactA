      MODULE rve_db
        USE param_db,ONLY: mnam
        IMPLICIT NONE
        SAVE

        ! Derived type for rigid (prescribed) velocities

        ! node in a list
        TYPE rve_nod
          INTEGER (kind=4) :: node        !node label
          REAL (kind=8) :: v(8)           !prescribed velocities
          TYPE (rve_nod), POINTER :: next   !next position
        END TYPE rve_nod

        ! set (list of nodes)
        TYPE rve_set
          CHARACTER (len=mnam) :: lc                !associated curve label
          INTEGER (kind=4) :: nrv                   !number of nodes in the set
          INTEGER(kind=4) :: dspflg                 !displacement flag
          REAL (kind=8) :: factor                   !curve factor
          TYPE (rve_nod), POINTER :: head,tail      !pointers to first and last
          TYPE (rve_set), POINTER :: next           !pointer to next set
        END TYPE rve_set

        ! global pointers
        TYPE (rve_set), POINTER :: headv,tailv   !pointers to first & last sets

      CONTAINS

        ! functions for rve nodes

        SUBROUTINE ini_rven (head, tail)
          !initialize a list of nodes

          !Dummy arguments
          TYPE (rve_nod), POINTER :: head, tail

          NULLIFY (head, tail)

        END SUBROUTINE ini_rven

        SUBROUTINE add_rven (new, head, tail)
          !This subroutine adds data to the end of the list of nodes
          !Dummy arguments
          TYPE (rve_nod), POINTER :: new, head, tail

          !Check if a list is empty
          IF (.NOT. ASSOCIATED (head)) THEN
            !list is empty, start it
            head => new
            tail => new
            NULLIFY (tail%next)

          ELSE
            !add a node to the list
            tail%next => new
            NULLIFY (new%next)
            tail => new

          END IF
          RETURN
        END SUBROUTINE add_rven

        ! functions for rve sets

        SUBROUTINE ini_rves (head, tail)
          !initialize a list

          !Dummy arguments
          TYPE (rve_set), POINTER :: head, tail

          NULLIFY (head, tail)

        END SUBROUTINE ini_rves

        SUBROUTINE add_rves (new, head, tail)
          !This subroutine adds a set to the end of the list of sets
          !Dummy arguments
          TYPE (rve_set), POINTER :: new, head, tail

          !Check if a list is empty
          IF (.NOT. ASSOCIATED (head)) THEN
            !list is empty, start it
            head => new
            tail => new
            NULLIFY (tail%next)

          ELSE
            !add a set to the list
            tail%next => new
            NULLIFY (new%next)
            tail => new

          END IF
          RETURN
        END SUBROUTINE add_rves

        SUBROUTINE srch_rves (head, anter, posic, lc, found)
          !This subroutine searches for a set associated to curve LC
          !Dummy arguments
          LOGICAL :: found         ! Answer
          CHARACTER (len=*) :: lc  ! curve label
          TYPE (rve_set), POINTER :: head, anter, posic

          found = .FALSE.                   !initializes
          NULLIFY (posic,anter)
          !Check if a list is empty
          IF (ASSOCIATED (head)) THEN
            posic => head                   !point to head of the list
            DO
              IF(TRIM(posic%lc) == TRIM(lc)) THEN   !compare with associated curve
                found = .TRUE.              !found
                EXIT                        !leave search
              END IF
              IF (ASSOCIATED(posic%next) ) THEN     !more sets to check
                anter => posic                      !keep previous
                posic => posic%next                 !point to next
              ELSE
                EXIT                        !no more sets => exit
              END IF
            END DO
          ENDIF
          IF (.NOT.found) NULLIFY (posic,anter)   !not found => nullify pointers
          RETURN
        END SUBROUTINE srch_rves

        SUBROUTINE del_rves (head, tail, anter, posic)
          !deletes a set pointed with posic
          TYPE (rve_set), POINTER :: head, tail, anter, posic

          IF (.NOT.ASSOCIATED (anter)) THEN  !if first set
            head => posic%next               !new head
          ELSE
            anter%next => posic%next         !skip posic
          END IF
          ! if posic == tail    (last set)
          IF (.NOT.ASSOCIATED (posic%next) ) tail => anter  !new last set
          CALL dalloc_rves (posic)           !release memory
          !NULLIFY (anter)                    !both anter & posic are null now
        END SUBROUTINE del_rves

        SUBROUTINE dalloc_rves (rves)
          ! deallocates a rve set (release memory)
          TYPE (rve_set), POINTER :: rves
          TYPE (rve_nod), POINTER :: rven, rvenaux

          rven => rves%head    !point to first node
          DO
            IF (.NOT.ASSOCIATED (rven) ) EXIT
            rvenaux => rven%next    !keep next pointer
            DEALLOCATE (rven)       !release memory of the node
            rven => rvenaux         !point to next
          END DO

          DEALLOCATE (rves)         !release rest of the vars. lc nrv & factor
          RETURN
        END SUBROUTINE dalloc_rves

        SUBROUTINE store_rve (head,vel,nods,ndofn)
          !This subroutine stores and prints the data

          !Dummy argument
          TYPE (rve_nod), POINTER :: head
          INTEGER (kind=4) :: ndofn, nods(:)
          REAL (kind=8) :: vel(:,:)

          !Local variables
          TYPE (rve_nod), POINTER :: ptr
          INTEGER :: n

          IF (ASSOCIATED(head))THEN
            ptr => head

            n = 0
            DO
              n = n + 1
              nods(n) = ptr%node
              vel(1:ndofn,n) = ptr%v(1:ndofn)
              ptr => ptr%next
              IF (.NOT.ASSOCIATED(ptr)) EXIT
            END DO
          ENDIF
        END SUBROUTINE store_rve
      END MODULE rve_db

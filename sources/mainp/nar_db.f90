      MODULE nar_db
        IMPLICIT NONE
        SAVE

        ! Derived type for a list containing restricted nodes
        ! (nodes on "arista")
        TYPE nar_nod
          INTEGER (kind=4) :: slave,mast01,mast02
          TYPE (nar_nod), POINTER :: next
        END TYPE nar_nod
        !
        ! array containing information for direct access
        INTEGER (kind=4), POINTER :: nardat(:,:)

      CONTAINS
        SUBROUTINE ini_nar (head, tail)
          !initialize a list

          !Dummy arguments
          TYPE (nar_nod), POINTER :: head, tail

          NULLIFY (head, tail)

        END SUBROUTINE ini_nar

        SUBROUTINE add_nar (new, head, tail)
          !This subroutine adds data to the end of the list
          !Dummy arguments
          TYPE (nar_nod), POINTER :: new, head, tail

          !Check if a list is empty
          IF (.NOT. ASSOCIATED (head)) THEN
            !list is empty, start it
            head => new
            tail => new
            NULLIFY (tail%next)

          ELSE
            !add to the list
            tail%next => new
            NULLIFY (new%next)
            tail => new

          ENDIF
        END SUBROUTINE add_nar

        SUBROUTINE store_nar (head)
          ! Transfer data from list (head) into array (nardat)
          ! and release memory of list

          !Dummy argument
          TYPE (nar_nod), POINTER :: head

          !Local variables
          TYPE (nar_nod), POINTER :: ptr,prev
          INTEGER :: n

          ptr => head    !point to first

          n = 0          !initializes array counter
          DO
            IF (.NOT.ASSOCIATED(ptr)) EXIT
            n = n + 1                 !increase counter
            nardat(1,n) = ptr%slave   !transfer data
            nardat(2,n) = ptr%mast01
            nardat(3,n) = ptr%mast02
            prev => ptr               !keep present data
            ptr => ptr%next           !point to next
            DEALLOCATE (prev)         !release memory already transfered
          END DO

          RETURN
        END SUBROUTINE store_nar
      END MODULE nar_db

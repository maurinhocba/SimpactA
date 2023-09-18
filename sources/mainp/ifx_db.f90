MODULE ifx_db
  IMPLICIT NONE
  SAVE

  ! Derived type for fixity data
  ! list
  TYPE ifx_nod
    INTEGER (kind=4) :: ifix(8)          !1:label 2:8 fix codes
    TYPE (ifx_nod), POINTER :: next
  END TYPE ifx_nod
  INTEGER (kind=4) :: nd1=0,nifx=0       !NDOFN+(1)   Number of Fixities
  TYPE (ifx_nod), POINTER :: ihead,itail

CONTAINS

  SUBROUTINE ini_ifx (head, tail)
    !initialize a dependent nodes list

    !Dummy arguments
    TYPE (ifx_nod), POINTER :: head, tail

    NULLIFY (head, tail)

  END SUBROUTINE ini_ifx

  SUBROUTINE add_ifx (new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (ifx_nod), POINTER :: new, head, tail

    !Check if a list is empty
    new%ifix(nd1+2:8) = 0
    IF (.NOT. ASSOCIATED (head)) THEN
      !list is empty, start it
      head => new
      tail => new
      NULLIFY (tail%next)

    ELSE
      !add data to the list
      tail%next => new
      NULLIFY (new%next)
      tail => new

    END IF
  END SUBROUTINE add_ifx

  SUBROUTINE srch_ifx (head, anter, posic, node, found)
    !This subroutine searches for a ifx data for a "node"
    !Dummy arguments
    LOGICAL :: found
    INTEGER (kind=4) :: node
    TYPE (ifx_nod), POINTER :: head, anter, posic

    found = .FALSE.                         !initializes flag and pointers
    NULLIFY (posic,anter)                   !
    !Check if a list is empty
    IF (ASSOCIATED (head)) THEN
      posic => head                         !point to first
      DO
        IF(posic%ifix(1) == node) THEN      !check label
          found = .TRUE.                    ! O.K.
          EXIT
        END IF
        IF (ASSOCIATED(posic%next) ) THEN   !if there are more points
          anter => posic                    !keep previous position
          posic => posic%next               !point to next
        ELSE
          EXIT                              !no more points EXIT
        END IF
      END DO
    END IF
    IF (.NOT.found) NULLIFY (posic,anter)   !not found, nullify pointers
    RETURN
  END SUBROUTINE srch_ifx

  SUBROUTINE del_ifx (head, tail, anter, posic)

    !This subroutine deletes data pointed with posic and nullifies posic & anter

    TYPE (ifx_nod), POINTER :: head, tail, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN   !see if we are at the beginning
      head => posic%next             !then modifies head pointer
    ELSE
      anter%next => posic%next       !else skip posic
    END IF
             ! if posic == tail
    IF (.NOT.ASSOCIATED (posic%next) ) tail => anter
    DEALLOCATE (posic)               !release memory
    NULLIFY (anter)                  !both posic and anter are null now
    RETURN
  END SUBROUTINE del_ifx

  SUBROUTINE dump_ifx
  !dumps data
  IMPLICIT NONE

    !Local variables
    TYPE (ifx_nod), POINTER :: ifx
    !INTEGER :: i,n,node

    WRITE (50,ERR=9999) nd1,nifx       !nd1 = NDOFN+(1) , nifx = Nodes with resctrictions
    IF (nifx > 0) THEN        !if restrictions exist

      ifx => ihead             !point to first

      DO
        WRITE (50,ERR=9999) ifx%ifix(1:nd1+1)     !writes label and codes
        ifx => ifx%next                  !point to next
        IF (.NOT.ASSOCIATED(ifx)) EXIT
      END DO
    END IF

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE dump_ifx

  SUBROUTINE rest_ifx
    !restores variables from disk
    !Local variables
    TYPE (ifx_nod), POINTER :: ifx
    INTEGER :: i

    !Initialize empty list
    CALL ini_ifx(ihead,itail)

    READ (51) nd1,nifx    !nd1 = NDOFN+(1) , nifx = Nodes with resctrictions
                          !read list
    DO i = 1,nifx
      ALLOCATE (ifx)                   !get memory
      READ (51) ifx%ifix(1:nd1+1)      !read data
      CALL add_ifx( ifx, ihead, itail )  !add to list
    END DO

    RETURN
  END SUBROUTINE rest_ifx

END MODULE ifx_db

      MODULE ele22_db
        USE param_db,ONLY: mnam
        USE mat_dba, ONLY: section, sect_search, mater, psecs
        IMPLICIT NONE

        ! Derived type for a set of HEAT3 (ELE22) elements
        ! 3D thermal elements
        TYPE ele22_set
          CHARACTER (len=mnam) :: sname   ! set name
          INTEGER (kind=4) :: nelem    ! number of elements
          INTEGER (kind=4) :: nreqs=0  ! number of GP for hist. output
          INTEGER (kind=4) :: ngaus=2  ! integration order
          INTEGER (kind=4) :: nnode=8  ! number of nodes
          INTEGER (kind=4) :: narch    ! number of output unit
          INTEGER (kind=4) :: lcur=0   ! number of curve controlling heat gen
          REAL (kind=8) :: posgp(2)
          REAL (kind=8) :: weigp(2)
          REAL (kind=8) :: hgs
          INTEGER (kind=4), POINTER :: lnods(:,:) ! element connectivity array
          INTEGER (kind=4), POINTER :: matno(:)   ! material number
          INTEGER (kind=4), POINTER :: ngrqs(:) !
          REAL (kind=8), POINTER    :: shape(:,:)
          REAL (kind=8), POINTER    :: deriv(:,:,:)
          TYPE (ele22_set), POINTER :: next
        END TYPE ele22_set
        TYPE (ele22_set), POINTER :: head
        TYPE (ele22_set), POINTER :: tail

      CONTAINS
        SUBROUTINE ini_ele22 (head, tail)
          !initialize a list of ELE22 sets

          TYPE (ele22_set), POINTER :: head, tail

          NULLIFY (head, tail)

        END SUBROUTINE ini_ele22

        SUBROUTINE add_ele22 (new, head, tail)
          !This subroutine adds data to the end of the list
          !Dummy arguments
          TYPE (ele22_set), POINTER :: new, head, tail

          !Check if a list is empty
          IF (.NOT. ASSOCIATED (head)) THEN
            !list is empty, start it
            head => new
            tail => new
            NULLIFY (tail%next)

          ELSE
            !add a segment to the list
            tail%next => new
            NULLIFY (new%next)
            tail => new

          ENDIF
        END SUBROUTINE add_ele22

        SUBROUTINE srch_ele22 (head, anter, posic, name, found)
          !This subroutine searches for a set named "name"
          !Dummy arguments
          LOGICAL :: found
          CHARACTER (len=*) :: name ! set name
          TYPE (ele22_set), POINTER :: head, anter, posic

          found = .FALSE.
          NULLIFY (posic,anter)
          !Check if a list is empty
          IF (ASSOCIATED (head)) THEN
            posic => head
            DO
              IF(TRIM(posic%sname) == TRIM(name)) THEN
                found = .TRUE.
                EXIT
              END IF
              IF (ASSOCIATED(posic%next) ) THEN
                anter => posic
                posic => posic%next
              ELSE
                EXIT
              END IF
            END DO
          ENDIF
          IF (.NOT.found) NULLIFY (posic,anter)
        END SUBROUTINE srch_ele22

        SUBROUTINE del_ele22 (head, anter, posic)
          !This subroutine deletes a set pointed with posic
          TYPE (ele22_set), POINTER :: head, anter, posic

          IF (.NOT.ASSOCIATED (anter)) THEN
            head => posic%next
          ELSE
            anter%next => posic%next
          ENDIF
          NULLIFY (posic,anter)
        END SUBROUTINE del_ele22

  INCLUDE 'acvd22.fi'
  INCLUDE 'code22.fi'
  INCLUDE 'comm22.fi'
  INCLUDE 'elmd22.fi'
  INCLUDE 'gaus22.fi'
  INCLUDE 'gauste.fi'
  INCLUDE 'heat22.fi'
  INCLUDE 'luma22.fi'
  INCLUDE 'mase22.fi'
  INCLUDE 'poin22.fi'
  INCLUDE 'updl22.fi'
      END MODULE ele22_db


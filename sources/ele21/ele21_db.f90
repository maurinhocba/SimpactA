MODULE ele21_db
  USE param_db,ONLY: mnam
  USE mat_dba, ONLY: section, sect_search, mater, psecs
  IMPLICIT NONE

  ! Derived type for a set of ELE21 elements
  ! 2D thermal elements  (plane and axisymmetric)
  TYPE ele21_set
    CHARACTER (len=mnam) :: sname  ! set name
    INTEGER (kind=4) :: nelem    ! number of elements
    INTEGER (kind=4) :: nreqs=0  ! number of GP for hist. output
    INTEGER (kind=4) :: ngaus=2  ! integration order
    INTEGER (kind=4) :: nnode=4 ! maximum number of nodes
    INTEGER (kind=4) :: narch    ! number of output unit
    INTEGER (kind=4) :: lcur=0   ! number of curve controlling heat generation
    REAL (kind=8) :: posgp(2)
    REAL (kind=8) :: weigp(2)
    REAL (kind=8) :: hgs
    INTEGER (kind=4), POINTER :: lnods(:,:) ! element connectivity array
    INTEGER (kind=4), POINTER :: matno(:)   ! material number
    INTEGER (kind=4), POINTER :: ngrqs(:) !
    REAL (kind=8), POINTER    :: shape(:,:)
    REAL (kind=8), POINTER    :: deriv(:,:,:)
    TYPE (ele21_set), POINTER :: next
  END TYPE ele21_set
  TYPE (ele21_set), POINTER :: head
  TYPE (ele21_set), POINTER :: tail

CONTAINS

  SUBROUTINE ini_ele21 (head, tail)
    !initialize a list of ELE21 sets

    TYPE (ele21_set), POINTER :: head, tail

    NULLIFY (head, tail)

  END SUBROUTINE ini_ele21

  SUBROUTINE add_ele21 (new, head, tail)
    !This subroutine adds data to the end of the list
    !Dummy arguments
    TYPE (ele21_set), POINTER :: new, head, tail

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
  END SUBROUTINE add_ele21

  SUBROUTINE srch_ele21 (head, anter, posic, name, found)
    !This subroutine searches for a set named "name"
    !Dummy arguments
    LOGICAL :: found
    CHARACTER (len=*) :: name ! set name
    TYPE (ele21_set), POINTER :: head, anter, posic

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
  END SUBROUTINE srch_ele21

  SUBROUTINE del_ele21 (head, anter, posic)
    !This subroutine deletes a set pointed with posic
    TYPE (ele21_set), POINTER :: head, anter, posic

    IF (.NOT.ASSOCIATED (anter)) THEN
      head => posic%next
    ELSE
      anter%next => posic%next
    ENDIF
    NULLIFY (posic,anter)
  END SUBROUTINE del_ele21

  INCLUDE 'acvd21.fi'
  INCLUDE 'code21.fi'
  INCLUDE 'comm21.fi'
  INCLUDE 'deriv1.fi'
  INCLUDE 'elmd21.fi'
  INCLUDE 'gaus21.fi'
  INCLUDE 'heat21.fi'
  INCLUDE 'luma21.fi'
  INCLUDE 'poin21.fi'
  INCLUDE 'shape1t.fi'
  INCLUDE 'updl21.fi'

END MODULE ele21_db

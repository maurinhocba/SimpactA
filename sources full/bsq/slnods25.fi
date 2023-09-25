SUBROUTINE slnods25(heade,lnods,nelem)
!***********************************************************************
!*** Extracts conectivities into the nodal set
!***********************************************************************
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4), POINTER :: lnods(:,:)   !Conectivities of the set
  INTEGER(kind=4), INTENT (IN) :: nelem
  TYPE(ele25),POINTER:: heade      !pointer to first element
  !--- Local Variables
  INTEGER(kind=4), PARAMETER :: nnode = 4
  INTEGER(kind=4):: ielem,nn
  TYPE(ele25),POINTER:: e


  IF( ASSOCIATED(lnods) )THEN
    nn    = SIZE(lnods,1)
    ielem = SIZE(lnods,2)
    IF( ielem /= nelem .OR. nn /= nnode)THEN
      DEALLOCATE(lnods)
      ALLOCATE(lnods(nnode,nelem))
    END IF
  ELSE
    ALLOCATE(lnods(nnode,nelem))
  END IF
  ielem = 0                           !Initiates counter
  e => heade                          !Point to first element for each element in the set
  DO
    IF (.NOT.ASSOCIATED(e)) EXIT      !exit from the loop
    ielem = ielem + 1

    lnods(1:4,ielem) = e%lnods(1:4)   !transfer data

    e => e%next                       !point to next element
  END DO

RETURN
END SUBROUTINE slnods25

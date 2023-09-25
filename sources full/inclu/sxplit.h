SUBROUTINE sxplit(flag,velfi,eqres,str_type,lambd,chkin,dl2)
!********************************************************************
!*** time stepping routine
!********************************************************************
IMPLICIT NONE

  !-- Dummy variables
  REAL(kind=8), INTENT (IN OUT) :: dl2       !norm of incremental displacements (squared)
  REAL(kind=8), INTENT (IN OUT) :: velfi     !factor for prescribed displacement increment
  REAL(kind=8), INTENT (IN OUT) :: eqres(:), &  !non-equilibrated nodal forces
                                   lambd        !present load factor
  INTEGER(kind=4), INTENT(IN) :: str_type
  LOGICAL, INTENT(IN OUT) :: flag,chkin

END SUBROUTINE sxplit

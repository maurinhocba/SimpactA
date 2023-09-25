      SUBROUTINE dampin(neq,damp,veloc,dtime)
!****************************************************************
!
!     includes damping diminishing the velocities
!
!****************************************************************
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: neq
      REAL (kind=8),INTENT(IN) :: dtime,damp(:)
      REAL (kind=8),INTENT(IN OUT) :: veloc(:)

      END SUBROUTINE dampin

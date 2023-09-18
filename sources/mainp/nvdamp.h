      SUBROUTINE nvdamp(neq,damp,veloc,ddisp)
!****************************************************************
!
!     includes non-viscous local damping (changing residual forces)
!
!****************************************************************
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: neq
      REAL (kind=8),INTENT(IN) :: damp(:),veloc(:)
      REAL (kind=8),INTENT(IN OUT) :: ddisp(:)

      END SUBROUTINE nvdamp

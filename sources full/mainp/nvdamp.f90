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

      INTEGER (kind=4) :: ieq

      DO ieq=1,neq
        ! F_d = - a |F| sgn(v)
        ! below changed sign, ddisp are negative reverted residual forces
        ddisp(ieq) = ddisp(ieq) + damp(ieq)*SIGN(ddisp(ieq),veloc(ieq))
      END DO

      RETURN
      END SUBROUTINE nvdamp

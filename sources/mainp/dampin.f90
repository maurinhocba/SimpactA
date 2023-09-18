      SUBROUTINE dampin(neq,damp,veloc,dtime)
!****************************************************************
!
!     includes damping diminishing the velocities
!
!****************************************************************
!$ USE omp_lib
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: neq
      REAL (kind=8),INTENT(IN) :: dtime,damp(:)
      REAL (kind=8),INTENT(IN OUT) :: veloc(:)

      INTEGER (kind=4) :: ieq
      REAL    (kind=8) :: factor

 !$OMP DO PRIVATE(factor)
      DO ieq=1,neq
        factor = 1d0 - 2d0*damp(ieq)*dtime
        veloc(ieq) = factor*veloc(ieq)
      END DO
 !$OMP END DO

      RETURN
      END SUBROUTINE dampin

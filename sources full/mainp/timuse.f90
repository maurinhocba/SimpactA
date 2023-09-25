      SUBROUTINE timuse(cpu)

!     READ SYSTEM CPU TIME
      USE IFPORT
      IMPLICIT NONE
      REAL (kind=8),INTENT(OUT) :: cpu
      REAL (kind=8), SAVE :: cpuold = 0d0

       !CALL date_and_time (values=dt)     !too much data perhaps
       !CALL cpu_time( cpu )               !not good for multiple threads
       CALL clockx( cpu )  !wall-clock time in microseconds
       cpu = cpu*1d-6      !transform to seconds
       IF ( cpu < cpuold ) cpu = NINT(cpuold/86400)*86400 + cpu  !check change of day
       cpuold = cpu        !update

      RETURN
      END SUBROUTINE timuse

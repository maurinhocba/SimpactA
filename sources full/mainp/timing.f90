      SUBROUTINE timing (k,ind)

      ! initializes and add CPU times for different tasks

      USE outp_db, ONLY : time
      IMPLICIT NONE
      INTEGER (kind=4),INTENT(IN) :: k,ind

      REAL (kind=8) :: tcpu

     ! 1   :Initial computations
     ! 2   :Restart input-output
     ! 4   :Contact data input
     ! 5   :Gaussian variables intializat
     ! 6   :Data input
     ! 7   :Initial conditions
     ! 8   :Initial stress state
     ! 9   :Nodal load vector comput
     !10   :Lumped mass computation
     !11   :Mesh refinament evaluation
     !12   :Dumping of springback data
     !13   :Time integration
     !14   :Residual forces computation
     !15   :Results output
     !16   :Contact computations
     !17   :Time increment evaluation
     !18   :File closing
     !19   :Total process time
     ! timeo :DYNAMIC ANALYSIS

      CALL timuse (tcpu)                          !get clock time
      IF (ind == 1) THEN                       !start task
        time(k+20) = tcpu                         !initial time
      ELSE                                     !end task
        time(k) = time(k) + tcpu - time(k+20)     !adds elapsed time
      END IF
      RETURN

      END SUBROUTINE timing

SUBROUTINE timerst(init,tmrst,actrst)
!------------------------------------------------------------------------------------
! determines if a restart file must be dumped
! also initializes DTIMO for comparison with TMRST
!------------------------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN) :: tmrst   ! period for restart file dumping
  LOGICAL,INTENT(IN) :: init     !.TRUE. if a new period starts now
  LOGICAL,INTENT(OUT) :: actrst  !.TRUE. if a restart file must be dumped
  !--- Local variables
  REAL(kind=8),SAVE :: dtimo=0d0 !set when INIT = .TRUE.
  REAL(kind=8) :: stim

  !--- compute elapsed time or initializes
  IF (init) THEN
    CALL timuse(dtimo)   !Initialize 'dtimo'
    stim = 0d0             !Initialize the elapsed CPU time [sec]
  ELSE
    CALL timuse(stim)    !compute the CPU time
    stim = stim - dtimo    !elapsed CPU time since last dumping
    !--- determines if a restart file must be dumped
    IF (stim > DBLE(tmrst)) THEN  !elapsed time is larger than TMRST
      actrst = .TRUE.             !restart file must be dumped
      CALL timuse(dtimo)        !Initialize 'dtimo'
    ELSE
      actrst = .FALSE.            !do not dump
    END IF
  END IF

RETURN
END SUBROUTINE timerst

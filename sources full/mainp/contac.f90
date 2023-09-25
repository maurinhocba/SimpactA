SUBROUTINE contac(task,iwrit,dtcal,ttime,velnp,maxve)

  !     main contac routine

  USE ctrl_db, ONLY: mconts, nconp, numct, ctype, cactive, top, bottom
  USE outp_db, ONLY: nreqc,nprqc
  IMPLICIT NONE

  !        Dummy arguments
  CHARACTER(len=*),INTENT(IN) :: task  !task to perform
  INTEGER (kind=4),INTENT(IN) :: iwrit
  INTEGER (kind=4),INTENT(IN), OPTIONAL :: maxve
  REAL (kind=8),INTENT(IN), OPTIONAL :: dtcal,ttime,velnp(:,:)
  !        Local variables
  INTEGER (kind=4) :: i

  INTERFACE
    INCLUDE 'conta1.h'
    INCLUDE 'conta2.h'
    INCLUDE 'conta3.h'
    INCLUDE 'dbea3d.h'
  END INTERFACE

  IF( TRIM(task) == 'DUMPIN' .OR.  &
      TRIM(task) == 'RESTAR' .OR.  &
      TRIM(task) == 'UPDLON')THEN

    IF(cactive(1)) CALL conta1(task,dtcal,ttime,iwrit,velnp,maxve)
    IF(cactive(2)) CALL conta2(task,dtcal,ttime,iwrit,velnp,maxve)
    IF(cactive(3)) CALL conta3(task,dtcal,ttime,iwrit,velnp,maxve)
    IF(cactive(4)) CALL dbea3d(task,ttime,iwrit)

    IF(TRIM(task) == 'UPDLON')THEN
      !  check
      numct = 0                 !initializes
      DO i=1,mconts                 !for each contact algorithm
        cactive(i) = nconp(i) > 0
        IF( cactive(i) ) THEN       !if there are contact pairs using this alg.
          numct = numct + 1         !increase the number of cont algors used
          ctype (numct) = i         !keep for this set the element type
        END IF
      END DO
      IF( numct == 0 ) THEN
        nreqc = 0
        IF( ASSOCIATED(nprqc) )DEALLOCATE(nprqc)
        bottom = .FALSE. ; top = .FALSE.
      END IF
    END IF

  ELSE

    DO i = 1, numct

      SELECT CASE (ctype(i))

      CASE (1)
        CALL conta1(task,dtcal,ttime,iwrit,velnp,maxve)
      CASE (2)
        CALL conta2(task,dtcal,ttime,iwrit,velnp,maxve)
      CASE (3)
        CALL conta3(task,dtcal,ttime,iwrit,velnp,maxve)
      CASE (4)
        CALL dbea3d(task,ttime,iwrit)

      END SELECT
    END DO
  END IF
  RETURN
END SUBROUTINE contac

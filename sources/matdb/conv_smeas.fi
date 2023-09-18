SUBROUTINE conv_smeas (curve,npts,young,straim,stresm)

  !  CONVersion between Stress MEASures

  IMPLICIT NONE
  CHARACTER (len=*),INTENT(IN) :: straim, stresm
  INTEGER (kind=4), INTENT(IN) :: npts
  REAL (kind=8),INTENT(IN) :: young
  REAL (kind=8),INTENT(INOUT) :: curve(:,:)

  INTEGER (kind=4) :: n
  REAL (kind=8) :: dx,dy,eel

  ! if strain measure used in the input is not eff.plastic strain
  IF (TRIM(straim) == 'LOGARI') THEN
    ! logarithmic strain (total) were input
    DO n = 1, npts
      eel = curve(2,n)/young
      curve(1,n) = curve(1,n) - eel
    END DO
  ELSEIF (TRIM(straim) == 'ENGINE') THEN
    ! logarithmic strain (total) were input
    DO n = 1, npts
      eel = curve(2,n)/young
      curve(1,n) = LOG(curve(1,n)+1d0) - eel
    END DO
  END IF
  ! error in the calculation of elastic strains due to different
  ! stress measures is neglected in the calculations above

  ! if stress measure used in the input is not Cauchy stress
  IF (TRIM(stresm) == 'NOMINA') THEN
    ! nominal stresses
    DO n = 1, npts
      curve(2,n) = curve(2,n) *EXP(curve(1,n))
    END DO
  END IF

  ! calculate slopes
  DO n = 1, npts-1
    dx = curve(1,n+1) - curve(1,n)
    dy = curve(2,n+1) - curve(2,n)
    curve(3,n) = dy/dx
  END DO
  curve(3,npts) = curve(3,npts-1)

  RETURN

END SUBROUTINE conv_smeas

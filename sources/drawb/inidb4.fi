SUBROUTINE inidb4 ( nstatold, nstini, fdat, fdaini, nrt)
  !-----------------------------------------------------------------
  ! saving initial intersection
  !-----------------------------------------------------------------

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER (kind=4) :: nrt, nstatold (:), nstini (:)
  REAL (kind=8) :: fdat (:, :), fdaini (:, :)
  !-----------------------------------------------------------------

  !   l o c a l   v a r i a b l e s

  INTEGER :: i, k

  k = 0
  DO i = 1, nrt
    IF ( nstatold (i) == 1 ) THEN
      k = k + 1
      nstini (k) = i
      fdaini (1:2,k) = fdat (1:2,i)
    END IF
  END DO

  RETURN
END SUBROUTINE inidb4

LOGICAL FUNCTION nonear (px, x, ndime,toli)

!  checks point p against extremities of x()

  IMPLICIT NONE

  INTEGER(kind=4) :: ndime
  REAL (kind=8) :: px(ndime), x(ndime,3)
  REAL (kind=8), OPTIONAL :: toli
  ! here there is problema with the tolerance, it should consider coordinates dimensions
  REAL (kind=8), PARAMETER :: tol1=0.0001

  INTEGER(kind=4) :: i
  REAL(kind=8) :: xmax,xmin,tol

  nonear = .FALSE.
  IF( PRESENT(toli))THEN
    tol = toli 
  ELSE
    tol = tol1
  END IF
  DO i=1,ndime
    xmin = MINVAL(x(i,1:3))       !minimum value for this coordinate
    xmax = MAXVAL(x(i,1:3))       !maximum value for this coordinate

    IF( px(i) < xmin-tol .OR. px(i) > xmax+tol )THEN
      nonear = .TRUE.
      EXIT
    END IF
  END DO

  RETURN

END FUNCTION nonear

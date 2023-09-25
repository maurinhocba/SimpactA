FUNCTION det2 (x1,x2,x3,ndime,keep)

  !  calculates determinant (area)

  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: keep
  INTEGER (kind=4), INTENT(IN) :: ndime
  REAL (kind=8), INTENT(IN) :: x1(ndime), x2(ndime), x3(ndime)
  REAL (kind=8) :: det2

  !local variables
  REAL (kind=8) ::  r(ndime), s(ndime), t(ndime)
  REAL (kind=8), SAVE :: tn(3)     !triangle normal for 3-D problems

  r = x2 - x1    !first side of the triangle
  s = x3 - x1    !second side of the triangle

  IF(ndime == 2)THEN                 !for 2-d problems
    det2 = (r(1)*s(2) - r(2)*s(1))

  ELSE                               ! (ndime == 3)
    CALL vecpro( r, s, t)            ! normal vector
    IF( keep )THEN                   ! compute normal vector
      CALL vecuni(ndime,t,det2)   ! unit vector and twice the element area
      tn = t                         ! keep positive normal to compare
    ELSE
      det2 = DOT_PRODUCT(t,tn)       ! proyect over positive normal (existent)
    END IF
  END IF

  RETURN

END FUNCTION det2

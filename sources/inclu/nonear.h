      LOGICAL FUNCTION nonear (px, x, ndime,toli)

      !  checks point p against extremities of x()

      IMPLICIT NONE

      INTEGER(kind=4) :: ndime
      REAL (kind=8) :: px(ndime), x(ndime,3)
      REAL (kind=8), OPTIONAL :: toli

      END FUNCTION nonear

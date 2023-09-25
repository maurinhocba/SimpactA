      SUBROUTINE vol2d(dvol,x,xcen,ntype)
!---------------------------------------------------
!***  calculates volume in 2D
!---------------------------------------------------

      IMPLICIT NONE

      INTEGER(kind=4):: ntype
      REAL(kind=8):: dvol,x(2,2),xcen(2)

      REAL(kind=8) :: x1(2),x2(2)
      REAL (kind=8), PARAMETER  :: twopi=6.283185307179586

      X1(:) = X(:,2)-X(:,1)
      X2(:) = XCEN(:)-X(:,1)
      dvol = x1(1)*x2(2) - x1(2)*x2(1)
      dvol = dvol/2.d0
      IF(ntype == 3) dvol=dvol*twopi*(x(1,1)+x(1,2))/3d0

      END SUBROUTINE vol2d

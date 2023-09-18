      SUBROUTINE tetr(dvol,x,xcen)
!---------------------------------------------------
!***  calculates volume of a tetrahedral
!---------------------------------------------------

      IMPLICIT NONE

      REAL(kind=8):: dvol,x(3,3),xcen(3)

      REAL(kind=8) :: x1(3),x2(3),x3(3),xp(3)

      X1(:) = X(:,2)-X(:,1)
      X2(:) = X(:,3)-X(:,1)
      X3(:) = X(:,1)-XCEN(:)
      CALL vecpro(x1,x2,xp)
      dvol = DOT_PRODUCT(x3,xp)
      dvol = dvol/6.d0

      END SUBROUTINE tetr

      SUBROUTINE pyrmd4(dvol,x,xcen)
!---------------------------------------------------
!***  calculates volume of a pyramid with a quad base (1GP scheme)
!---------------------------------------------------

      IMPLICIT NONE

      REAL(kind=8):: dvol,x(3,4),xcen(3)

      REAL(kind=8) :: x1(3),x2(3),x3(3),xp(3),xcq(3)

      x1(:) = -x(:,1)+x(:,2)+x(:,3)-x(:,4)
      x2(:) = -x(:,1)-x(:,2)+x(:,3)+x(:,4)
      xcq(:) = .25d0 * (x(:,1)+x(:,2)+x(:,3)+x(:,4))
      x3(:) = xcq(:)-xcen(:)
      CALL vecpro(x1,x2,xp)
      dvol = DOT_PRODUCT(x3,xp)
      dvol = dvol/12.d0

      END SUBROUTINE pyrmd4

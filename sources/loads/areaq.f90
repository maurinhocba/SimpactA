      SUBROUTINE areaq(x,area,vn)
!---------------------------------------------------
!***  calculates area of a quad base (1GP scheme)
!---------------------------------------------------

      IMPLICIT NONE

      REAL(kind=8):: area,x(3,4),vn(3)

      REAL(kind=8) :: x1(3),x2(3)

      x1(:) = (-x(:,1)+x(:,2)+x(:,3)-x(:,4))/2d0
      x2(:) = (-x(:,1)-x(:,2)+x(:,3)+x(:,4))/2d0
      CALL vecpro(x1,x2,vn)
      CALL vecuni(3,vn,area)

      END SUBROUTINE areaq

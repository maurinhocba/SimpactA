      FUNCTION hmin(x)
!---------------------------------------------------
!***  calculates minimum height of a tetrahedron
!---------------------------------------------------

      IMPLICIT NONE
      REAL(kind=8):: hmin, x(3,4)

      REAL(kind=8):: vol

      REAL(kind=8) :: x1(3),x2(3),x3(3),xp(3),h(4),a

      x1(:) = x(:,2)-x(:,1)
      x2(:) = x(:,3)-x(:,1)
      x3(:) = x(:,4)-x(:,1)
      CALL vecpro(x1,x2,xp)
      vol = DOT_PRODUCT(x3,xp)
      vol = vol/6.d0
      a = SQRT(DOT_PRODUCT(xp,xp)) / 2d0
      h(1) = 3d0 * vol / a

      x1(:) = x(:,3)-x(:,2)
      x2(:) = x(:,4)-x(:,2)
      CALL vecpro(x1,x2,xp)
      a = SQRT(DOT_PRODUCT(xp,xp)) / 2d0
      h(2) = 3d0 * vol / a

      x1(:) = x(:,4)-x(:,1)
      x2(:) = x(:,2)-x(:,1)
      CALL vecpro(x1,x2,xp)
      a = SQRT(DOT_PRODUCT(xp,xp)) / 2d0
      h(3) = 3d0 * vol / a

      x1(:) = x(:,3)-x(:,1)
      x2(:) = x(:,4)-x(:,1)
      CALL vecpro(x1,x2,xp)
      a = SQRT(DOT_PRODUCT(xp,xp)) / 2d0
      h(4) = 3d0 * vol / a

      hmin = MIN (h(1), h(2), h(3), h(4))

      END FUNCTION hmin

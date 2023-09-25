      SUBROUTINE area2d(x,area,vn,ntype)
!---------------------------------------------------
!***  calculates volume of a pyramid with a quad base (1GP scheme)
!---------------------------------------------------

      IMPLICIT NONE

      INTEGER(kind=4):: ntype
      REAL(kind=8):: area,x(2,2),vn(2)

      REAL(kind=8) :: x1(2)
      REAL (kind=8), PARAMETER  :: twopi=6.283185307179586

      X1(:) = X(:,2)-X(:,1)
      vn(1) = -x1(2)
      vn(2) =  x1(1)
      CALL vecuni(2,vn,area)
      IF(ntype == 3) area=area*twopi*(x(1,1)+x(1,2))/2d0

      END SUBROUTINE area2d

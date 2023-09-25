       SUBROUTINE jacob5(cartd,deriv,jac,elcod,nnode,ielem,istop)
!*********************************************************************
!
!**** this SUBROUTINE evaluates the jacobian matrix and the cartesian
!     shape FUNCTION derivatives
!
!*********************************************************************
      USE lispa0, ONLY : lures
      IMPLICIT NONE

      INTEGER (kind=4), PARAMETER :: ndime = 3
      INTEGER (kind=4) nnode,ielem,istop
      REAL    (kind=8) cartd(nnode,ndime),deriv(nnode,ndime),jac,       &
     &                 elcod(ndime,nnode)

      INTEGER (kind=4) inode
      REAL    (kind=8) jaci(ndime,ndime), jacm(ndime,ndime)

!     create jacobian matrix jacm

      jacm = MATMUL(elcod,deriv)

!     calculate determinant and inverse of jacobian matrix

      jac =  jacm(1,1)*jacm(2,2)*jacm(3,3)                              &
     &      +jacm(1,3)*jacm(2,1)*jacm(3,2)                              &
     &      +jacm(3,1)*jacm(1,2)*jacm(2,3)                              &
     &      -jacm(3,1)*jacm(2,2)*jacm(1,3)                              &
     &      -jacm(3,3)*jacm(1,2)*jacm(2,1)                              &
     &      -jacm(1,1)*jacm(2,3)*jacm(3,2)
      IF(jac <= 0) THEN
        WRITE(*,600,ERR=9999) ielem
        WRITE(*,900,ERR=9999)
        WRITE(*,910,ERR=9999) (inode,elcod(1:ndime,inode),inode=1,nnode)
        WRITE(lures,600,ERR=9999) ielem
        WRITE(lures,900,ERR=9999)
        WRITE(lures,910,ERR=9999) (inode,elcod(1:ndime,inode),inode=1,nnode)
        iSTOP = 1
        RETURN
      END IF

      jaci(1,1) =  (jacm(2,2)*jacm(3,3)-jacm(2,3)*jacm(3,2))/jac
      jaci(2,1) = -(jacm(2,1)*jacm(3,3)-jacm(3,1)*jacm(2,3))/jac
      jaci(3,1) =  (jacm(2,1)*jacm(3,2)-jacm(2,2)*jacm(3,1))/jac
      jaci(1,2) = -(jacm(1,2)*jacm(3,3)-jacm(1,3)*jacm(3,2))/jac
      jaci(2,2) =  (jacm(1,1)*jacm(3,3)-jacm(3,1)*jacm(1,3))/jac
      jaci(3,2) = -(jacm(1,1)*jacm(3,2)-jacm(1,2)*jacm(3,1))/jac
      jaci(1,3) =  (jacm(1,2)*jacm(2,3)-jacm(1,3)*jacm(2,2))/jac
      jaci(2,3) = -(jacm(1,1)*jacm(2,3)-jacm(1,3)*jacm(2,1))/jac
      jaci(3,3) =  (jacm(1,1)*jacm(2,2)-jacm(2,1)*jacm(1,2))/jac

!     calculate cartesian derivatives

      cartd = MATMUL(deriv,jaci)

      RETURN
  600 FORMAT(//,' PROGRAM halted in SUBROUTINE jacob5',/,11x,           &
     &' zero or negative volume',/,10x,' element number ',i5)
  900 FORMAT(//,5x,'coordinates of element nodes')
  910 FORMAT(5x,i5,3e15.8)
 9999 CALL runen2('')
      END SUBROUTINE jacob5

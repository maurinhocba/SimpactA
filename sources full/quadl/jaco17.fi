       SUBROUTINE jaco17(cartd,deriv,detj,elcod,elem,nnode,istop,jaci)
!*********************************************************************
!
!**** this SUBROUTINE evaluates the jacobian matrix and the cartesian
!     shape function derivatives
!
!*********************************************************************
      USE lispa0
      IMPLICIT NONE

      INTEGER (kind=4), PARAMETER :: ndime=2
      INTEGER (kind=4) elem,nnode,istop
      REAL    (kind=8) cartd(nnode,ndime),deriv(nnode,ndime),detj,      &
     &                 elcod(ndime,nnode),jaci(ndime,ndime)

      REAL    (kind=8)  jacm(ndime,ndime)

!     ***  create jacobian matrix jacm

      jacm = MATMUL(elcod,deriv) !jacm(i,j) = d[x(i)]/d[xita(j)]

!     ***  calculate determinant and inverse of jacobian matrix

      detj = jacm(1,1)*jacm(2,2) - jacm(1,2)*jacm(2,1)
      IF(detj <= 0) THEN
        WRITE(*,"(//, '  PROGRAM halted in SUBROUTINE jacob3')")
        WRITE(*,"(11x,' zero or negative area')")
        WRITE(*,"(11x,' element number ',i5)") elem
        WRITE(*,"(2e15.4)") elcod,detj
        WRITE(lures,"(//, '  PROGRAM halted in SUBROUTINE jacob3')",ERR=9999)
        WRITE(lures,"(11x,' zero or negative area')",ERR=9999)
        WRITE(lures,"(11x,' element number ',i5)",ERR=9999) elem
        WRITE(lures,"(2e15.4)",ERR=9999) elcod,detj
        iSTOP = 1
        RETURN
      END IF
      jaci(1,1) =  jacm(2,2)/detj       !jaci(i,j) = d[xita(i)]/d[x(j)]
      jaci(1,2) = -jacm(1,2)/detj
      jaci(2,1) = -jacm(2,1)/detj
      jaci(2,2) =  jacm(1,1)/detj

!     ***  calculate cartesian derivatives

      cartd = MATMUL(deriv,jaci)       ! cartd(i,j)= d[N(i)]/d[x(j)]

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE jaco17

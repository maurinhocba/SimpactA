      SUBROUTINE shape1t (elcod,deriv,jacm,detj,elem,nnode,ndime,istop)

      USE lispa0
      IMPLICIT NONE

      INTEGER(kind=4) :: elem, nnode, ndime, istop
      REAL    (kind=8) deriv(nnode,ndime),detj,elcod(ndime,nnode), &
     &                 jacm(ndime,ndime)


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
        WRITE(lures,"(2e15.4)") elcod,detj
        iSTOP = 1
        RETURN
      END IF

      RETURN
 9999 CALL runen2('')
      END SUBROUTINE shape1t

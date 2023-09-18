      SUBROUTINE deriv1 (elcod,shape,deriv,bmat,jacm,detj,rad,ntype,nel,nnode,ndime)

      IMPLICIT NONE

!     **************************************************************
!     EVALUATION OF THE GLOBAL DERIVATIVE OPERATOR (B) AT A POINT (R,S)
!     FOR A QUADRILATERAL ELEMENT HAVING PLANAR OR AXISYMMETRIC GEOMETRY
!     **************************************************************

      INTEGER*4  ntype,nel,nnode,ndime !,K
      REAL (kind=8) :: deriv(nnode,ndime),detj,elcod(ndime,nnode)
      REAL*8  shape(nnode),bmat(ndime,nnode),jacm(ndime,ndime),RAD
      REAL*8  jaci(ndime,ndime)

!     COMPUTE INVERSE OF THE JACOBIAN MATRIX

      jaci(1,1) =  jacm(2,2)/detj       !jaci(i,j) = d[xita(i)]/d[x(j)]
      jaci(1,2) = -jacm(1,2)/detj
      jaci(2,1) = -jacm(2,1)/detj
      jaci(2,2) =  jacm(1,1)/detj

!     EVALUATE GLOBAL DERIVATIVE OPERATOR ( B-MATRIX )

      bmat = TRANSPOSE(MATMUL(deriv,jaci))
!      DO 10 K=1,nnode
!      bmat(1,K) = jaci(1,1)*deriv(K,1) + jaci(1,2)*deriv(K,2)
!   10 bmat(2,K) = jaci(2,1)*deriv(K,1) + jaci(2,2)*deriv(K,2)

      RAD = 1.0D0
      IF(ntype /= 3) RETURN

!     COMPUTE THE RADIUS AT POINT (R,S) FOR AXISYMMETRIC SOLID

      rad = DOT_PRODUCT(elcod(1,:),shape)

      IF(rad >  1.0E-5) RETURN
      WRITE(6,3000) NEL
      STOP

 3000 FORMAT(//50H **ERROR**  ZERO RADIUS ENCOUNTERED IN ELEMENT NO.,I5)

  500 RETURN
      END SUBROUTINE deriv1

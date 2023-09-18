SUBROUTINE LINEQU(X1, Y1, X2, Y2, A, B, C)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE kind_param, ONLY:  DOUBLE
!-----------------------------------------------------------------------------
!
!   CALCULATES COEFFICIENTS A, B AND C OF
!   THE EQUATION OF THE LINE: A*X + B*Y = C
!   PASSING THROUGH THE POINTS (X1,Y1) AND (X2,Y2)
!
!-----------------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  REAL(DOUBLE) , INTENT(IN) :: X1
  REAL(DOUBLE) , INTENT(IN) :: Y1
  REAL(DOUBLE) , INTENT(IN) :: X2
  REAL(DOUBLE) , INTENT(IN) :: Y2
  REAL(DOUBLE) , INTENT(INOUT) :: A
  REAL(DOUBLE) , INTENT(OUT) :: B
  REAL(DOUBLE) , INTENT(OUT) :: C
END SUBROUTINE LINEQU

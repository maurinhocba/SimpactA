 SUBROUTINE eige17(stran,r1,r2,lb)

 ! compute eigen decomposition from Shear log strains

 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT (IN) :: stran(3) !(IN) Shear log strains (twice)
 REAL(kind=8), INTENT (OUT) :: r1,r2, & !components of first eigevector
                               lb(2)    !eigenvalues

 ! local variables
 REAL (kind=8) :: c,d,r

 !compute eigenvalues
 c = (stran(1)+stran(2))/2d0              !center of circle
 d = (stran(1)-stran(2))/2d0              !semi-difference
 r = SQRT(d**2+stran(3)**2)               !circle radius
 lb(1) = c+r                              !first (maximum) eigenvalue
 lb(2) = c-r                              !second (minimum) eigenvalue
 !compute eigenvectors
 IF( r > 0d0 )THEN                        !check
   c = ATAN2(stran(3),d)/2d0              !compute angle
   r1 = COS(c)                            !first  component of eigenvector
   r2 = SIN(c)                            !second component of eigenvector
 ELSE
   r1 = 1d0                               !local direction is the vector
   r2 = 0d0
 END IF
 RETURN
 END SUBROUTINE eige17

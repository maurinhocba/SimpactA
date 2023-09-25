 SUBROUTINE eige18(str,r,lb,flag)

 ! compute eigenvalues and eigenvector of symmetric tensor STR

 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), INTENT (IN) :: str(6)     !Symmetric tensor
 REAL(kind=8), INTENT (OUT) :: r(3,3), & !components of eigevectors
                               lb(3)     !eigenvalues of STR
 INTEGER (kind=4), INTENT(OUT) :: flag   !flag

 ! local variables
 INTEGER (kind=4) k
 REAL (kind=8) i1,i2,i3,i13,phi,la,a11,a12,a13,a22,a23,a33,d11,d22,d33
 REAL (kind=8),PARAMETER :: s23  = 1.15470053837925, pi6  = 0.523598775598299,  &
                            pi23 = 2.09439510239320, pi43 = 4.18879020478639,   &
                            tol = 1d-16

 !
 !     compute principal values
 !

 i13  = (str(1) + str(2) + str(3))/3d0              !I1/3
 !  deviator
 d11 = str(1) - i13
 d22 = str(2) - i13
 d33 = str(3) - i13
 a12 = str(4)
 a13 = str(5)
 a23 = str(6)
 !  check if the tensor is not null
 !i2  = s23*SQRT(-d11*d22 - d11*d33 - d22*d33 + a12**2 + a13**2 + a23**2)  !2 SQRT(-I2/3)
 i2 = -d11*d22 - d11*d33 - d22*d33 + a12**2 + a13**2 + a23**2
 
 IF(i2 > tol )THEN     !if radius is large enough
   i2 = s23*SQRT(i2)
   i3  = d11*d22*d33 + 2d0*a12*a13*a23 -d11*a23**2 -d22*a13**2 -d33*a12**2      !i3
   phi = -4d0*i3/(i2**3)       !sino
   IF( phi > 1d0 )THEN         !check range for round errors
     phi = pi6
   ELSE IF (phi < -1d0) THEN
     phi = -pi6
   ELSE                        !normally
     phi = ASIN(phi)/3d0
   END IF
   lb(1) = i2*SIN(phi+pi23)         !largest eigenvalue
   lb(2) = i2*SIN(phi)              !mediun eigenvalue
   lb(3) = i2*SIN(phi+pi43)         !smallest eigenvalue

   !     compute principal directions

   DO k=1,3,2       !for largest and smallest eigenvalues
     !  diagonal terms
     a11 = d11 - lb(k)
     a22 = d22 - lb(k)
     a33 = d33 - lb(k)
     !  determinants
     i1 =  a22*a33 - a23**2
     i2 =  a11*a33 - a13**2
     i3 =  a11*a22 - a12**2
     !  choose the system with the largest determinant
     IF( ABS(i1) >= ABS(i2) )THEN      !compare first and second
       IF( ABS(i1) >= ABS(i3) )THEN    !compare first and third
         ! first is the largest, use 2nd and 3rd rows
         i2 = ( -a12*a33 + a13*a23 )/i1
         i3 = ( -a22*a13 + a23*a12 )/i1
         la = SQRT(1d0+i2**2+i3**2)
         r(1,k) = 1d0/la
         r(2,k) =  i2/la
         r(3,k) =  i3/la
       ELSE
         ! third is the largest, use 1st and 2nd rows
         i1 = ( -a13*a22 + a23*a12 )/i3
         i2 = ( -a11*a23 + a12*a13 )/i3
         la = SQRT(i1**2+i2**2+1d0)
         r(1,k) =  i1/la
         r(2,k) =  i2/la
         r(3,k) = 1d0/la
       END IF
     ELSE IF( ABS(i2) >= ABS(i3) )THEN
       ! second is the largest, use 1st and 3rd rows
       i1 = ( -a12*a33 + a23*a13 )/i2
       i3 = ( -a11*a23 + a13*a12 )/i2
       la = SQRT(i1**2+1d0+i3**2)
       r(1,k) =  i1/la
       r(2,k) = 1d0/la
       r(3,k) =  i3/la
     ELSE
       ! third is the largest, use 1st and 2nd rows
       i1 = ( -a13*a22 + a23*a12 )/i3
       i2 = ( -a11*a23 + a12*a13 )/i3
       la = SQRT(i1**2+i2**2+1d0)
       r(1,k) =  i1/la
       r(2,k) =  i2/la
       r(3,k) = 1d0/la
     END IF
   END DO
   !  Compute second eigenvector r2 = r3 x r1
   r(1,2) = r(2,3)*r(3,1) - r(2,1)*r(3,3)
   r(2,2) = r(3,3)*r(1,1) - r(3,1)*r(1,3)
   r(3,2) = r(1,3)*r(2,1) - r(1,1)*r(2,3)
   lb = lb + i13           !add mean value
   flag = 0                !eigen-decomposition O.K.
 ELSE
   flag = 1                !no eigen-decomposition performed
   r = RESHAPE((/ 1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0 /),(/3,3/))
   lb = str(1:3)
 END IF
 RETURN
 END SUBROUTINE eige18

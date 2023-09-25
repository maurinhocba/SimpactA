 SUBROUTINE bmem13(a,b,cd,a0,gamma,bmem,t,mems)

 !*** membrane matrix  Linear Strain Triangle
 IMPLICIT NONE

 LOGICAL, INTENT(IN) :: mems(3)
 REAL (kind=8), INTENT(IN) :: a(3),b(3),cd(4,2,3),t(3,15),a0(3),gamma(3)
 REAL (kind=8), INTENT(OUT) :: bmem(3,6,3)

 INTEGER (kind=4) :: i,j,n,n1
 REAL(kind=8) :: d1,d2

 bmem = 0d0        !initializes
 n = 10            !initialize pointer
 DO i=1,3          !for each side
   IF( mems(i) )THEN  !linear strain approach
     n1 = n+1              !position of second derivative
     d1 = ABS(a0(i)+2d0*gamma(i))    !angle at boundary
     IF( d1 <= alp1 )THEN            !smooth  use side (cuadratic) metric tensor
        DO j=1,3          !node
          bmem(:,j,1) = bmem(:,j,1) + cd(j,1,i)*t(:,n)
          bmem(:,j,2) = bmem(:,j,2) + cd(j,2,i)*t(:,n1)
          bmem(:,j,3) = bmem(:,j,3) + cd(j,2,i)*t(:,n) + cd(j,1,i)*t(:,n1)
        END DO
        j = i+3          ! extra node
        bmem(:,j,1) = bmem(:,j,1) + cd(4,1,i)*t(:,n)
        bmem(:,j,2) = bmem(:,j,2) + cd(4,2,i)*t(:,n1)
        bmem(:,j,3) = bmem(:,j,3) + cd(4,2,i)*t(:,n) + cd(4,1,i)*t(:,n1)

     ELSE IF( d1 >= alp2 )THEN       !folded, use main triangle metric tensor
       DO j=1,3          !node
         bmem(:,j,1) = bmem(:,j,1) - b(j)*t(:,1)
         bmem(:,j,2) = bmem(:,j,2) + a(j)*t(:,2)
         bmem(:,j,3) = bmem(:,j,3) + a(j)*t(:,1) - b(j)*t(:,2)
       END DO

     ELSE                            !interpolate
       d1 = (d1-alp1)/alp3           !function for main triangle
       d2 = 1d0 - d1                 !function for cuadratic interpolation
       DO j=1,3          !node
         bmem(:,j,1) = bmem(:,j,1)- d1*b(j)*t(:,1) + d2*cd(j,1,i)*t(:,n)
         bmem(:,j,2) = bmem(:,j,2)+ d1*a(j)*t(:,2) + d2*cd(j,2,i)*t(:,n1)
         bmem(:,j,3) = bmem(:,j,3)+ d1*(a(j)*t(:,1) - b(j)*t(:,2)) &
                                  + d2*(cd(j,2,i)*t(:,n) + cd(j,1,i)*t(:,n1))
       END DO
       j = i+3          ! extra node
       bmem(:,j,1) = bmem(:,j,1) + d2*cd(4,1,i)*t(:,n)
       bmem(:,j,2) = bmem(:,j,2) + d2*cd(4,2,i)*t(:,n1)
       bmem(:,j,3) = bmem(:,j,3) + d2*(cd(4,2,i)*t(:,n) + cd(4,1,i)*t(:,n1))

     END IF
   ELSE   ! constant strain triangle
     DO j=1,3          !node
       bmem(:,j,1) = bmem(:,j,1) - b(j)*t(:,1)
       bmem(:,j,2) = bmem(:,j,2) + a(j)*t(:,2)
       bmem(:,j,3) = bmem(:,j,3) + a(j)*t(:,1) - b(j)*t(:,2)
     END DO
   END IF
   n = n+2
 END DO
 RETURN
 END SUBROUTINE bmem13

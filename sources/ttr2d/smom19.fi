 SUBROUTINE smom19(Area1,emass)
 !=======================================================================
 ! This routine evaluates the element smoothing matrix and inverts it
 ! for element 19
 !=======================================================================
 IMPLICIT NONE

   !--- Dummy variables
   REAL(kind=8),INTENT(IN):: Area1
   REAL(kind=8),INTENT(OUT):: emass(3,3)
   !--- Local variables
   INTEGER(kind=4):: i, j
   REAL(kind=8):: deter, denom, a(3,3)

   !Consistent smoothing matrix
   a(1,1) = Area1/6d0
   a(1,2) = Area1/12d0
   a(1,3) = Area1/12d0
   a(2,2) = Area1/6d0
   a(2,3) = Area1/12d0
   a(3,3) = Area1/6d0
   DO i=1,3
     DO j=1,i-1
       a(i,j) = a(j,i)
     END DO
   END DO


   deter = a(1,1)*a(2,2)*a(3,3) + a(2,1)*a(3,2)*a(1,3) + a(3,1)*a(2,1)*a(2,3)   &
          -a(1,3)*a(2,2)*a(3,1) - a(2,3)*a(3,2)*a(1,1) - a(3,3)*a(1,2)*a(2,1)
   denom = 1d0/deter

   emass(1,1) = denom*( a(2,2)*a(3,3) - a(3,2)*a(2,3) )
   emass(2,1) = denom*( a(3,1)*a(2,3) - a(2,1)*a(3,3) )
   emass(3,1) = denom*( a(2,1)*a(3,2) - a(3,1)*a(2,2) )
   emass(1,2) = denom*( a(3,2)*a(1,3) - a(1,2)*a(3,3) )
   emass(2,2) = denom*( a(1,1)*a(3,3) - a(3,1)*a(1,3) )
   emass(3,2) = denom*( a(3,1)*a(1,2) - a(1,1)*a(3,2) )
   emass(1,3) = denom*( a(1,2)*a(2,3) - a(2,2)*a(1,3) )
   emass(2,3) = denom*( a(2,1)*a(1,3) - a(1,1)*a(2,3) )
   emass(3,3) = denom*( a(1,1)*a(2,2) - a(2,1)*a(1,2) )

 RETURN
 END SUBROUTINE smom19

 SUBROUTINE reduce_mass(mass,maxav,neq,maxa)
 !
 !  reduce mass matrix size checking null values
 !
 IMPLICIT NONE
 INTEGER(kind=4), INTENT(IN) :: neq
 INTEGER(kind=4), INTENT(IN OUT) :: maxav(:)
 INTEGER(kind=4), INTENT(OUT) :: maxa
 REAL(kind=8), INTENT(IN OUT) :: mass(:)

 INTEGER(kind=4) :: i,j,j1,k,h,l
 REAL(kind=8)  :: diag
 REAL(kind=8), PARAMETER  :: tol=1e-3

 i=2              !second column
 j  = maxav(i)    !original diagonal position
 j1 = maxav(i+1)  !original diagonal position  (next eq.)
 DO   ! loop
   diag = mass(j1)
   k = maxav(i)    !new  diagonal position
   h = 0          !initializes column height
   DO l=j1-1,j+1,-1              !for each value in column
     IF( ABS(mass(l))/diag > tol )THEN       !find first non zero value
       h = l - j                    !height
       EXIT
     END IF
   END DO
   mass(k:k+h) = mass(j:j+h)     !pass new values
   maxav(i+1)  = maxav(i) + h + 1 !new position of next diagonal value
   j = j1                        !keep old position of next diagonal value
   IF( i == neq ) EXIT
   i = i+1                       !update column
   j1 = maxav(i+1)                !keep old position of next column
 END DO
 PRINT *, maxa, maxav(neq+1)
 maxa = maxav(neq+1)
 RETURN
 END SUBROUTINE reduce_mass

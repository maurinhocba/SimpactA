 SUBROUTINE toar04n(npoin)

 ! Generate nodal connectivities and space for B matrices

 IMPLICIT NONE
 INTEGER (kind=4) :: npoin

 ! local variables
 INTEGER (kind=4) i,j

 !_____________________________


 DO i=1,npoin
   j = 1
   DO
     IF( ne(j,i) == 0 )EXIT
     j = j+1
   END DO
   ncc(1,i) = j              !size
   IF( j > 1 )THEN
     ALLOCATE( nlnod(i)%p(j) )
     nlnod(i)%p(1) = i              !the node
     nlnod(i)%p(2:j) = ne(1:j-1,i)  !rest of nodes
     ALLOCATE( ncard(i)%p(j,3) )    !B matrix
     ncard(i)%p = 0d0               !initializes
   END IF
 END DO

 RETURN
 END SUBROUTINE toar04n

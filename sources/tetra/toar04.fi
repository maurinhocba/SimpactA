 SUBROUTINE toar04(heade,lside)

 ! Create topological array

 IMPLICIT NONE
 TYPE (ele04), POINTER :: heade   !pointer to first element in the set
 LOGICAL :: lside

 ! local variables
 INTEGER (kind=4) i,j,k,l,k1,k2
 TYPE (ele04), POINTER  :: e                            !element pointers

 !_____________________________

 e => heade             !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT !if all element processed EXIT loop
   DO i=1,4                        !for each node
     j = e%lnods(i)            !node number
     IF( ncc(2,j) == 0 ) ncc(2,j) = e%matno    !Associated section
     DO k=i+1,4
       l = e%lnods(k)            !node number
       DO k1=1,29
         IF( ne(k1,j) == l )EXIT        !already
         IF( ne(k1,j) == 0 )THEN        !not yed
           ne(k1,j) = l                 !add to list
           DO k2=1,29
             IF( ne(k2,l) == 0 )THEN
               ne(k2,l) = j               !add to list
               EXIT
             END IF
           END DO
           EXIT
         END IF
       END DO
     END DO
   END DO
   e => e%next               !point to next element to process
 END DO

 lside = .TRUE.              !topological array created

 RETURN
 END SUBROUTINE toar04

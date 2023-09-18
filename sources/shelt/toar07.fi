 SUBROUTINE toar07(heade,lside,nelem)

 ! Create topological array
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem !number of elements
 TYPE (ele07), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT(OUT) :: lside

 ! local variables
 INTEGER (kind=4) :: i,ii,iii,j,jj
 TYPE (ele07), POINTER  :: e, ej               !element pointers

 !_____________________________
 e => heade             !point to first element
 DO i=1,nelem-1         !for each element
   DO ii=1,3              !for each element side
     IF(e%lnods(ii+6) > 0 )CYCLE  !extra node already paired
     iii = e%lnods(ii+3)  ! side node
     ! search element with same side node
     ej => e%next         !point to next element
     search : DO j=i+1,nelem
       DO jj=1,3
         IF( ej%lnods(jj+3) == iii )THEN
           e%lnods(ii+6) = ej%lnods(jj)
           ej%lnods(jj+6) = e%lnods(ii)
           EXIT search
         END IF
       END DO
       ej => ej%next
     END DO search
   END DO
   e => e%next
  END DO
 lside = .TRUE.              !topological array created
 RETURN
 END SUBROUTINE toar07

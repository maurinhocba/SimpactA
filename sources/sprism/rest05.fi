 SUBROUTINE rest05 (nelem, nnode, nreqs, heade, taile, ngrqs, ngaus)
 !******************************************************************
 !
 !***  restore an element set at restart for shell element CST-BST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs, & !Number of elements for requested output
                                 nnode, & !number of nodes per elements in the set
                                 ngaus, & !number of integration points
                                 nelem    !number of elements in the set
 INTEGER (kind=4), POINTER  :: ngrqs(:) !elements for requested output
 TYPE (ele05), POINTER  :: heade,taile  !pointers to first and last element
 ! local variables
 INTEGER (kind=4) :: i,j,ielem,size1,size2
 TYPE (ele05), POINTER :: e

 DO ielem=1,nelem    !for each element

   CALL new_ele05e (e)    !reserve memory for data
   ALLOCATE (e%lnods(nnode),e%stint(6,ngaus))  !get memory
   READ (51) e%numel,e%matno,e%lnods,size1,size2
   READ (51) e%angle,e%dvol,e%cartd,e%nfdas,e%jacin,e%alpha
   IF( nnode == 12 )THEN
     ALLOCATE( e%cdq(4,2,3,2) )
   ELSE
     ALLOCATE( e%cdq(3,2,1,2) )
   END IF
   READ (51) e%cdq
   READ (51) e%stint
   IF( size1 > 0 )THEN
     ALLOCATE (e%gausv(size1,size2))             !get space for internal vars
     READ (51) ((e%gausv(i,j),i=1,size1),j=1,size2) !read internal variables
   END IF

   CALL add_ele05e (e, heade, taile)           !add element to the list

 END DO
 IF ( nreqs > 0 )THEN
   ALLOCATE(ngrqs(nreqs))
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF

 RETURN
 END SUBROUTINE rest05

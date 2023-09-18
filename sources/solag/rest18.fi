 SUBROUTINE rest18 (nelem, nreqs, heade, taile, ngrqs, ngaus,shell)
 !******************************************************************
 !
 !***  restore an element set at restart for shell element CST-BST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs, & !Number of elements for requested output
                                 ngaus, & !number of integration points
                                 nelem    !number of elements in the set
 INTEGER (kind=4), POINTER  :: ngrqs(:) !elements for requested output
 TYPE (ele18), POINTER  :: heade,taile  !pointers to first and last element
 LOGICAL ::  shell         !.TRUE. assumed strain for transverse shear

 INTEGER (kind=4) :: i,j,ielem,size1,size2
 TYPE (ele18), POINTER :: e

 DO ielem=1,nelem    !for each element

   CALL new_ele18e (e)      !get memory for an element
   ALLOCATE (e%cartd(nnode,3,ngaus),e%dvol(ngaus),e%stint(6,ngaus))  !get memory
   READ (51) e%numel,e%matno,e%lnods,size1,size2
   READ (51) e%dvol,e%angle,e%cartd,e%stint
   IF( size1 > 0 )THEN
     ALLOCATE (e%gausv(size1,size2))             !get space for internal vars
     READ (51) ((e%gausv(i,j),i=1,size1),j=1,size2) !read internal variables
   END IF
   IF( shell )THEN
     ALLOCATE(e%nfdas(8,4),e%jacin(2,2,ngaus))
     READ (51) e%nfdas,e%jacin                   !read values
   END IF
   CALL add_ele18e (e, heade, taile)           !add element to the list
 END DO
 IF ( nreqs > 0 )THEN
   ALLOCATE(ngrqs(nreqs))
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF

 RETURN
 END SUBROUTINE rest18

 SUBROUTINE rest04 (nelem, nreqs, heade, taile, ngrqs)
 !******************************************************************
 !
 !***  restore an element set at restart for shell element CST-BST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs,  & !Number of elements for requested output
                                 nelem     !number of elements in the set

 INTEGER (kind=4), POINTER  :: ngrqs(:) !elements for requested output
 TYPE (ele04), POINTER  :: heade,taile  !pointers to first and last element

 INTEGER (kind=4) :: i,ielem,size1
 TYPE (ele04), POINTER :: e


 DO ielem=1,nelem    !for each element

   ALLOCATE (e)      !get memory for an element
   READ (51) e%numel,e%matno,e%lnods,size1
   READ (51) e%dvol,e%angle,e%cd,e%stint
   IF( size1 > 0 )THEN
     ALLOCATE (e%gausv(size1))             !get space for internal vars
     READ (51) (e%gausv(i),i=1,size1)      !read internal variables
   ELSE
     NULLIFY( e%gausv )
   END IF

   CALL add_ele04e (e, heade, taile)           !add element to the list

 END DO
 IF ( nreqs > 0 )THEN
   ALLOCATE(ngrqs(nreqs))
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF

 RETURN
 END SUBROUTINE rest04

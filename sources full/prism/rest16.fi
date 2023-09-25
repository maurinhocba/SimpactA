 SUBROUTINE rest16 (nelem, nreqs, heade, taile, ngrqs, ngaus, quad, shell)
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
 TYPE (ele16), POINTER  :: heade,taile  !pointers to first and last element
 LOGICAL, INTENT(IN) :: quad,shell

 ! local variables
 INTEGER (kind=4) :: i,j,ielem,size1,size2,nn
 TYPE (ele16), POINTER :: e

 IF( quad )THEN
   nn = 12
 ELSE
  nn = nnode
 END IF

 DO ielem=1,nelem    !for each element

   CALL new_ele16e (e)    !reserve memory for data
   ALLOCATE (e%lnods(nn),e%cartd(nnode,3,ngaus),e%dvol(ngaus),e%stint(6,ngaus))  !get memory
   READ (51) e%numel,e%matno,e%lnods,size1,size2
   READ (51) e%dvol,e%angle,e%cartd,e%stint
   IF( quad )THEN
     ALLOCATE( e%cdq(4,2,3,2) )
     READ (51) e%cdq
   END IF
   IF( shell )THEN
     ALLOCATE(e%nfdas(6,3,ngaus),e%jacin(2,2,ngaus))
     READ (51) e%nfdas,e%jacin
     e%angle(2:3) = 0d0
   END IF
   IF( size1 > 0 )THEN
     ALLOCATE (e%gausv(size1,size2))             !get space for internal vars
     READ (51) ((e%gausv(i,j),i=1,size1),j=1,size2) !read internal variables
   END IF

   CALL add_ele16e (e, heade, taile)           !add element to the list

 END DO
 IF ( nreqs > 0 )THEN
   ALLOCATE(ngrqs(nreqs))
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF

 RETURN
 END SUBROUTINE rest16

 SUBROUTINE rest17 (nelem, nreqs, heade, taile, ngrqs, ngaus)
 !******************************************************************
 !
 !***  restore an element set at restart for shell element CST-BST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs, & !Number of elements for requested output
                                 nelem, & !number of elements in the set
                                 ngaus    !number of gauss points
 INTEGER (kind=4), POINTER  :: ngrqs(:) !elements for requested output
 TYPE (ele17), POINTER  :: heade,taile  !pointers to first and last element

 INTEGER (kind=4) :: i,j,nvare,ielem,nucom,osec
 TYPE (ele17), POINTER :: e
 TYPE (section), POINTER :: sec
 LOGICAL :: plast

 osec = -1
 nvare = 5           !number of internal variables per component
 DO ielem=1,nelem    !for each element

   ALLOCATE (e)      !get memory for an element
   READ (51) e%numel,e%matno,e%lnods
   ALLOCATE (e%dvol(ngaus),e%stint(4,ngaus))
   ALLOCATE (e%cartd(nnode,2,MAX(ngaus,2)))  !reserve space for variables
   READ (51) e%dvol,e%angle,e%cartd,e%stint
   IF( e%matno /= osec )THEN
     sec => psecs(e%matno)%p
     nucom = sec%iprop(2)              !number of material compon.
     plast = sec%mtbas%matdef(3) > 1   !plastic material
     osec = e%matno
   END IF
   IF( plast )THEN
     ALLOCATE (e%gausv(nvare*nucom,ngaus))             !get space for internal vars
     READ (51) ((e%gausv(i,j),i=1,nvare*nucom),j=1,ngaus) !read internal variables
   ELSE
     ALLOCATE ( e%gausv(1,1) )
   END IF

   CALL add_ele17e (e, heade, taile)           !add element to the list

 END DO
 IF ( nreqs > 0 )THEN
   ALLOCATE(ngrqs(nreqs))
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF
 RETURN
 END SUBROUTINE rest17

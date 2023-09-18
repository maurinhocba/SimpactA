 SUBROUTINE rest14 (nelem, nreqs, heade, taile, ngrqs, stint)
 !******************************************************************
 !
 !***  restore an element set at restart for shell element CST-BST (TLF)
 !
 !******************************************************************
 IMPLICIT NONE
 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nreqs, & !Number of elements for requested output
                                 nelem    !number of elements in the set
 INTEGER (kind=4), POINTER  :: ngrqs(:)   !elements for requested output
 REAL(kind=8), POINTER :: stint(:,:)         !forces, moments and shear
 TYPE (ele14), POINTER  :: heade,taile    !pointers to first and last element

 INTEGER (kind=4) :: i,j,nlayr,nvare,ielem
 TYPE (ele14), POINTER :: e

 DO ielem=1,nelem    !for each element

   CALL new_ele14e(e)    !get memory for an element
   READ (51) e%numel,e%matno,e%lnods,e%lside,nvare,nlayr
   READ (51) e%area1,e%lb,e%angle,e%a,e%b,e%ci,e%stra0,e%stra1,stint(:,ielem)
   IF( nvare > 0 )THEN
     ALLOCATE (e%gausv(nvare,nlayr))             !get space for internal vars
     READ (51) ((e%gausv(i,j),i=1,nvare),j=1,nlayr) !read internal variables
   END IF
   CALL add_ele14e (e, heade, taile)           !add element to the list
 END DO
 IF ( nreqs > 0 ) THEN
   ALLOCATE( ngrqs(nreqs) )
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF
 RETURN
 END SUBROUTINE rest14

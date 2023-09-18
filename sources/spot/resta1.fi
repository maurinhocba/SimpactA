 SUBROUTINE resta1 (nnode, nelem, nreqs, heade, taile, ngrqs,nvare)
 !******************************************************************
 !
 !*** restores an element set at restart
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) ::  nnode, nelem, nreqs,nvare
 INTEGER (kind=4), POINTER :: ngrqs(:)
 TYPE (ele01), POINTER :: heade,taile

 INTEGER (kind=4) :: i,ielem
 TYPE (ele01), POINTER :: e

 DO ielem=1,nelem
   ALLOCATE (e)      !get memory for an element
   ALLOCATE (e%lnods(nnode))      !get memory for an element
   READ (51) e%lnods , e%numel, e%matno, e%rot
   ALLOCATE (e%gausv(nvare))
   READ (51) e%gausv
   CALL add_ele01e (e, heade, taile)           !add element to the list
 END DO
 IF ( nreqs > 0 ) THEN
   ALLOCATE( ngrqs(nreqs) )
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF

 RETURN
 END SUBROUTINE resta1

 SUBROUTINE resta2 ( nelem, nreqs, heade, taile, ngrqs)
 !******************************************************************
 !
 !*** restores an element set at restart
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) ::  nelem, nreqs
 INTEGER (kind=4), POINTER :: ngrqs(:)
 TYPE (ele02), POINTER :: heade,taile

 INTEGER (kind=4) :: i,ielem,isec,osec
 TYPE (ele02), POINTER :: e
 TYPE (section), POINTER :: sec
 LOGICAL :: plast

 DO ielem=1,nelem
   ALLOCATE (e)      !get memory for an element
   READ (51) e%lnods,e%numel,e%matno
   isec = e%matno
   IF( isec /= osec )THEN
     sec => psecs(isec)%p
     plast = sec%mtbas%matdef(3) > 1
     osec = isec
   END IF
   IF( plast )THEN   !reserve space for Gauss variables
     ALLOCATE (e%gausv(8))
   ELSE
     ALLOCATE (e%gausv(4))
   END IF
   READ (51) e%gausv
   CALL add_ele02e (e, heade, taile)           !add element to the list
 END DO
 IF ( nreqs > 0 ) THEN
   ALLOCATE( ngrqs(nreqs) )
   READ (51) (ngrqs(i), i=1,nreqs)  !read elements for output
 END IF

 RETURN
 END SUBROUTINE resta2

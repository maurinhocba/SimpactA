 SUBROUTINE slnods17(heade,lnods,nelem,reorder)
 !***********************************************************************
 !*** Extracts conectivities into the nodal set
 !***********************************************************************
 USE npo_db, ONLY : coord
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4), POINTER :: lnods(:,:)   !Conectivities of the set
   INTEGER(kind=4),INTENT(IN):: nelem       !number of elements
   TYPE(ele17),POINTER:: heade      !pointer to first element
   LOGICAL, INTENT(IN) :: reorder

   !--- Local Variables
   INTEGER(kind=4):: ielem,nn,i, np(1)
   REAL (kind=8) :: l(2),x(2,3),ls(2)
   INTEGER (kind=4), PARAMETER :: la(2,2) = (/ 1,2, 1,3 /)
   TYPE(ele17),POINTER:: elm


   IF( ASSOCIATED(lnods) )THEN
     nn    = SIZE(lnods,1)
     ielem = SIZE(lnods,2)
     IF( ielem /= nelem .OR. nn /= nnode)THEN
       DEALLOCATE(lnods)
       ALLOCATE(lnods(nnode,nelem))
     END IF
   ELSE
     ALLOCATE(lnods(nnode,nelem))
   END IF
   ielem = 0                             !Initiates counter
   elm => heade                          !Point to first element for each element in the set
   DO
     IF (.NOT.ASSOCIATED(elm)) EXIT      !Exit from the loop
     ielem = ielem + 1
     IF( reorder )THEN
       ! compute shortest direction
       x = coord(:,elm%lnods((/ 1,2,4 /)))
       DO i=1,2
         ls = x(:,la(2,i)) - x(:,la(1,i))
         l(i) = DOT_PRODUCT(ls,ls)
       END DO
       np = MINLOC(l)              !shortest direction
       SELECT CASE(np(1))          !set first and last faces normal to shortest direction
       CASE(1)
         lnods(:,ielem) = elm%lnods((/2,3,4,1 /))
       CASE(2)
         lnods(:,ielem) = elm%lnods
       END SELECT
       elm%lnods = lnods(:,ielem)
     ELSE
       lnods(1:nnode,ielem) = elm%lnods   !transfer data
     END IF

     elm => elm%next                     !point to next element
   END DO

 RETURN
 END SUBROUTINE slnods17

 SUBROUTINE slnods20(heade,lnods,nelem,reorder)
 !***********************************************************************
 !*** Extracts conectivities into the nodal set
 !***********************************************************************
 USE npo_db, ONLY : coord
 IMPLICIT NONE
   !  nnode = 3
   !--- Dummy variables
   INTEGER(kind=4), POINTER :: lnods(:,:)   !Conectivities of the set
   INTEGER(kind=4),INTENT(IN):: nelem       !number of elements
   TYPE(ele20),POINTER:: heade      !pointer to first element
   LOGICAL, INTENT(IN) :: reorder

   !--- Local Variables
   INTEGER(kind=4):: ielem,nn,i, np(1)
   REAL (kind=8) :: l(3),x(2,3),ls(2)
   INTEGER (kind=4), PARAMETER :: la(2,3) = (/ 1,2, 2,3, 3,1 /)
   TYPE(ele20),POINTER:: elm


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
       x = coord(:,elm%lnods(1:3))
       DO i=1,3
         ls = x(:,la(2,i)) - x(:,la(1,i))
         l(i) = DOT_PRODUCT(ls,ls)
       END DO
       np = MINLOC(l)              !shortest direction
       SELECT CASE(np(1))          !set first and last faces normal to shortest direction
       CASE(1)
         lnods(1:3,ielem) = elm%lnods(1:3)
       CASE(2)
         lnods(1:3,ielem) = elm%lnods((/2,3,1 /))
       CASE(3)
         lnods(1:3,ielem) = elm%lnods((/3,1,2 /))
       END SELECT
       elm%lnods(1:3) = lnods(:,ielem)
     ELSE
       lnods(1:nnode,ielem) = elm%lnods(1:nnode)   !transfer data
     END IF

     elm => elm%next                     !point to next element
   END DO

 RETURN
 END SUBROUTINE slnods20

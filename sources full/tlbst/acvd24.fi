 SUBROUTINE acvd24(ifpre,heade,lside,nelem,nbs,bhead)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN OUT) :: lside !flag TRUE if topology already computed
 INTEGER (kind=4):: ifpre(:,:), & !(OUT) active DOFs
                    nelem         !(IN) number of elements
 INTEGER (kind=4):: nbs           !(OUT) number of branching sides
 TYPE (ele24), POINTER :: heade   !pointer to first element in the set
 TYPE (sideb), POINTER :: bhead   !pointer to first branching side

 ! local variables
 INTEGER (kind=4) inode,i
 TYPE (ele24), POINTER  :: e                     !element pointers


 !***  create topological array  lnods(4:6,nelem)

   IF(.NOT.lside) CALL toar24(heade,lside,nbs,bhead,nelem)

 !
 !*** activates degrees of freedom
 !
   e => heade                         !point to first element
   DO
     IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
     DO i=1,3
       IF( ifpre(1,e%lnods(i)) == 1 ) &   !to avoid to overwrite RIGID constraints
           ifpre(1:3,e%lnods(i)) = 0      !release DOFs
     END DO
     IF( therm )THEN
       DO inode = 1,3                      !for each node in the element
         iftmp(1:ndoft,e%lnods(inode)) = 0 !release thermal DOF
       END DO
     END IF
     e => e%next                      !point to next element
   END DO

 RETURN
 END SUBROUTINE acvd24

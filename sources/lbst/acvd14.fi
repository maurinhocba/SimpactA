 SUBROUTINE acvd14(ifpre,heade,nelem,lside)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN OUT) :: lside !flag TRUE if topology already computed
 INTEGER (kind=4):: ifpre(:,:)    !(OUT) active DOFs
 INTEGER (kind=4):: nelem         !number of elements
 TYPE (ele14), POINTER :: heade   !pointer to first element in the set

 ! local variables
 INTEGER (kind=4) inode
 TYPE (ele14), POINTER  :: e                     !element pointers

 !***  create topological array  lnods(4:6,nelem)

 IF (.NOT.lside) CALL toar14(lside,nelem,heade)
 !
 !*** activates degrees of freedom
 !
 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
   DO inode = 1,3                   !for each node in the element
     IF( ifpre(1,e%lnods(inode)) == 1 ) &   !to avoid to overwrite RIGID constraints
         ifpre(1:3,e%lnods(inode)) = 0      !release DOFs
   END DO
   e => e%next                      !point to next element
 END DO

 RETURN
 END SUBROUTINE acvd14

 SUBROUTINE acvd04(ifpre,heade,lside,nodvo)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN) :: nodvo     !flag TRUE if nodal volume approach
 LOGICAL, INTENT(IN OUT) :: lside !flag TRUE if topology already computed
 INTEGER (kind=4):: ifpre(:,:)    !(IN OUT)
 TYPE (ele04), POINTER :: heade   !pointer to first element in the set

 ! local variables
 INTEGER (kind=4) inode
 TYPE (ele04), POINTER  :: e                     !element pointers

 !***  create topological array  lnods(5:8,nelem)

 IF(nodvo .AND. .NOT.lside) CALL toar04(heade,lside)

 !*** activates degrees of freedom

 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
   DO inode = 1,4                   !for each node in the element
     IF (ifpre(1,e%lnods(inode)) == 1 ) ifpre(1:3,e%lnods(inode)) = 0 !release DOFs
   END DO
   e => e%next                      !point to next element
 END DO

 RETURN
 END SUBROUTINE acvd04

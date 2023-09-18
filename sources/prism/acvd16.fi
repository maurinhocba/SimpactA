 SUBROUTINE acvd16(ifpre,heade,quad,lside,nelem)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: nelem   !(IN) number of element
 INTEGER (kind=4):: ifpre(:,:)    !(OUT) active DOFs
 TYPE (ele16), POINTER :: heade   !pointer to first element in the set
 LOGICAL :: quad,lside

 ! local variables
 INTEGER (kind=4) inode
 TYPE (ele16), POINTER  :: e                     !element pointers

 !***  create topological array  lnods(7:12,nelem)

 IF(quad )THEN       !for quadratic approach
   IF(.NOT.lside) CALL toar16(lside,nelem,heade)  !generate extended connectivities
 END IF

 !*** activates degrees of freedom
 !
 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
   DO inode = 1,nnode               !for each node in the element
     IF (ifpre(1,e%lnods(inode)) == 1 ) ifpre(1:3,e%lnods(inode)) = 0 !release DOFs
   END DO
   e => e%next                      !point to next element
 END DO

 RETURN
 END SUBROUTINE acvd16

 SUBROUTINE acvd17(ifpre,heade)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 USE esets_db, ONLY : gnods, gelem !INTENT(OUT) gnods(maxnn,melen)
 USE ctrl_db, ONLY : lumped        !INTENT(IN)
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: ifpre(:,:)    !(OUT) active DOFs
 TYPE (ele17), POINTER :: heade   !pointer to first element in the set

 ! local variables
 INTEGER (kind=4) inode
 TYPE (ele17), POINTER  :: e                     !element pointers

 !
 !*** activates degrees of freedom
 !
   e => heade                         !point to first element
   DO
     IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
     DO inode = 1,nnode               !for each node in the element
       IF( ifpre(1,e%lnods(inode)) == 1 ) &
           ifpre(1:2,e%lnods(inode)) = 0      !release DOFs
     END DO
     IF( .NOT.lumped )THEN
       gelem = gelem + 1                !update total number of elements
       gnods(1:4,gelem) = e%lnods(1:4)  !pass to global array
     END IF
     e => e%next                      !point to next element
   END DO

 RETURN
 END SUBROUTINE acvd17

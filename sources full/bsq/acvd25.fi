 SUBROUTINE acvd25(heade,lside,nelem,ifpre,nbs,bhead)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************

 IMPLICIT NONE

 ! dummy arguments

 LOGICAL, INTENT(IN OUT) :: lside !flag TRUE if topology already computed
 INTEGER (kind=4):: nelem,      & !(IN) number of elements
                    ifpre(:,:)    !(OUT) active DOFs
 INTEGER (kind=4):: nbs           !(OUT) number of branching sides
 TYPE (ele25), POINTER :: heade   !pointer to first element in the set
 TYPE (sideb), POINTER :: bhead   !pointer to first branching side

 ! local variables

 TYPE (ele25), POINTER  :: e                     !element pointers
 INTEGER (kind=4) inode

 !***  create topological array  lnods(5:12,nelem)

 IF(.NOT.lside) CALL nghb25(lside,nelem,nbs,heade,bhead)

 !*** activates degrees of freedom


 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
   DO inode = 1,4                   !for each node in the element
     IF( ifpre(1,e%lnods(inode)) == 1 ) & !to avoid to overwrite RIGID constraints
       ifpre(1:3,e%lnods(inode)) = 0      !release DOFs
   END DO
   !WRITE(55,"(12i5)")e%lnods
   e => e%next                      !point to next element
 END DO

 RETURN
 END SUBROUTINE acvd25

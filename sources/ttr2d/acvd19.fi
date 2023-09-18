 SUBROUTINE acvd19(ifpre,heade,lside,swapc,eside)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN OUT) :: lside !flag TRUE if topology already computed
 LOGICAL, INTENT(IN) :: swapc     !flag TRUE if elements with one neighbor will be swaped
 INTEGER (kind=4):: ifpre(:,:), & !(OUT) active DOFs
                    eside(:,:)    !(OUT) adjacent neighbor elements
 TYPE (ele19), POINTER :: heade   !pointer to first element in the set

 ! local variables
 INTEGER (kind=4)  inode
 TYPE (ele19), POINTER  :: e                     !element pointers

 !***  create topological array  lnods(4:6,nelem)

 IF(.NOT.lside) CALL toar19(heade,lside,swapc,eside)
 !
 !*** activates degrees of freedom
 !
 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
   DO inode = 1,3                   !for each node in the element
     IF( ifpre(1,e%lnods(inode)) == 1 ) &
       ifpre(1:2,e%lnods(inode)) = 0 !release DOFs
   END DO
   IF( therm )THEN
     DO inode = 1,3                      !for each node in the element
       iftmp(1:ndoft,e%lnods(inode)) = 0 !release thermal DOF
     END DO
   END IF
   e => e%next                      !point to next element
 END DO

 RETURN
 END SUBROUTINE acvd19

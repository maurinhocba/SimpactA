 SUBROUTINE updl18(heade,taile,oldlb,nelem)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: oldlb(:)      !(IN) previous strategy nodal labels
 INTEGER (kind=4):: nelem         !(OUT) number of elements
 TYPE (ele18), POINTER :: heade,taile   !pointer to first & last element in the set


 ! local variables
 INTEGER (kind=4) i,lab,chnode,ln(nnode),deleted
 TYPE (ele18), POINTER  :: e,ea                     !element pointers

 !     old local nodes are changed into new ones
 e => heade   !point to first element
 NULLIFY(ea)  !previous element is nothing
 nelem = 0    !initializes element count
 deleted = 0  !initializes number of deleted elements
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT    !if all element processed EXIT loop
   DO i=1,nnode                      !for each node
     lab  = oldlb(e%lnods(i))      !original label
     ln(i) = chnode(lab)           !new internal number
   END DO
   IF( ANY(ln == 0) )THEN            !if any node deleted
     deleted = deleted + 1           !increase number of deleted elements
     CALL del_ele18e (heade, taile, ea, e)  !delete element from the list
   ELSE
     nelem = nelem + 1               !increase number of remaining elements
     e%lnods = ln                    !modify inner connectivities
     ea => e                         !update previous element
     e => e%next                     !point to next element
   END IF
 END DO
 IF( deleted > 0 )  WRITE(lures,"(i10,' elements deleted', i10,' remaining')")deleted,nelem

 RETURN
 END SUBROUTINE updl18

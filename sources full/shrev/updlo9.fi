 SUBROUTINE updlo9(nnode,heade,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nnode,    & !number of nodes per element
                                 oldlb(:)    !previous strategy nodal labels
 TYPE (ele09), POINTER :: heade   !pointer to first element in the set

 INTEGER (kind=4) i,lab,chnode
 TYPE (ele09), POINTER  :: e                     !element pointers

!     old local nodes are changed into new ones

 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
   DO i=1,nnode
     lab  = oldlb(e%lnods(i))
     e%lnods(i) = chnode(lab)
   END DO
   e => e%next                      !point to next element
 END DO

 RETURN
 END SUBROUTINE updlo9

 SUBROUTINE updlo6(heade,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: oldlb(:)      !(IN) previous strategy nodal labels
 TYPE (ele06), POINTER :: heade   !pointer to first element in the set

 INTEGER (kind=4) i,lab,chnode
 TYPE (ele06), POINTER  :: e                     !element pointers

!     old local nodes are changed into new ones

 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
   DO i=1,nnode
     lab  = oldlb(e%lnods(i))
     e%lnods(i) = chnode(lab)
   END DO
 END DO

 RETURN
 END SUBROUTINE updlo6

 SUBROUTINE updl20(heade,oldlb,lside,origl)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN OUT) :: lside !flag TRUE if topology already computed
 LOGICAL, INTENT(IN OUT) :: origl !flag TRUE if original labels
 INTEGER (kind=4):: oldlb(:)      !(IN) previous strategy nodal labels
 TYPE (ele20), POINTER :: heade   !pointer to first element in the set

 ! local variables
 INTEGER (kind=4) i,lab,chnode
 TYPE (ele20), POINTER  :: e                     !element pointers

 !     old local nodes are changed into new ones
 e => heade   !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT    !if all element processed EXIT loop
   DO i=1,3                      !for each node
     IF( origl )THEN
       lab  = e%lnods(i)         !original label
     ELSE
       lab  = oldlb(e%lnods(i))  !original label
     END IF
     e%lnods(i) = chnode(lab)      !new internal number
   END DO
   e => e%next                       !point to next element
 END DO
 lside = .FALSE.

 RETURN
 END SUBROUTINE updl20

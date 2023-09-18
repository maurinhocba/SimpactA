 SUBROUTINE updl16(elset,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: oldlb(:)      !(IN) previous strategy nodal labels
 TYPE (ele16_set), POINTER :: elset   !pointer to first element in the set
 TYPE (ele16), POINTER :: heade   !pointer to first element in the set

 ! local variables
 INTEGER (kind=4) i,lab,chnode,deleted,nel
 TYPE (ele16), POINTER  :: e,anter               !element pointers

 !     old local nodes are changed into new ones
 e => elset%head   !point to first element
 NULLIFY( anter )
 nel = 0
 deleted = 0
 outer : DO
   IF (.NOT.ASSOCIATED (e) ) EXIT    !if all element processed EXIT loop
   DO i=1,nnode                      !for each node
     lab  = oldlb(e%lnods(i))      !original label
     IF( lab > 0 )THEN
       e%lnods(i) = chnode(lab)      !new internal number
     ELSE
       CALL del_ele16e (elset%head,elset%tail, anter, e )  !deletes element
       deleted = deleted + 1
       CYCLE outer
     END IF
   END DO
   nel = nel + 1
   anter => e
   e => e%next                       !point to next element

 END DO outer

 IF( deleted == 0  .AND. elset%quad )THEN

   e => elset%head   !point to first element
   DO nel=1,elset%nelem
     DO i=7,12                          !for each node
       IF(e%lnods(i)/=0) THEN          !if node exist
         lab  = oldlb(e%lnods(i))      !original label
         e%lnods(i) = chnode(lab)      !new internal number
       END IF
     END DO
     e => e%next                       !point to next element
   END DO

 ELSE

   elset%nelem = nel                   !new number of elements
   elset%lface = .FALSE.               !recompute connectivities

 END IF

 RETURN
 END SUBROUTINE updl16

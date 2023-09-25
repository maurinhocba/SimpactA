 SUBROUTINE updl24(elset,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: oldlb(:)      !(IN) previous strategy nodal labels
 TYPE (ele24_set), POINTER :: elset   !pointer to set

 ! local variables
 INTEGER (kind=4) i,j,k,l,lab,chnode,deleted,nel
 TYPE (ele24), POINTER  :: e,anter                  !element pointers
 TYPE (sideb), POINTER :: sb      !pointer to a branching side

 !     old local nodes are changed into new ones

 e => elset%head   !point to first element
 NULLIFY( anter )
 nel = 0
 deleted = 0
 OUTER : DO
   IF (.NOT.ASSOCIATED (e) ) EXIT    !if all element processed EXIT loop
   DO i=1,3                          !for each node
     IF(e%lnods(i)/=0) THEN          !if node exist
       lab  = oldlb(e%lnods(i))      !original label
       e%lnods(i) = chnode(lab)      !new internal number
     ELSE
       CALL del_ele24e (elset%head,elset%tail, anter, e )  !deletes element
       deleted = deleted + 1
       CYCLE OUTER
     END IF
   END DO
   nel = nel + 1
   anter => e
   e => e%next                       !point to next element
 END DO OUTER

 IF( deleted == 0 )THEN

   e => elset%head   !point to first element
   DO nel=1,elset%nelem
     DO i=4,6                           !for each node
       IF(e%lnods(i)/=0) THEN          !if node exist
         lab  = oldlb(e%lnods(i))      !original label
         e%lnods(i) = chnode(lab)      !new internal number
       END IF
     END DO
     e => e%next                       !point to next element
   END DO

   !     old local nodes are changed into new ones
   IF( elset%nbs > 0 )THEN
     sb => elset%bhead
     DO i=1,elset%nbs
       DO j=-1,sb%nn                     !for each node
         k = ABS(sb%lnods(j))            !node
         l = sb%lnods(j)/k               !sign
         lab  = oldlb(k)                 !original label
         sb%lnods(j) = l*chnode(lab)     !new internal number
       END DO
     END DO
   END IF

 ELSE

   elset%nelem = nel                   !new number of elements
   elset%lside = 0                     !recompute connectivities

 END IF

 RETURN
 END SUBROUTINE updl24

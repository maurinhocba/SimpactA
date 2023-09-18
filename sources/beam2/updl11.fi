 SUBROUTINE updl11(elset,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: oldlb(:)      !(IN) previous strategy nodal labels
 TYPE (ele11_set), POINTER :: elset

 ! local variables
 INTEGER (kind=4) i,j,lab,chnode,nel,deleted
 TYPE (ele11), POINTER  :: e,anter             !element pointers
 TYPE (nodei), POINTER :: nb,ni                !pointer to branching node

   !     old local nodes are changed into new ones
   deleted = 0       !initializes
   e => elset%head   !point to first element
   NULLIFY (anter)
   nel = 0
   outer : DO
     IF (.NOT.ASSOCIATED (e) ) EXIT    !if all element processed EXIT loop
     DO i=2,3                          !for each node
       IF(e%lnods(i) > 0) THEN           !if node exist
         lab = oldlb(e%lnods(i))      !original label
         lab = chnode(lab)      !new internal number
         IF( lab > 0 )THEN
           e%lnods(i) = chnode(lab)      !new internal number
         ELSE
           CALL del_ele11e (elset%head, elset%tail, anter, e)
           deleted = deleted + 1
           CYCLE outer
         END IF
       END IF
     END DO
     anter => e
     e => e%next                       !point to next element
     nel = nel + 1
   END DO outer

   IF( deleted > 0 )THEN
     elset%nelem = nel
     ! deallocate branching nodes database
     IF( elset%nbn > 0 )THEN
       ni => elset%nhead
       DO i=1,elset%nbn
         nb => ni
         DEALLOCATE(nb%lnods,nb%alph0,nb%fc,nb%bb,nb)
         ni => nb%next
         DEALLOCATE(nb)
       END DO
     END IF
     elset%lside = .FALSE.
   ELSE
     !     old local nodes are changed into new ones
     e => elset%head   !point to first element
     DO i=1,nel
       DO j=1,4,3                        !for each node
         IF(e%lnods(j) > 0) THEN           !if node exist
           lab = oldlb(e%lnods(j))      !original label
           e%lnods(j) = chnode(lab)            !new internal number
         END IF
       END DO
       e => e%next                       !point to next element
     END DO
     ni => elset%nhead
     DO i=1,elset%nbn
       DO j=0,ni%nn
         lab  = oldlb(ni%lnods(i))      !original label
         ni%lnods(i) = chnode(lab)      !new internal number
       END DO
       ni => ni%next
     END DO
   END IF

 RETURN
 END SUBROUTINE updl11

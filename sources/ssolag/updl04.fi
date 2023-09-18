 SUBROUTINE updl04(elset,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: oldlb(:)      !(IN) previous strategy nodal labels
 TYPE (ele04_set), POINTER :: elset   !pointer to set


 ! local variables
 INTEGER (kind=4) i,lab,chnode,deleted,nel
 TYPE (ele04), POINTER  :: e,anter                  !element pointers

 !     old local nodes are changed into new ones

 e => elset%head   !point to first element
 NULLIFY( anter )
 nel = 0
 deleted = 0
 OUTER : DO
   IF (.NOT.ASSOCIATED (e) ) EXIT    !if all element processed EXIT loop
   DO i=1,nnode                      !for each node
     IF(e%lnods(i)/=0) THEN          !if node exist
       lab  = oldlb(e%lnods(i))      !original label
       e%lnods(i) = chnode(lab)      !new internal number
     ELSE
       CALL del_ele04e (elset%head,elset%tail, anter, e )  !deletes element
       deleted = deleted + 1
       CYCLE OUTER
     END IF
   END DO
   nel = nel + 1
   anter => e
   e => e%next                       !point to next element
 END DO OUTER

 IF( deleted /= 0  ) elset%nelem = nel    !new number of elements

 RETURN
 END SUBROUTINE updl04

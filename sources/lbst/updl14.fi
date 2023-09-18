 SUBROUTINE updl14(elset,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: oldlb(:)      !(IN) previous strategy nodal labels
 TYPE (ele14_set), POINTER :: elset   !pointer to set

 ! local variables
 INTEGER (kind=4) i,lab,chnode,deleted,nel
 TYPE (ele14), POINTER  :: e,anter                  !element pointers


 !     old local nodes are changed into new ones

 e => elset%head   !point to first element
 NULLIFY( anter )
 nel = 0
 deleted = 0
 OUTER : DO
   IF (.NOT.ASSOCIATED (e) ) EXIT    !if all element processed EXIT loop
   DO i=1,3                          !for each node
      IF( elset%origl ) THEN  ! Original label
        lab = e%lnods(i)
      ELSE                    ! previously changed
        lab = oldlb(e%lnods(i))  ! old label recovery
      END IF
      lab = chnode(lab) ! update internal label
     IF( lab /= 0) THEN             !if node exist
       e%lnods(i) = lab              !new internal number
     ELSE
       CALL del_ele14e (elset%head,elset%tail, anter, e )  !deletes element
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

 ELSE

   elset%nelem = nel                   !new number of elements
   elset%lside = 0                     !recompute connectivities

 END IF

 RETURN
 END SUBROUTINE updl14

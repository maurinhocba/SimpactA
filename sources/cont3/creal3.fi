 SUBROUTINE creal3(sname)
 !
 !     RENUMBER nodes in surface contact Data_Base
 !
 IMPLICIT NONE
 ! routine arguments
 CHARACTER(len=*),INTENT(IN):: sname            !surface name

 ! local variables
 INTEGER (kind=4) :: i,n,isurf,lab,chnode
 INTEGER (kind=4), PARAMETER :: nnseg = 3
 TYPE (surf3_db), POINTER :: surf


 ! Updates surface connectivities

 surf => shead                               ! initializes
 DO isurf=1,nsurf                            ! for each pair
   IF( TRIM(surf%sname) /= TRIM(sname) )THEN         ! if not the modified surface
     IF( ASSOCIATED(surf%lcnod) )THEN
       DO i=1,surf%ncnod
         lab  = oldlb(surf%lcnod(i))
         surf%lcnod(i) = chnode(lab)
       END DO
     END IF
     IF( ASSOCIATED(surf%lcseg) )THEN
       DO n=1,surf%nsegm
         DO i=1,nnseg
           lab  = oldlb(surf%lcseg(i,n))
           surf%lcseg(i,n) = chnode(lab)
         END DO
       END DO
     END IF
     IF( ASSOCIATED(surf%lcseb) )THEN
       DO n=1,surf%nsegm
         DO i=1,nnseg
           lab  = oldlb(surf%lcseb(i,n))
           surf%lcseb(i,n) = chnode(lab)
         END DO
       END DO
     END IF

   END IF
   surf => surf%next                         ! next pair
 END DO

 RETURN
 END SUBROUTINE creal3

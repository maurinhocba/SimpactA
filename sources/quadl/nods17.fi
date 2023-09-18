 SUBROUTINE nods17(npoin, numpo, heade, nodset, label)
 !***********************************************************************
 ! Extracts node numbers (labels) into the nodal set
 !***********************************************************************
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4):: npoin, numpo, label(:)
   INTEGER(kind=4),POINTER:: nodset(:)
   TYPE(ele17),POINTER:: heade
   !--- Local variables
   INTEGER(kind=4):: i, k, n
   INTEGER(kind=4),ALLOCATABLE:: iwork(:)
   TYPE(ele17),POINTER:: elm

   !Marks used nodes
   ALLOCATE (iwork(npoin))
   iwork(1:npoin) = 0


    elm => heade
    DO
      IF (.NOT.ASSOCIATED(elm)) EXIT
      iwork(elm%lnods(1:4)) = 1
      elm => elm%next
    END DO

   !Create a nodal set
   numpo = SUM(iwork)
   ALLOCATE (nodset(numpo))
   k = 0
   DO i=1,npoin
     IF (iwork(i) == 0) CYCLE
     k = k + 1
     n = label(i)
     nodset(k) = n
   END DO

   DEALLOCATE (iwork)

 RETURN
 END SUBROUTINE nods17

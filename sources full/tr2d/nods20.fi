 SUBROUTINE nods20(npoin, numpo, heade, nodset, label)
 !***********************************************************************
 ! Extracts node numbers (labels) into the nodal set
 !***********************************************************************
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4):: npoin, numpo, label(:)
   INTEGER(kind=4),POINTER:: nodset(:)
   TYPE(ele20),POINTER:: heade
   !--- Local variables
   INTEGER(kind=4):: i, k, n
   INTEGER(kind=4),ALLOCATABLE:: iwork(:)
   TYPE(ele20),POINTER:: elm

   !Marks used nodes
   ALLOCATE (iwork(npoin))
   iwork(1:npoin) = 0

   IF (.NOT.r_elm_zone) THEN
     elm => heade
     DO
       IF (.NOT.ASSOCIATED(elm)) EXIT
       iwork(elm%lnods(1:3)) = 1
       elm => elm%next
     END DO
   ELSE ! if zone remeshing then
     elm => heade
     DO
       IF (.NOT.ASSOCIATED(elm)) EXIT
       IF(elm%repla) iwork(elm%lnods(1:3)) = 1
       elm => elm%next
     END DO
   END IF

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
 END SUBROUTINE nods20

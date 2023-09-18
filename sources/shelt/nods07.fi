SUBROUTINE nods07(nelem, heade, label)
!***********************************************************************
!
!*** extracts node numbers (labels) into the nodal set
!
!***********************************************************************
USE ctrl_db, ONLY :npoin
USE meshmo_db,ONLY:  nodset,numpo
IMPLICIT NONE
  !--- Dummy variables
  INTEGER(kind=4):: nelem, label(:)
  TYPE(ele07),POINTER:: heade      !pointer to first element
  !--- Local Variables
  INTEGER(kind=4):: i, k, n, ielem
  INTEGER(kind=4),ALLOCATABLE:: iwork(:)
  TYPE(ele07),POINTER:: e

  !Marks used nodes
  ALLOCATE( iwork(npoin) )     !auxiliar array

  iwork = 0                      !initializes
  e => heade                     !point to first element
  DO ielem=1,nelem               !for each element in the set
    iwork(e%lnods(1:6)) = 1      !set element nodes as existent
    e => e%next                  !point to next element
  END DO

  !Create a nodal set
  numpo = SUM(iwork)             !number of points in the set
  ALLOCATE (nodset(numpo))       !reserve memory for the labels
  k = 0                          !initializes counter to NODSET
  DO i=1,npoin                   !for each node in the mesh
    IF (iwork(i) == 1) THEN      !if node is in the set
      k = k + 1                  !increase counter in NODSET
      n = label(i)               !nodal label
      nodset(k) = n              !store nodal label
    END IF
  END DO

  DEALLOCATE (iwork)             !release auxiliar array

RETURN
END SUBROUTINE nods07

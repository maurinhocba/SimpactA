      SUBROUTINE nods24 (npoin, nelem, numpo, nodset, heade, label)
     &
!***********************************************************************
!
!*** extracts node numbers (labels) into the nodal set
!
!***********************************************************************
      IMPLICIT NONE
!     Arguments
      INTEGER (kind=4) :: npoin, nelem, numpo, label(:)
      INTEGER (kind=4), POINTER :: nodset(:)
      TYPE (ele24), POINTER  :: heade      !pointer to first element

!     Local Variables
      INTEGER (kind=4) i,k,n,ielem
      INTEGER (kind=4), ALLOCATABLE :: iwork(:)
      TYPE (ele24), POINTER :: e

!     marks used nodes

      ALLOCATE ( iwork(npoin) )     !auxiliar array

      iwork = 0                      !initializes
      e => heade                     !point to first element
      DO ielem=1,nelem               !for each element in the set
        iwork(e%lnods(1:3)) = 1      !set element nodes as existent
        e => e%next                  !point to next element
      END DO

!     create a nodal set

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
      END SUBROUTINE nods24

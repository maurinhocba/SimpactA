      SUBROUTINE boun13 ( heade, nelem)
!******************************************************************
!
!*** get the boundary definition from the element set
!
!******************************************************************
      IMPLICIT NONE

      INTEGER (kind=4), INTENT(IN) :: nelem !number of elements in the set
      TYPE (ele13), POINTER :: heade        !pointer to first element

      INTEGER (kind=4) :: ie
      TYPE (ele13), POINTER :: elem
      INTEGER(kind=4), ALLOCATABLE :: lnods(:,:)

      INTERFACE
        INCLUDE 'surf03.h'
      END INTERFACE

      ALLOCATE( lnods(3,nelem) ) ! auxiliar space to store connectivities

      elem => heade              ! point to first element
      DO ie=1,nelem              ! for every element in the list
        lnods(1:3,ie) = elem%lnods(1:3)  ! assign connectivities
        elem => elem%next                ! point to next element in the list
      END DO

      CALL surf03 ( lnods, 3, nelem)     ! compute boundary

      DEALLOCATE ( lnods )               ! release auxiliar space

      RETURN
      END SUBROUTINE boun13

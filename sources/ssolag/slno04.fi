      SUBROUTINE slno04 (nelem, heade, lnods )
!***********************************************************************
!
!*** extracts connectivities into the nodal set
!
!***********************************************************************
      IMPLICIT NONE
!     Arguments
      INTEGER (kind=4) :: nelem
      INTEGER (kind=4), POINTER :: lnods(:,:)
      TYPE (ele04), POINTER :: heade

!     Local Variables
      INTEGER (kind=4) ielem,nn
      TYPE (ele04), POINTER :: e

      IF( ASSOCIATED(lnods) )THEN
        nn    = SIZE(lnods,1)
        ielem = SIZE(lnods,2)
        IF( ielem /= nelem .OR. nn /= nnode)THEN
          DEALLOCATE(lnods)
          ALLOCATE(lnods(nnode,nelem))
        END IF
      ELSE
        ALLOCATE(lnods(nnode,nelem))
      END IF
      e => heade
      DO ielem=1,nelem
        lnods(:,ielem) = e%lnods(1:nnode)
        e => e%next
      END DO
      RETURN
      END SUBROUTINE slno04

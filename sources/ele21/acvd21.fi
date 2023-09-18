      SUBROUTINE acvd21 (nelem, mxnods, lm)
!***********************************************************************
!
!*** updates internal node numbering and activates degrees of freedom
!
!***********************************************************************
      USE npo_db, ONLY : iftmp
      USE ctrl_db, ONLY : therm
      IMPLICIT NONE
!     Arguments
      INTEGER (kind=4) :: nelem,mxnods,lm(:,:)

!     Local Variables
      INTEGER (kind=4) ielem

!*** activates degrees of freedom

      IF( therm )THEN
        DO ielem=1,nelem
          iftmp(1,lm(1:mxnods,ielem)) = 0
        END DO
      END IF

      RETURN
      END SUBROUTINE acvd21

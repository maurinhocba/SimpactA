SUBROUTINE acvd22 ( nelem, nnode, lnods)

!*** updates internal node numbering and activates degrees of freedom

  USE npo_db, ONLY : iftmp
  USE ctrl_db, ONLY : therm
  IMPLICIT NONE
! Arguments
  INTEGER (kind=4) :: nelem,nnode,lnods(:,:)

! Local Variables
  INTEGER (kind=4) ielem

!activates degrees of freedom

  IF( therm )THEN
    DO ielem=1,nelem
      iftmp(1,lnods(1:nnode,ielem)) = 0
    END DO
  END IF

  RETURN
END SUBROUTINE acvd22

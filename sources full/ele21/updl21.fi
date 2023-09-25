      SUBROUTINE updl21 (nelem, mxnods, lm, oldlb)
!***********************************************************************
!
!*** updates internal node numbering
!
!***********************************************************************
      IMPLICIT NONE
!     Arguments
      INTEGER (kind=4) :: nelem,mxnods,oldlb(:),lm(:,:)

!     Local Variables
      INTEGER (kind=4) ielem,i,lab,chnode

!     old local nodes are changed into new ones

      DO ielem=1,nelem
        DO i=1,mxnods
          lab  = oldlb(lm(i,ielem))
          lm(i,ielem) = chnode(lab)
        END DO
      END DO

      RETURN
      END SUBROUTINE updl21

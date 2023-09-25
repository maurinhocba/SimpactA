SUBROUTINE updl22 ( nelem, nnode, lnods, oldlb)

!*** updates internal node numbering and activates degrees of freedom

  IMPLICIT NONE
! Arguments
  INTEGER (kind=4) :: nelem,nnode,lnods(:,:),oldlb(:)

! Local Variables
  INTEGER (kind=4) ielem,i,lab,chnode

! old local nodes are changed into new ones

  DO ielem=1,nelem
    DO i=1,nnode
      lab  = oldlb(lnods(i,ielem))
      lnods(i,ielem) = chnode(lab)
    END DO
  END DO

  RETURN
END SUBROUTINE updl22

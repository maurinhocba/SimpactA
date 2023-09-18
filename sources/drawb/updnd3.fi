SUBROUTINE updnd3(oldlb,nmn,nrtm,msr,irect,keep)
!***********************************************************************
!
!*** updates internal node numbering for 3D drawbeads
!
!***********************************************************************
  IMPLICIT NONE
! Arguments
  INTEGER (kind=4), INTENT(IN) :: oldlb(:), & ! original labels
                                  nmn         ! number of nodes in drawbead line
  INTEGER (kind=4), INTENT(IN OUT) ::  msr(:),   & ! drawbead connectivities
                                       nrtm,     & ! number of segments in surface
                                       irect(:,:)  ! surface connectivities
  LOGICAL, INTENT(OUT) :: keep

  ! Local Variables
  INTEGER (kind=4) i,j,k,lab,chnode,ln(4)

  ! old local nodes are changed into new ones
  ! drawbead connectivity
  DO i=1,nmn
    lab  = oldlb(msr(i))
    msr(i) = chnode(lab)
  END DO
  keep = ALL( msr > 0 )
  ! surface connectivity
  k = 0
  DO j=1,nrtm
    DO i=1,4
      lab  = oldlb(irect(i,j))
      ln(i) = chnode(lab)
      irect(i,j) = chnode(lab)
    END DO
    IF( ALL(ln(1:3) > 0) )THEN
      k = k+1
      irect(:,k) = ln
    END IF
  END DO
  nrtm = k
  keep = keep .AND. nrtm > 0

  RETURN
END SUBROUTINE updnd3

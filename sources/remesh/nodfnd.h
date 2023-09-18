   SUBROUTINE nodfnd(lnbnd,nlbnd,cdbnd,oldlb,ndnum,coorn,found)
  !**************************************************************
  !  This subroutine check the boundary nodes in zone remeshing
  !  and save the old & new nodes label for paste zone elem-patch
  !**************************************************************
  USE ctrl_db,ONLY: ndime
  IMPLICIT NONE
  ! dummy variables
  INTEGER(kind=4),POINTER       :: lnbnd(:), &  ! old boundary nodes
                                   nlbnd(:), &  ! segments in boundary lines
                                   oldlb(:)     ! old labels 
  INTEGER(kind=4),INTENT(INOUT) :: ndnum        ! new label of node
  REAL(kind=8),POINTER          :: cdbnd(:,:)   ! old boundary coordinates
  REAL(kind=8),INTENT(IN)       :: coorn(ndime) ! new node coordinates
  LOGICAL, INTENT(INOUT)        :: found      

  END SUBROUTINE nodfnd

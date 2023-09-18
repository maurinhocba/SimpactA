SUBROUTINE verifd (nreqd,nprqd,numpo,ndime,numpn,coord,coorn, &
                   nodset,label)

! This routine modifies points for requested output
! if a point belongs to a remeshed mesh,
! it is replaced by the nearest point in the new mesh

IMPLICIT NONE
! dummy arguments
INTEGER (kind=4), INTENT(IN) :: numpo,ndime,numpn, &
                                nodset(:),label(:)
REAL (kind=8), INTENT(IN) :: coord(:,:),coorn(:,:)
INTEGER (kind=4), INTENT(IN OUT) :: nreqd,nprqd(:)

END SUBROUTINE verifd


SUBROUTINE  verife(nreqd,nprqd,maxnn,rzone,nodlb,oldlb)

! This routine modifies points for requested output

IMPLICIT NONE
! dummy arguments
INTEGER (kind=4), INTENT(IN) :: maxnn
INTEGER (kind=4), INTENT(IN OUT) :: nreqd,nprqd(:)
INTEGER (kind=4), INTENT(IN) :: nodlb(:,:),oldlb(:)
LOGICAL         , INTENT(IN) :: rzone

END SUBROUTINE verife  

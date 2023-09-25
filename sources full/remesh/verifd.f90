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

! local variables
INTEGER (kind=4) :: i,j,k,n,m,nl
REAL (kind=8) :: xn(ndime),d(ndime),dmin,dd
LOGICAL :: found,modif

  modif = .FALSE.           !initializes
  DO i=1,nreqd              !for each requested point for output
    found = .FALSE.         !initializes flag
    n = nprqd(i)            !old internal number
    nl = label(n)           !label in the original mesh
    DO j=1,numpo            !for each node in the deleted SET
      found = ( nl == nodset(j) )    !check labels
      IF( found )THEN                !if node will be deleted
        modif = .TRUE.
        xn = coord(:,n)              !original coordinates (change j by n) wbc
        k  = 0                       !initializes nearest node
        dmin = HUGE(1d0)             !initializes minimum distance
        DO m=1,numpn                 !for each node in the new mesh
          d = coorn(:,m) - xn        !distance vector
          dd = DOT_PRODUCT(d,d)      !distance (squared)
          IF( dd < dmin )THEN        !compare with previous value
            k = m                    !new nearest node
            dmin = dd                !new distance
            IF( dmin == 0d0 )EXIT    !if node coincide, exit search
          END IF
        END DO
        nprqd(i) = - k               !store with minus sign to be modified
        EXIT
      END IF
    END DO
  END DO
  IF( modif ) nreqd = -nreqd

  RETURN

END SUBROUTINE verifd


SUBROUTINE  verife(nreqd,nprqd,maxnn,rzone,nodlb,oldlb)

! This routine modifies points for requested output

IMPLICIT NONE
! dummy arguments
INTEGER (kind=4), INTENT(IN) :: maxnn
INTEGER (kind=4), INTENT(IN OUT) :: nreqd,nprqd(:)
INTEGER (kind=4), INTENT(IN) :: nodlb(:,:),oldlb(:)
LOGICAL         , INTENT(IN) :: rzone

INTEGER (kind=4) :: i,chnode,ndnum(1),nreqo

IF(nreqd < 0) nreqd = -nreqd
DO i=1,nreqd
  IF(nprqd(i) < 0)THEN
    nreqo = -nprqd(i)+maxnn !new node label
    IF(rzone)THEN ! if zone remesh: check if it's a border node
      ndnum = MAXVAL(nodlb(2,:),MASK= nodlb(1,:)==nreqo)
      IF(ndnum(1) > 0) nreqo = ndnum(1) !border node label
    END IF
  ELSE
    nreqo =  oldlb(nprqd(i)) !old node label
  END IF
  !IF(nprqd(i) < 0) nprqd(i) = chnode(nreqo)
  nprqd(i) = chnode(nreqo) !new position in node array
END DO
RETURN

END SUBROUTINE verife

  SUBROUTINE nodfnd(lnbnd,nlbnd,cdbnd,oldlb,ndnum,coorn,found)
  !**************************************************************
  !
  !  This subroutine check the boundary nodes in zone remeshing
  !  and save the old & new nodes label for paste zone elem-patch
  !
  !**************************************************************
  USE ctrl_db,ONLY: ndime
  IMPLICIT NONE
  ! dummy variables
  INTEGER(kind=4),POINTER       :: lnbnd(:), &  ! index boundary nodes
                                   nlbnd(:), &  ! closed lines in boundary
                                   oldlb(:)     ! old node labels
  INTEGER(kind=4),INTENT(INOUT) :: ndnum        ! new label of node
  REAL(kind=8),POINTER          :: cdbnd(:,:)   ! old boundary coordinates
  REAL(kind=8),INTENT(IN)       :: coorn(ndime) ! new node coordinates
  LOGICAL, INTENT(INOUT)        :: found
  ! local variables
  INTEGER(kind=4)  :: i,j,i1,i2
  REAL(kind=8)     :: xmin,xmax,tol,xdis(ndime),rdis

  REAL(kind=8), PARAMETER ::  eps =1.0d-6  !

  i1    = 0      ! initializes
  i2    = 0      !
  found = .FALSE.
  L10: DO j = 1,SIZE(nlbnd) ! for each distorted zone
    i1 = i2 + 1        ! lower segment in boundary line
    i2 = i2 + nlbnd(j) ! upper           "
    ! check if node is enclosed by boundary bounding box
    DO i = 1,ndime
      xmin = MINVAL(cdbnd(i,i1:i2)) !minimum value for this coordinate
      xmax = MAXVAL(cdbnd(i,i1:i2)) !maximum           "
      tol  = 0.1d0*ABS(xmax-xmin) !tolerance for bounding box
      IF( coorn(i) < (xmin - tol)  .OR. &
          coorn(i) > (xmax + tol) ) CYCLE L10 ! node is a outlier
    END DO

    ! check if node coincide with the boundary nodes
    DO i = i1,i2  ! for each node in boundary
      xdis(:) = coorn(:) - cdbnd(:,i)
      rdis    = SQRT(xdis(1)*xdis(1) + xdis(2)*xdis(2))
      IF( ABS(rdis) < eps )THEN
        ndnum = oldlb(lnbnd(i))  !restore the old label
        found = .TRUE.   !mark as found
        EXIT L10       !and exit
      END IF
    END DO
  END DO L10

  RETURN
  END SUBROUTINE nodfnd

      SUBROUTINE slno18 (nelem, nnode, heade, lnods, reorder )
!***********************************************************************
!
!*** extracts connectivities into the nodal set
!
!***********************************************************************
      USE npo_db, ONLY : coord
      IMPLICIT NONE
!     Arguments
      INTEGER (kind=4) :: nelem, nnode
      INTEGER (kind=4), POINTER :: lnods(:,:)
      TYPE (ele18), POINTER :: heade
      LOGICAL, INTENT(IN) :: reorder

!     Local Variables
      INTEGER (kind=4) ielem,nn,i, np(1)
      TYPE (ele18), POINTER :: e
      REAL (kind=8) :: l(3),x(3,4),ls(3)
      INTEGER (kind=4), PARAMETER :: la(2,3) = (/ 1,2, 2,3, 2,4 /)

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
        IF( reorder )THEN
          ! compute shortest direction
          x = coord(:,e%lnods((/ 1,2,3,6 /)))
          DO i=1,3
            ls = x(:,la(2,i)) - x(:,la(1,i))
            l(i) = DOT_PRODUCT(ls,ls)
          END DO
          np = MINLOC(l)              !shortest direction
          SELECT CASE(np(1))          !set first and last faces normal to shortest direction
          CASE(1)
            lnods(:,ielem) = e%lnods((/1,4,8,5,2,3,7,6 /))
          CASE(2)
            lnods(:,ielem) = e%lnods((/4,3,7,8,1,2,6,5 /))
          CASE(3)
            lnods(:,ielem) = e%lnods
          END SELECT
          e%lnods(:) = lnods(:,ielem)
        ELSE
          lnods(:,ielem) = e%lnods
        END IF
        e => e%next
      END DO
      RETURN
      END SUBROUTINE slno18

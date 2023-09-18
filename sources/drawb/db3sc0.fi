SUBROUTINE db3sc0 (x, xl, nrtm, ldbs, irect, &
 &                   nstatold, ndbpoin, fdat, nnode, ncrosi)
  !-------------------------------------------------------------
  !  search for the intersection of a drawbead line with a sheet
  !-------------------------------------------------------------

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT (IN) :: nrtm
  INTEGER :: ndbpoin
  INTEGER :: nnode
  INTEGER :: ncrosi
  INTEGER :: ldbs (:,:)
  INTEGER :: irect (:,:)
  INTEGER :: nstatold (:)
  REAL (kind=8) :: x (:,:)
  REAL (kind=8) :: xl (:, :)
  REAL (kind=8) :: fdat (:, :)

  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s

  INTEGER :: i
  !-----------------------------------------------

  ncrosi = 0
  DO i = 1, nrtm
!#if UNIX
!        CALL chqua0 (x, xl, nstatold(i), ldbs(:, i), &
!     &    irect(:, i), ndbpoin, fdat(:, i), nnode)
!#else
    CALL chqua0 (x, xl, nstatold(i), ldbs(1:2, i), &
      irect(1:4,i), ndbpoin, fdat(1:3, i), nnode)
!#endif
    ! count intersected elements
    IF (nstatold(i) == 1) ncrosi = ncrosi + 1
  END DO
  RETURN
END SUBROUTINE db3sc0

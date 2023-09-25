SUBROUTINE inidb3 ( irect, msr, iwork, nseg, lmsr,         &
 &    coora,xl, ldbs, nstatold, fdat, ndbpoin, nrt, nmn, nnode,&
 &       ncrosi)
  !-----------------------------------------------------------------
  ! initialization of the drawbead-line-and-sheet database
  !-----------------------------------------------------------------

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT (IN) :: ndbpoin, & !number of nodes along DB line
                          nrt,     & !number of elements to check
                          nmn,     & !number of nodes in the mesh
                          nnode      !number of nodes in connect. 
  INTEGER :: ncrosi
  INTEGER, INTENT (IN) :: irect (:,:)  !element connectivities
  INTEGER :: msr (:)
  INTEGER :: iwork (:)
  INTEGER :: nseg (:)
  INTEGER, INTENT(OUT) :: lmsr (:)
  INTEGER :: ldbs (:, :)
  INTEGER :: nstatold (:)
  REAL (kind=8), INTENT (IN) :: coora (:,:)  !present coordinates
  REAL (kind=8), INTENT (IN) :: xl (:, :)    !DB line coordinates
  REAL (kind=8) :: fdat (:,:)

  !-----------------------------------------------------------------

  !  create the database of sheet segments surrounding every node

  lmsr = 0
  CALL linksd (irect, lmsr, msr, nseg, nrt, nmn, iwork)

  !  search for the intersection of a drawbead line with a sheet

  CALL db3sc0 (coora, xl, nrt, ldbs, irect, nstatold, ndbpoin, fdat,  &
 &             nnode, ncrosi)

  RETURN
END SUBROUTINE inidb3

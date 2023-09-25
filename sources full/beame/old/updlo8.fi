SUBROUTINE updlo8(nelem,nnode,heade,oldlb)

  !*** updates internal node numbering

  IMPLICIT NONE
  INTEGER (kind=4) :: nelem,nnode,oldlb(:)
  TYPE (ele08), POINTER :: heade

  INTEGER (kind=4) i,ielem,lab,chnode
  TYPE (ele08), POINTER :: e

  ! old local nodes are changed into new ones

  e => heade
  DO ielem=1,nelem
    DO i=1,nnode
      lab        = oldlb(e%lnods(i))
      e%lnods(i) = chnode(lab)
    END DO
    e => e%next
  END DO

  RETURN

END SUBROUTINE updlo8

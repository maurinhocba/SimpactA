SUBROUTINE updlo8(nelem,nnode,heade,taile,oldlb)

  !*** updates internal node numbering

  IMPLICIT NONE
  INTEGER (kind=4) :: nelem,nnode,oldlb(:)
  TYPE (ele08), POINTER :: heade,taile

  INTEGER (kind=4) i,ielem,lab,chnode,deleted,nel
  TYPE (ele08), POINTER :: e,anter

  ! old local nodes are changed into new ones

  e => heade
  deleted = 0       !initializes
  NULLIFY (anter)
  nel = 0
  outer : DO ielem=1,nelem
    DO i=1,nnode
      lab        = oldlb(e%lnods(i))
      IF( lab > 0 )THEN
        e%lnods(i) = chnode(lab)
      ELSE
        CALL del_ele08e (heade, taile, anter, e)
        deleted = deleted + 1
        CYCLE outer
      END IF
    END DO
    e => e%next
  END DO outer
  IF( deleted > 0 )nelem=nel

  RETURN

END SUBROUTINE updlo8

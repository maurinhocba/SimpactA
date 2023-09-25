MODULE vreal_db
IMPLICIT NONE

  TYPE vrea_db   !Definition: List of nodes (edges, level structure, etc)
    REAL(kind=8):: val
    TYPE (vrea_db),POINTER:: next
  END TYPE vrea_db

CONTAINS

!-----------------------------------------------------------------------
  SUBROUTINE new_rea(new)
!-----------------------------------------------------------------------
  !initialize a list of edges
  TYPE(vrea_db),POINTER:: new

  ALLOCATE(new)
  new%val = 0d0
  NULLIFY(new%next)

  RETURN
  END SUBROUTINE new_rea

!-----------------------------------------------------------------------
  SUBROUTINE add_rea(new,head,last)
!-----------------------------------------------------------------------
  !This subroutine adds data to the end of the list
  TYPE(vrea_db),POINTER:: new, head, last

  IF (.NOT.ASSOCIATED(head)) THEN
    !list is empty, start it
    head => new
    last => new
  ELSE
    !Check if a list is empty
    last%next => new
    last => new
  END IF

  RETURN
  END SUBROUTINE add_rea

!-----------------------------------------------------------------------
  SUBROUTINE dalloc_rea(head,last)
!-----------------------------------------------------------------------
  !Delete a list
  TYPE(vrea_db),POINTER:: head,last
  !Local variables
  TYPE(vrea_db),POINTER:: aux

  IF (.NOT.ASSOCIATED(head)) RETURN

  !Search the starting point
  NULLIFY(last)
  DO
    IF (.NOT.ASSOCIATED(head)) EXIT
    aux => head%next
    DEALLOCATE(head)
    head => aux
  END DO

  RETURN
  END SUBROUTINE dalloc_rea

END MODULE vreal_db

MODULE vstru_db
IMPLICIT NONE

  TYPE vstr_db   !Definition: List of nodes (edges, level structure, etc)
    INTEGER(kind=4):: dim
    REAL(kind=8),POINTER:: vec(:)
    TYPE (vstr_db),POINTER:: next
  END TYPE vstr_db

CONTAINS

!-----------------------------------------------------------------------
  SUBROUTINE new_vec(dim,new)
!-----------------------------------------------------------------------
  !initialize a list of edges
  INTEGER(kind=4),INTENT(IN):: dim
  TYPE(vstr_db),POINTER:: new

  ALLOCATE(new)
  new%dim = dim
  ALLOCATE(new%vec(dim))
  new%vec(1:dim) = 0d0
  NULLIFY(new%next)

  RETURN
  END SUBROUTINE new_vec

!-----------------------------------------------------------------------
  SUBROUTINE add_vec(new,head,last)
!-----------------------------------------------------------------------
  !This subroutine adds data to the end of the list
  TYPE(vstr_db),POINTER:: new, head, last

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
  END SUBROUTINE add_vec

!-----------------------------------------------------------------------
  SUBROUTINE del_vec(elm)
!-----------------------------------------------------------------------
  !Delete a list
  TYPE(vstr_db),POINTER:: elm

  IF (.NOT.ASSOCIATED(elm)) RETURN

  DEALLOCATE(elm%vec)
  DEALLOCATE(elm)

  RETURN
  END SUBROUTINE del_vec

!-----------------------------------------------------------------------
  SUBROUTINE dalloc_vec(head,last)
!-----------------------------------------------------------------------
  !Delete a list
  TYPE(vstr_db),POINTER:: head,last
  !Local variables
  TYPE(vstr_db),POINTER:: aux

  IF (.NOT.ASSOCIATED(head)) RETURN

  !Search the starting point
  NULLIFY(last)
  DO
    IF (.NOT.ASSOCIATED(head)) EXIT
    aux => head%next
    CALL del_vec(head)
    head => aux
  END DO

  RETURN
  END SUBROUTINE dalloc_vec

END MODULE vstru_db

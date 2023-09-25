MODULE vstrg_db
IMPLICIT NONE

  TYPE strg_db   !Definition: List of nodes (edges, level structure, etc)
    INTEGER(kind=4):: len
    CHARACTER(len=1),POINTER:: strng(:)
    TYPE (strg_db),POINTER:: next
  END TYPE strg_db

CONTAINS

!-----------------------------------------------------------------------
  SUBROUTINE new_strg(len,new)
!-----------------------------------------------------------------------
  !initialize a list of edges
  INTEGER(kind=4),INTENT(IN):: len
  TYPE(strg_db),POINTER:: new

  ALLOCATE(new)
  new%len = len
  ALLOCATE(new%strng(len))
  new%strng(1:len) = ' '
  NULLIFY(new%next)

  RETURN
  END SUBROUTINE new_strg

!-----------------------------------------------------------------------
  SUBROUTINE add_strg(new,head,last)
!-----------------------------------------------------------------------
  !This subroutine adds data to the end of the list
  TYPE(strg_db),POINTER:: new, head, last

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
  END SUBROUTINE add_strg

!-----------------------------------------------------------------------
  SUBROUTINE dalloc_strg(head,last)
!-----------------------------------------------------------------------
  !Delete a list
  TYPE(strg_db),POINTER:: head,last
  !Local variables
  TYPE(strg_db),POINTER:: aux

  IF (.NOT.ASSOCIATED(head)) RETURN

  !Search the starting point
  NULLIFY(last)
  DO
    IF (.NOT.ASSOCIATED(head)) EXIT
    aux => head%next
    DEALLOCATE(head%strng)
    DEALLOCATE(head)
    head => aux
  END DO

  RETURN
  END SUBROUTINE dalloc_strg

!-----------------------------------------------------------------------
  SUBROUTINE wr_strg(elm,string)
!-----------------------------------------------------------------------
  !Write a string in an element list
  CHARACTER(len=*),INTENT(IN):: string
  TYPE(strg_db),POINTER:: elm
  !Local variables
  INTEGER(kind=4):: i, lng

!  lng = LEN_TRIM(string)
  lng = MIN0(SIZE(elm%strng),LEN_TRIM(string))
  DO i=1,lng
    elm%strng(i) = string(i:i)
  END DO

  RETURN
  END SUBROUTINE wr_strg

!-----------------------------------------------------------------------
  SUBROUTINE rd_strg(elm,string)
!-----------------------------------------------------------------------
  !Read a string from an element list
  CHARACTER(len=*),INTENT(OUT):: string
  TYPE(strg_db),POINTER:: elm
  !Local variables
  INTEGER(kind=4):: i, lng

  string = ''

  lng = elm%len
  DO i=1,elm%len
    string(i:i) = elm%strng(i)
  END DO

  RETURN
  END SUBROUTINE rd_strg

!-----------------------------------------------------------------------
  SUBROUTINE search_strg(string,head,iel,elm)
!-----------------------------------------------------------------------
  !Read a string from an element list
  INTEGER(kind=4),INTENT(OUT):: iel
  CHARACTER(len=*),INTENT(IN):: string
  TYPE(strg_db),POINTER:: head, elm
  !Local variables
  INTEGER(kind=4):: i, lng
  LOGICAL:: found

  found = .FALSE.
  lng = LEN_TRIM(string)

  iel = 0
  elm => head
  Loop1 : DO
    IF (.NOT.ASSOCIATED(elm)) EXIT Loop1
    iel = iel + 1
    IF (lng <= elm%len) THEN
      DO i=1,lng
        IF (string(i:i) == elm%strng(i)) CYCLE
        elm => elm%next
        CYCLE Loop1
      END DO
      IF (lng < elm%len) THEN
        DO i=lng+1,elm%len
          IF (elm%strng(i) == ' ') CYCLE
          elm => elm%next
          CYCLE Loop1
        END DO
      END IF
      found = .TRUE.
      EXIT Loop1
    END IF
    elm => elm%next
  END DO Loop1

  IF (.NOT.found) THEN
    iel = 0
    NULLIFY(elm)
  END IF

  RETURN
  END SUBROUTINE search_strg

END MODULE vstrg_db

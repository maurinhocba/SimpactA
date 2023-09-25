SUBROUTINE delnod(iwrit, ndi)

  ! delete nodes

  USE ndinf_db

  IMPLICIT NONE
  INTEGER (kind=4),  INTENT (IN)    :: iwrit
  TYPE (list_ndata), INTENT (INOUT) :: ndi
END SUBROUTINE delnod

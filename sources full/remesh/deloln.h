SUBROUTINE deloln (iwrit, numpo, nodset, maxnn, ndi)

  !     delete old nodes (nodes used in an old mesh)

  USE ndinf_db

  IMPLICIT NONE
  INTEGER (kind=4), INTENT (IN) :: iwrit, numpo, nodset(:)
  INTEGER (kind=4), INTENT (OUT) :: maxnn
  TYPE (list_ndata), INTENT (INOUT) :: ndi

END SUBROUTINE deloln

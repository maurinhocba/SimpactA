FUNCTION vertex (coo_b, n1, n2, n3)


  IMPLICIT NONE
  LOGICAL :: vertex
  !        INPUT
  INTEGER (kind=4) :: n1, n2, n3   !number of nodes in contour line
  REAL    (kind=8) :: coo_b(:,:)   !(ndime,nelem) coordinates of contour nodes

END FUNCTION  vertex

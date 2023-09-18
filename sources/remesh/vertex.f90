FUNCTION vertex (coo_b, n1, n2, n3)

  ! check if angle < 165 degs - THIS IS MORE NATURAL IN GiD
  ! OLD: check if angle < 160 degs

  USE ctrl_db, ONLY: ndime

  IMPLICIT NONE
  LOGICAL :: vertex
  !        INPUT
  INTEGER (kind=4) :: n1, n2, n3   !number of nodes in contour line
  REAL    (kind=8) :: coo_b(:,:)   !(ndime,nelem) coordinates of contour nodes

  ! local
  REAL    (kind=8) :: t1(ndime),t2(ndime),xl,cosine
  !REAL    (kind=8), PARAMETER :: climit = -0.939692621d0 ! cos (160)
  REAL    (kind=8), PARAMETER :: climit = -0.965925826289068d0 ! cos (165)


  !=============================

  vertex = .FALSE.

  t1 = coo_b(1:ndime,n1) - coo_b(1:ndime,n2)   !vector n2-n1
  t2 = coo_b(1:ndime,n3) - coo_b(1:ndime,n2)   !vector n2-n3
  CALL vecuni (ndime, t1, xl)                  !normalization of t1
  CALL vecuni (ndime, t2, xl)                  !normalization of t2
  cosine = DOT_PRODUCT(t1, t2)              !cosine of the angle betw. t1, t2
  IF (cosine > climit) vertex = .TRUE.         !angle < 160

  RETURN

END FUNCTION  vertex

SUBROUTINE fixtem(iwrit,npoin,iftmp,label)

  !***  APPLY fixities

  USE c_input
  USE ift_db
  USE curv_db, ONLY: getcun
  IMPLICIT NONE
  INTEGER (kind=4), INTENT(IN) :: iwrit,npoin,label(:)
  INTEGER (kind=4), INTENT(IN OUT) :: iftmp(:,:)

END SUBROUTINE fixtem

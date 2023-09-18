 SUBROUTINE stgp06(ngaus,s,e,d)
 !*****************************************************************************
 !
 !*****evaluates total resultant stresses for shell element
 !     for linear elastic isotropic material
 !
 !****************************************************************************
 IMPLICIT NONE

 !                        routine parameters

 INTEGER (kind=4), INTENT(IN) :: ngaus

 REAL (kind=8), INTENT(IN) :: e(:,:),d(:)
 REAL (kind=8), INTENT(OUT) :: s(:,:)

 END SUBROUTINE stgp06

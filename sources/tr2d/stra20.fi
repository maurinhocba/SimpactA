 SUBROUTINE stra20(t,ntype,lb3,stran)

 !     Compute metric tensors, element CST (TLF)
 IMPLICIT NONE

 REAL (kind=8), INTENT(IN) :: t(2,2)  !derivatives of the element configuration
 INTEGER (kind=4), INTENT(IN) :: ntype  !problem type
 REAL (kind=8), INTENT(IN OUT)  :: lb3      !hoop stretch
 REAL (kind=8), INTENT(OUT) :: stran(4) !metric tensor


 ! metric tensor C = U^2
 stran(1) = DOT_PRODUCT(t(:,1),t(:,1))
 stran(2) = DOT_PRODUCT(t(:,2),t(:,2))
 stran(3) = DOT_PRODUCT(t(:,1),t(:,2))
 SELECT CASE (ntype)
 CASE (1)  !   plane stress
   stran(4) = 1d0/(stran(1)*stran(2) - stran(3)**2)  !isochoric
   lb3 = SQRT(stran(4))
 CASE (2)  !   plane strain
   stran(4) = 1d0
 CASE (3)  !   axilsymmetric
   stran(4) = lb3*lb3
 END SELECT

 RETURN
 END SUBROUTINE stra20

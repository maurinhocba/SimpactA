 SUBROUTINE bmem14(a,b,bmem,t,nnode)

 !*** membrane matrix  Constant Strain Triangle

 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: nnode
 REAL (kind=8), INTENT(IN) :: a(3),b(3),t(3,2)
 REAL (kind=8), INTENT(OUT) :: bmem(3,nnode,3)
 INTEGER (kind=4) :: j

 DO j=1,3          !node
   bmem(1:3,j,1) =  -b(j)*t(1:3,1)
   bmem(1:3,j,2) =   a(j)*t(1:3,2)
   bmem(1:3,j,3) =   a(j)*t(1:3,1)-b(j)*t(1:3,2)
 END DO

 RETURN
 END SUBROUTINE bmem14

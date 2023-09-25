 SUBROUTINE asstr6(ngamm,x,t,gamma)
 !*****************************************************************************
 !
 !*****natural coordinates shear deformation at assumed strain points
 !     for 4-node cuadrilateral shell element
 !*****************************************************************************
 IMPLICIT NONE

 !                        routine parameters

 INTEGER (kind=4), INTENT(IN) :: ngamm

 REAL (kind=8), INTENT(IN) ::  x(:,:),t(:,:)
 REAL (kind=8), INTENT(OUT) :: gamma(:)

 !                        local variables

 INTEGER (kind=4) j

 REAL (kind=8) dx(3,ngamm),tt(3,ngamm)

 tt(:,1) = t(:,1) + t(:,4)
 tt(:,2) = t(:,2) + t(:,1)
 tt(:,3) = t(:,3) + t(:,2)
 tt(:,4) = t(:,4) + t(:,3)

 dx(:,1) = x(:,4) - x(:,1)
 dx(:,2) = x(:,2) - x(:,1)
 dx(:,3) = x(:,3) - x(:,2)
 dx(:,4) = x(:,3) - x(:,4)

 DO j=1,ngamm
   gamma(j) = DOT_PRODUCT(tt(:,j),dx(:,j))/4d0
 END DO

 RETURN
 END SUBROUTINE asstr6

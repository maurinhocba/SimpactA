 SUBROUTINE psib06(psia,psib,jas)
 !*****************************************************************************
 !
 !*****natural coordinates shear deformation at assumed strain points
 !     for 4-node cuadrilateral shell element
 !*****************************************************************************
 IMPLICIT NONE

 !                        routine parameters


 REAL (kind=8), INTENT(IN) ::  psia(:,:),jas(:,:)
 REAL (kind=8), INTENT(OUT) :: psib(:)

 !                        local variables

 psib(1) = DOT_PRODUCT( jas(:,1),(psia(:,1) + psia(:,4)) ) /2d0
 psib(2) = DOT_PRODUCT( jas(:,2),(psia(:,2) + psia(:,1)) ) /2d0
 psib(3) = DOT_PRODUCT( jas(:,3),(psia(:,3) + psia(:,2)) ) /2d0
 psib(4) = DOT_PRODUCT( jas(:,4),(psia(:,4) + psia(:,3)) ) /2d0

 RETURN
 END SUBROUTINE psib06

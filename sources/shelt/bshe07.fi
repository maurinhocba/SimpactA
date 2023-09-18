 SUBROUTINE bshe07(x,t,lb,bs,bf,jas,nangl)
 !*****************************************************************************
 !
 !*****evaluates b matrix for shear at the assumed strain points
 !
 !****************************************************************************
 IMPLICIT NONE
 !                        routine parameters

 REAL (kind=8), INTENT(IN) :: x(:,:),t(:,:),lb(:,:,:)
 REAL (kind=8), INTENT(IN), POINTER :: jas(:,:),nangl(:,:)
 REAL (kind=8), INTENT(OUT) :: bs(:,:,:),bf(:,:,:)

 !                        local variables
 INTEGER (kind=4) :: i,j,k
 REAL    (kind=8) dx(3)

 ! contributions to vertex nodes (mid-surface displacements)
 bs(1:3,1,1) = 0d0    !node 1
 bs(1:3,2,1) =-t(:,1) !node 2
 bs(1:3,3,1) = t(:,1) !node 3
 bs(1:3,1,2) = t(:,2) !node 1
 bs(1:3,2,2) = 0d0    !node 2
 bs(1:3,3,2) =-t(:,2) !node 3
 bs(1:3,1,3) =-t(:,3) !node 1
 bs(1:3,2,3) = t(:,3) !node 2
 bs(1:3,3,3) = 0d0    !node 3
 ! contributions to mid-side nodes (director)
 dx = x(:,3)-x(:,2)  !side 2-3
 bs(4,1,1) = -DOT_PRODUCT(dx,lb(:,2,1))
 bs(5,1,1) =  DOT_PRODUCT(dx,lb(:,1,1))
 bs(4:5,2:3,1) = 0d0
 dx = x(:,1)-x(:,3)  !side 3-1
 bs(4:5,1,2) = 0d0
 bs(4,2,2) = -DOT_PRODUCT(dx,lb(:,2,2))
 bs(5,2,2) =  DOT_PRODUCT(dx,lb(:,1,2))
 bs(4:5,3,2) = 0d0
 dx = x(:,2)-x(:,1)  !side 1-2
 bs(4:5,1:2,3) = 0d0
 bs(4,3,3) = -DOT_PRODUCT(dx,lb(:,2,3))
 bs(5,3,3) =  DOT_PRODUCT(dx,lb(:,1,3))

 IF( ASSOCIATED(jas) )THEN
   DO i=1,3      !for each side
     DO j=2,3       !for each node on the side
        k=hh(j,i)      !node
        bf(1,k,i) = (nangl(1,k)*jas(1,i)-nangl(2,k)*jas(2,i))/2d0
        bf(2,k,i) = (nangl(1,k)*jas(2,i)+nangl(2,k)*jas(1,i))/2d0
     END DO
     bf(:,hh(1,i)-3,i) = 0d0
   END DO
 END IF

 RETURN
 END SUBROUTINE bshe07

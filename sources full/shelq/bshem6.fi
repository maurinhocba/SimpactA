 SUBROUTINE bshem6(x,t,lb,bs,ndofe,jas,nangl)
 !*****************************************************************************
 !
 !***** evaluates b matrix for shear at the assumed strain points
 !      for 4-node cuadrilateral shell element
 !****************************************************************************
 IMPLICIT NONE
 !                        routine parameters
 INTEGER(kind=4) :: ndofe
 REAL (kind=8), INTENT(IN) :: x(:,:),t(:,:),lb(:,:,:)
 REAL (kind=8), INTENT(IN), POINTER :: jas(:,:),nangl(:,:)
 REAL (kind=8), INTENT(OUT) :: bs(:,:,:)
 !                        local variables
 REAL (kind=8) :: tt(3),dx(3)

 !  point A
 tt = t(:,1) + t(:,4)
 dx = x(:,4) - x(:,1)
 bs(1:3,1,1) = -tt/4d0
 bs(1:3,4,1) =  tt/4d0
 bs(4,1,1)=-DOT_PRODUCT(dx,lb(:,2,1))/4d0
 bs(5,1,1)= DOT_PRODUCT(dx,lb(:,1,1))/4d0
 bs(4,4,1)=-DOT_PRODUCT(dx,lb(:,2,4))/4d0
 bs(5,4,1)= DOT_PRODUCT(dx,lb(:,1,4))/4d0
 !  point B
 tt = t(:,1) + t(:,2)
 dx = x(:,2) - x(:,1)
 bs(1:3,1,2) = -tt/4d0
 bs(1:3,2,2) =  tt/4d0
 bs(4,1,2)=-DOT_PRODUCT(dx,lb(:,2,1))/4d0
 bs(5,1,2)= DOT_PRODUCT(dx,lb(:,1,1))/4d0
 bs(4,2,2)=-DOT_PRODUCT(dx,lb(:,2,2))/4d0
 bs(5,2,2)= DOT_PRODUCT(dx,lb(:,1,2))/4d0
 !  point C
 tt = t(:,2) + t(:,3)
 dx = x(:,3) - x(:,2)
 bs(1:3,2,3) = -tt/4d0
 bs(1:3,3,3) =  tt/4d0
 bs(4,2,3)=-DOT_PRODUCT(dx,lb(:,2,2))/4d0
 bs(5,2,3)= DOT_PRODUCT(dx,lb(:,1,2))/4d0
 bs(4,3,3)=-DOT_PRODUCT(dx,lb(:,2,3))/4d0
 bs(5,3,3)= DOT_PRODUCT(dx,lb(:,1,3))/4d0
 ! point D
 tt = t(:,4) + t(:,3)
 dx = x(:,3) - x(:,4)
 bs(1:3,4,4) = -tt/4d0
 bs(1:3,3,4) =  tt/4d0
 bs(4,4,4)=-DOT_PRODUCT(dx,lb(:,2,4))/4d0
 bs(5,4,4)= DOT_PRODUCT(dx,lb(:,1,4))/4d0
 bs(4,3,4)=-DOT_PRODUCT(dx,lb(:,2,3))/4d0
 bs(5,3,4)= DOT_PRODUCT(dx,lb(:,1,3))/4d0
 IF( ndofe == 7 )THEN
   ! node A (-1, 0)  1 -> 4
   bs(6,1,1) = (nangl(1,1)*jas(1,1)-nangl(2,1)*jas(2,1))*0.5d0
   bs(7,1,1) = (nangl(1,1)*jas(2,1)+nangl(2,1)*jas(1,1))*0.5d0
   bs(6,4,1) = (nangl(1,4)*jas(1,1)-nangl(2,4)*jas(2,1))*0.5d0
   bs(7,4,1) = (nangl(1,4)*jas(2,1)+nangl(2,4)*jas(1,1))*0.5d0
   ! node B ( 0,-1)  1 -> 2
   bs(6,1,2) = (nangl(1,1)*jas(1,2)-nangl(2,1)*jas(2,2))*0.5d0
   bs(7,1,2) = (nangl(1,1)*jas(2,2)+nangl(2,1)*jas(1,2))*0.5d0
   bs(6,2,2) = (nangl(1,2)*jas(1,2)-nangl(2,2)*jas(2,2))*0.5d0
   bs(7,2,2) = (nangl(1,2)*jas(2,2)+nangl(2,2)*jas(1,2))*0.5d0
   ! node C ( 1, 0) 2 - > 3
   bs(6,2,3) = (nangl(1,2)*jas(1,3)-nangl(2,2)*jas(2,3))*0.5d0
   bs(7,2,3) = (nangl(1,2)*jas(2,3)+nangl(2,2)*jas(1,3))*0.5d0
   bs(6,3,3) = (nangl(1,3)*jas(1,3)-nangl(2,3)*jas(2,3))*0.5d0
   bs(7,3,3) = (nangl(1,3)*jas(2,3)+nangl(2,3)*jas(1,3))*0.5d0
   ! node D ( 0, 1) 4 -> 3
   bs(6,4,4) = (nangl(1,4)*jas(1,4)-nangl(2,4)*jas(2,4))*0.5d0
   bs(7,4,4) = (nangl(1,4)*jas(2,4)+nangl(2,4)*jas(1,4))*0.5d0
   bs(6,3,4) = (nangl(1,3)*jas(1,4)-nangl(2,3)*jas(2,4))*0.5d0
   bs(7,3,4) = (nangl(1,3)*jas(2,4)+nangl(2,3)*jas(1,4))*0.5d0
 END IF
 RETURN
 END SUBROUTINE bshem6

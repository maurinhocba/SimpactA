 SUBROUTINE shape9(weigh,posgp,shape,deriv,nnode,ngaus)
 !***********************************************************************
 !
 !     gauss points, weights, shape functions and derivatives for
 !     one-dimensinal elements (beams and 2-d shells)
 !***********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nnode,ngaus
 REAL (kind=8), INTENT(OUT) :: weigh(:),posgp(:),shape(:,:),deriv(:,:)

 !     gauss points weight and possitions

 SELECT CASE (ngaus)
 CASE(1)
   weigh(1) =  2d0
   posgp(1) =  0d0
   SELECT CASE (nnode)
   CASE(2) !   two nodes, one gauss point
     shape(1,1) =  0.5d0
     shape(2,1) =  0.5d0
     deriv(1,1) = -0.5d0
     deriv(2,1) =  0.5d0
   END SELECT
 CASE(2)
   weigh(1:2) =  1d0
   posgp(1) = -1d0/SQRT(3d0)
   posgp(2) =  1d0/SQRT(3d0)
   SELECT CASE (nnode)
   CASE (2) !     two nodes, two gauss points
     shape(1,1) =  0.5d0*(1d0-posgp(1))
     shape(2,1) =  0.5d0*(1d0+posgp(1))
     shape(1,2) =  0.5d0*(1d0-posgp(2))
     shape(2,2) =  0.5d0*(1d0+posgp(2))
     deriv(1,1) = -0.5d0
     deriv(2,1) =  0.5d0
     deriv(1,2) = -0.5d0
     deriv(2,2) =  0.5d0
   CASE (3)  !     three nodes, two gauss points
     shape(1,1) = -0.5d0*(1d0-posgp(1))*posgp(1)
     shape(2,1) =  1d0-posgp(1)**2
     shape(3,1) =  0.5d0*(1d0+posgp(1))*posgp(1)
     shape(1,2) = -0.5d0*(1d0-posgp(2))*posgp(2)
     shape(2,2) =  1d0-posgp(2)**2
     shape(3,2) =  0.5d0*(1d0+posgp(2))*posgp(2)
     deriv(1,1) =  posgp(1)-0.5d0
     deriv(2,1) = -posgp(1)*2d0
     deriv(3,1) =  posgp(1)+0.5d0
     deriv(1,2) =  posgp(2)-0.5d0
     deriv(2,2) = -posgp(2)*2d0
     deriv(3,2) =  posgp(2)+0.5d0
   END SELECT

 CASE(3)     !Three Gauss points
   weigh(1) =  5d0/9d0
   weigh(2) =  8d0/9d0
   weigh(3) =  5d0/9d0
   posgp(1) = -0.774596669241483d0
   posgp(2) =  0d0
   posgp(3) =  0.774596669241483d0
   SELECT CASE (nnode)
   CASE (3)  !     three nodes, three gauss points
     shape(1,1) = -0.5d0*(1d0-posgp(1))*posgp(1)
     shape(2,1) =  1d0-posgp(1)**2
     shape(3,1) =  0.5d0*(1d0+posgp(1))*posgp(1)
     shape(1,2) = -0.5d0*(1d0-posgp(2))*posgp(2)
     shape(2,2) =  1d0-posgp(2)**2
     shape(3,2) =  0.5d0*(1d0+posgp(2))*posgp(2)
     shape(1,3) = -0.5d0*(1d0-posgp(3))*posgp(3)
     shape(2,3) =  1d0-posgp(3)**2
     shape(3,3) =  0.5d0*(1d0+posgp(3))*posgp(3)
     deriv(1,1) =  posgp(1)-0.5d0
     deriv(2,1) = -posgp(1)*2d0
     deriv(3,1) =  posgp(1)+0.5d0
     deriv(1,2) =  posgp(2)-0.5d0
     deriv(2,2) = -posgp(2)*2d0
     deriv(3,2) =  posgp(2)+0.5d0
     deriv(1,3) =  posgp(3)-0.5d0
     deriv(2,3) = -posgp(3)*2d0
     deriv(3,3) =  posgp(3)+0.5d0
   END SELECT
 END SELECT
 RETURN
 END SUBROUTINE shape9

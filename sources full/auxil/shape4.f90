 SUBROUTINE shape4 (nnode, shape, deriv ,xita, eta, zeta )
 !********************************************************************
 !
 !***  calculates shape functions and their derivatives for 3d
 !     prismatic elements
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER(Kind=4) :: nnode
 REAL (kind=8)    deriv(nnode,3),shape(nnode),xita,eta,zeta

 REAL (kind=8) l1,l2,zzet,bf

   l1 = (1d0-zeta)/2d0  !linear functions L1
   l2 = (1d0+zeta)/2d0  !                 L2
   zzet = 1d0-xita-eta  !third area coordinate

   IF( nnode == 6 )THEN
     !*** shape functions

     shape(1) = zzet * l1
     shape(2) = xita * l1
     shape(3) = eta  * l1
     shape(4) = zzet * l2
     shape(5) = xita * l2
     shape(6) = eta  * l2

     !*** and derivatives

     deriv(1,1) = -l1
     deriv(2,1) =  l1
     deriv(3,1) =  0d0
     deriv(4,1) = -l2
     deriv(5,1) =  l2
     deriv(6,1) =  0d0

     deriv(1,2) = -l1
     deriv(2,2) =  0d0
     deriv(3,2) =  l1
     deriv(4,2) = -l2
     deriv(5,2) =  0d0
     deriv(6,2) =  l2

     deriv(1,3) = -zzet/2d0
     deriv(2,3) = -xita/2d0
     deriv(3,3) = -eta /2d0
     deriv(4,3) =  zzet/2d0
     deriv(5,3) =  xita/2d0
     deriv(6,3) =  eta /2d0

   ELSE IF( nnode == 15 )THEN
     bf = 1d0 - zeta**2                      !bubble function
     shape(1) = zzet*(2d0*zzet-1d0)*l1 - zzet*bf/2d0
     shape(2) = xita*(2d0*xita-1d0)*l1 - xita*bf/2d0
     shape(3) = eta *(2d0*eta -1d0)*l1 - eta *bf/2d0
     shape(4) = zzet*(2d0*zzet-1d0)*l2 - zzet*bf/2d0
     shape(5) = xita*(2d0*xita-1d0)*l2 - xita*bf/2d0
     shape(6) = eta *(2d0*eta -1d0)*l2 - eta *bf/2d0
     shape(7) = 4d0*zzet*xita*l1
     shape(8) = 4d0*xita*eta *l1
     shape(9) = 4d0*eta *zzet*l1
     shape(10) = zzet*bf
     shape(11) = xita*bf
     shape(12) = eta *bf
     shape(13) = 4d0*zzet*xita*l2
     shape(14) = 4d0*xita*eta *l2
     shape(15) = 4d0*eta *zzet*l2

     deriv(1,1) = (-4d0*zzet+1d0)*l1 + bf/2d0
     deriv(2,1) =  (4d0*xita-1d0)*l1 - bf/2d0
     deriv(3,1) =  0d0
     deriv(4,1) = (-4d0*zzet+1d0)*l2 + bf/2d0
     deriv(5,1) =  (4d0*xita-1d0)*l2 - bf/2d0
     deriv(6,1) =  0d0
     deriv(7,1) =  4d0*(zzet-xita)*l1
     deriv(8,1) =  4d0*eta *       l1
     deriv(9,1) = -4d0*eta *       l1
     deriv(10,1) = -bf
     deriv(11,1) =  bf
     deriv(12,1) =  0d0
     deriv(13,1) =  4d0*(zzet-xita)*l2
     deriv(14,1) =  4d0*eta *       l2
     deriv(15,1) = -4d0*eta *       l2

     deriv(1,2) = (-4d0*zzet+1d0)*l1 + bf/2d0
     deriv(2,2) =  0d0
     deriv(3,2) =  (4d0*eta -1d0)*l1 - bf/2d0
     deriv(4,2) = (-4d0*zzet+1d0)*l2 + bf/2d0
     deriv(5,2) =  0d0
     deriv(6,2) =  (4d0*eta -1d0)*l2 - bf/2d0
     deriv(7,2) = -4d0*xita*       l1
     deriv(8,2) =  4d0*xita*       l1
     deriv(9,2) =  4d0*(-eta+zzet)*l1
     deriv(10,2) = -bf
     deriv(11,2) =  0d0
     deriv(12,2) =  bf
     deriv(13,2) = -4d0*xita*      l2
     deriv(14,2) =  4d0*xita*      l2
     deriv(15,2) = 4d0*(-eta+zzet)*l2

     deriv(1,3) = -zzet*(2d0*zzet-1d0)/2d0 + zzet*zeta
     deriv(2,3) = -xita*(2d0*xita-1d0)/2d0 + xita*zeta
     deriv(3,3) = -eta *(2d0*eta -1d0)/2d0 + eta *zeta
     deriv(4,3) =  zzet*(2d0*zzet-1d0)/2d0 + zzet*zeta
     deriv(5,3) =  xita*(2d0*xita-1d0)/2d0 + xita*zeta
     deriv(6,3) =  eta *(2d0*eta -1d0)/2d0 + eta *zeta
     deriv(7,3) = -2d0*zzet*xita
     deriv(8,3) = -2d0*xita*eta
     deriv(9,3) = -2d0*eta *zzet
     deriv(10,3) = -zzet*2d0*zeta
     deriv(11,3) = -xita*2d0*zeta
     deriv(12,3) = -eta *2d0*zeta
     deriv(13,3) = 2d0*zzet*xita
     deriv(14,3) = 2d0*xita*eta
     deriv(15,3) = 2d0*eta *zzet
   END IF

 RETURN
 END SUBROUTINE shape4

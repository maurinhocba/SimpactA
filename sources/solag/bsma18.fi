 SUBROUTINE bsma18 (gas,tas,nfdas,bs)
 ! compute Strain-displacement matrix for transverse shear
 ! at assumed strain points for solid element 18 (8 nodes)
 IMPLICIT NONE
 REAL (kind=8) :: tas(3,4),gas(3,4),nfdas(8,4)
 REAL (kind=8) :: bs(24,4)

 INTEGER (kind=4) :: i,j,k,l

 DO i=1,4  !for each assumed strain point
   k = 0
   DO j=1,8  !for each node
     DO l=1,3
       k = k+1
       bs(k,i) = gas(l,i)*nfdas(j,i)
     END DO
   END DO
 END DO
 tas = tas/4d0
 ! point A (1)     gas(:,1) = (-x(:,1)-x(:,5)+x(:,4)+x(:,8))/4d0
 bs( 1: 3,1) = bs( 1: 3,1) - tas(:,1)
 bs(10:12,1) = bs(10:12,1) + tas(:,1)
 bs(13:15,1) = bs(13:15,1) - tas(:,1)
 bs(22:24,1) = bs(22:24,1) + tas(:,1)
 ! point B (2)     gas(:,2) = (-x(:,1)-x(:,5)+x(:,2)+x(:,6))/4d0
 bs( 1: 3,2) = bs( 1: 3,2) - tas(:,2)
 bs( 4: 6,2) = bs( 4: 6,2) + tas(:,2)
 bs(13:15,2) = bs(13:15,2) - tas(:,2)
 bs(16:18,2) = bs(16:18,2) + tas(:,2)
 ! point C (3)     gas(:,3) = (-x(:,2)-x(:,6)+x(:,3)+x(:,7))/4d0
 bs( 4: 6,3) = bs( 4: 6,3) - tas(:,3)
 bs( 7: 9,3) = bs( 7: 9,3) + tas(:,3)
 bs(16:18,3) = bs(16:18,3) - tas(:,3)
 bs(19:21,3) = bs(19:21,3) + tas(:,3)
 ! poinr D (4)     gas(:,4) = (-x(:,4)-x(:,8)+x(:,3)+x(:,7))/4d0
 bs( 7: 9,4) = bs( 7: 9,4) + tas(:,4)
 bs(10:12,4) = bs(10:12,4) - tas(:,4)
 bs(19:21,4) = bs(19:21,4) + tas(:,4)
 bs(22:24,4) = bs(22:24,4) - tas(:,4)
 RETURN
 END SUBROUTINE bsma18

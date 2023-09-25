 SUBROUTINE bsma16 (gas,tas,nfdas,bs,sh)
 ! compute Strain-displacement matrix for transverse shear
 ! at assumed strain points for PRISM element 16 (6 nodes)
 IMPLICIT NONE
 REAL (kind=8) :: tas(3,3),gas(3,3),nfdas(6,3),sh(2)
 REAL (kind=8) :: bs(3,6,3)

 INTEGER (kind=4) :: i,j,l

 DO i=1,3  !for each assumed strain point
   DO j=1,6  !for each node
     DO l=1,3   !for each nodal DOF
       bs(l,j,i) = gas(l,i)*nfdas(j,i)
     END DO
   END DO
 END DO
 ! point (1)
 bs(:,2,1) = bs(:,2,1) - tas(:,1)*sh(1)
 bs(:,3,1) = bs(:,3,1) + tas(:,1)*sh(1)
 bs(:,5,1) = bs(:,5,1) - tas(:,1)*sh(2)
 bs(:,6,1) = bs(:,6,1) + tas(:,1)*sh(2)
 ! point (2)
 bs(:,1,2) = bs(:,1,2) + tas(:,2)*sh(1)
 bs(:,3,2) = bs(:,3,2) - tas(:,2)*sh(1)
 bs(:,4,2) = bs(:,4,2) + tas(:,2)*sh(2)
 bs(:,6,2) = bs(:,6,2) - tas(:,2)*sh(2)
 ! point (3)
 bs(:,1,3) = bs(:,1,3) - tas(:,3)*sh(1)
 bs(:,2,3) = bs(:,2,3) + tas(:,3)*sh(1)
 bs(:,4,3) = bs(:,4,3) - tas(:,3)*sh(2)
 bs(:,5,3) = bs(:,5,3) + tas(:,3)*sh(2)
 RETURN
 END SUBROUTINE bsma16

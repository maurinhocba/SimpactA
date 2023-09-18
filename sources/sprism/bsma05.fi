 SUBROUTINE bsma05 (gas,tas,nfdas,bs,ji,ii)
 ! compute Strain-displacement matrix for transverse shear
 ! at assumed strain points for SPRISM  solid-shell element 05 (6 nodes)
 IMPLICIT NONE
 REAL (kind=8) :: tas(3,3),gas(3,3),nfdas(6,3),ji(2,2)
 REAL (kind=8) :: bs(3,6,2)
 INTEGER (kind=4) :: ii

 INTEGER (kind=4) :: i,j,l
 REAL (kind=8) :: baux(3,2),bm(3,6,3)

 DO i=1,3  !for each assumed strain point
   DO j=1,6  !for each node
     DO l=1,3   !for each nodal DOF
       bm(l,j,i) = gas(l,i)*nfdas(j,i)
     END DO
   END DO
 END DO
 IF( ii == 1 ) THEN
   ! point (1)
   bm(:,2,1) = bm(:,2,1) - tas(:,1)
   bm(:,3,1) = bm(:,3,1) + tas(:,1)
   ! point (2)
   bm(:,1,2) = bm(:,1,2) + tas(:,2)
   bm(:,3,2) = bm(:,3,2) - tas(:,2)
   ! point (3)
   bm(:,1,3) = bm(:,1,3) - tas(:,3)
   bm(:,2,3) = bm(:,2,3) + tas(:,3)
 ELSE
   ! point (1)
   bm(:,5,1) = bm(:,5,1) - tas(:,1)
   bm(:,6,1) = bm(:,6,1) + tas(:,1)
   ! point (2)
   bm(:,4,2) = bm(:,4,2) + tas(:,2)
   bm(:,6,2) = bm(:,6,2) - tas(:,2)
   ! point (3)
   bm(:,4,3) = bm(:,4,3) - tas(:,3)
   bm(:,5,3) = bm(:,5,3) + tas(:,3)
 END IF
 DO j=1,6
   baux(:,1) =  (-bm(:,j,1) - bm(:,j,2)  + 2d0*bm(:,j,3))/3d0
   baux(:,2) =  ( bm(:,j,1) + bm(:,j,3)  - 2d0*bm(:,j,2))/3d0
   bs(1:3,j,1) = ji(1,1)*baux(:,1)+ ji(2,1)*baux(:,2)
   bs(1:3,j,2) = ji(1,2)*baux(:,1)+ ji(2,2)*baux(:,2)
 END DO
 RETURN
 END SUBROUTINE bsma05

 SUBROUTINE nodxy7(nelem,heade,coord,eule0,euler,ngaus)
 !*********************************************************************
 !
 !     interpolates coordinates for 6-node triangles
 !
 !*********************************************************************
 USE npo_db, ONLY : coora
 IMPLICIT NONE
 INTEGER (kind=4) nelem,ngaus
 REAL    (kind=8) coord(:,:),eule0(:,:),euler(:,:)
 TYPE(ele07), POINTER :: heade

 INTEGER (kind=4) i,k,l,n,m
 TYPE(ele07), POINTER :: e
 INTEGER(kind=4), PARAMETER :: kk(2,3) = (/ 2,3, 3,1, 1,2 /)

 e => heade
 DO n = 1,nelem
   IF( ngaus /= 1 )THEN
     DO k = 4,6
       i = e%lnods(k)
       IF(ALL(coord(:,i) == 0d0) ) THEN
         l = e%lnods(k-3)
         m = e%lnods(MOD(k,3)+1)
         coord(:,i) = (coord(:,l)+coord(:,m))/2d0
         coora(:,i) = coord(:,i)
       END IF
       IF(ABS(eule0(1,i))+ABS(eule0(2,i))+ABS(eule0(3,i)) == 0d0) THEN
         l = e%lnods(k-3)
         m = e%lnods(MOD(k,3)+1)
         eule0(1,i) = (eule0(1,l)+eule0(1,m))/2d0
         eule0(2,i) = (eule0(2,l)+eule0(2,m))/2d0
         eule0(3,i) = (eule0(3,l)+eule0(3,m))/2d0
         euler(1:3,i) = eule0(1:3,i)
       END IF
     END DO
   ELSE
     DO k = 4,6
       i = e%lnods(k)
       l = e%lnods(kk(1,k-3))
       m = e%lnods(kk(2,k-3))
       coord(:,i) = (coord(:,l)+coord(:,m))/2d0
     END DO
   END IF
   e => e%next
 END DO
 RETURN
 END SUBROUTINE nodxy7

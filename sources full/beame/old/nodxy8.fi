 SUBROUTINE nodxy8(ndime,nnode,lnods,coord,eule0,euler)
 !*********************************************************************
 !
 !     interpolates coordinates for 3-node beam elements
 !
 !*********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) ndime,nnode,lnods(nnode)
 REAL    (kind=8) coord(ndime,*),eule0(3,*),euler(9,*)

 INTEGER (kind=4) i,j,l,m

   m = lnods(nnode)
   IF(m /= 0) THEN
     i = lnods(2)
     l = lnods(1)
     IF(ALL(coord(1:3,i) == 0))THEN
       DO j=1,3
         coord(j,i)=(coord(j,l)+coord(j,m))/2d0
       END DO
     END IF
     IF(ABS(eule0(1,i))+ABS(eule0(2,i))+ABS(eule0(3,i)) == 0.0)THEN
       eule0(1,i) = (eule0(1,l)+eule0(1,m))/2d0
       eule0(2,i) = (eule0(2,l)+eule0(2,m))/2d0
       eule0(3,i) = (eule0(3,l)+eule0(3,m))/2d0
       euler(1:3,i) = eule0(1:3,i)
     END IF
   END IF
 RETURN
 END SUBROUTINE nodxy8

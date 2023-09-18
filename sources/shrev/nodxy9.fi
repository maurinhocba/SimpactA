 SUBROUTINE nodxy9(nelem,heade,coord,eule0,euler)
 !*********************************************************************
 !
 !     interpolates coordinates for 3-node beam/shell elements
 !
 !*********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nelem
 REAL (kind=8), INTENT(IN OUT) :: coord(:,:),eule0(:,:),euler(:,:)
 TYPE(ele09), POINTER :: heade

 !local variables
 INTEGER (kind=4) i,l,n,m
 TYPE(ele09), POINTER :: e

 e => heade
 DO n=1,nelem
   l = e%lnods(1)
   i = e%lnods(2)
   m = e%lnods(3)
   IF(ALL(coord(:,i) == 0d0)) coord(:,i) = (coord(:,l)+coord(:,m))/2d0
   IF(ALL(eule0(:,i) == 0d0))THEN
     eule0(:,i) = (eule0(:,l)+eule0(:,m))/2d0
     euler(:,i) = eule0(:,i)
   END IF
   e => e%next
 END DO
 RETURN
 END SUBROUTINE nodxy9

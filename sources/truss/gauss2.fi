 SUBROUTINE gauss2(ndime,nelem,heade,coord,flag,istop)
 !********************************************************************
 !
 !*** evaluates gauss points constants for truss element
 !
 !********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) :: ndime,nelem,istop
 REAL    (kind=8) :: coord(:,:)
 TYPE (ele02), POINTER :: heade
 LOGICAL flag

 INTEGER (kind=4) i,ielem
 REAL    (kind=8) x(ndime,nnode),l0
 TYPE (ele02), POINTER :: e

 e => heade
 DO ielem = 1,nelem
   !       gathers the initial coordinates and local coordinate transf.
   x  = coord(1:ndime,e%lnods)
   l0 = 0d0
   DO i = 1,ndime
     l0  = l0 + (x(i,2)-x(i,1))**2
   END DO
   l0 = SQRT(l0)
   IF(flag .AND. e%gausv(1) == 0d0) THEN
     e%gausv(1) = l0
     IF(l0 <= 0) THEN
       WRITE(lures,1000,ERR=9999) e%numel,e%gausv(1)
       istop = 1
       RETURN
     END IF
   END IF
   e => e%next
 END DO

  RETURN
 1000 FORMAT(' null length for truss element (1) ',i5,/                 &
&       ' jacobian = ', e12.4 )
 9999 CALL runen2('')
 END SUBROUTINE gauss2

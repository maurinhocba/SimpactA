 SUBROUTINE acvdf2(ndime,nelem,ifpre,heade)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE
 INTEGER(kind=4),INTENT(IN) :: ndime,nelem
 INTEGER(kind=4),INTENT(OUT) :: ifpre(:,:)
 TYPE(ele02), POINTER :: heade

 INTEGER (kind=4) ielem
 TYPE(ele02), POINTER :: e

 !*** activates degrees of freedom

 e => heade
 DO ielem=1,nelem
   ifpre(1:ndime,e%lnods) = 0
   e => e%next
 END DO

 RETURN
 END SUBROUTINE acvdf2

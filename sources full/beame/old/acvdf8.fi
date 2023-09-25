 SUBROUTINE acvdf8(nelem,nnode,ifpre,heade)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 USE esets_db, ONLY : gnods, gelem !INTENT(OUT) gnods(maxnn,melen)
 USE ctrl_db, ONLY : lumped        !INTENT(IN)
 USE npo_db, ONLY : naeul
 IMPLICIT NONE
 INTEGER (kind=4) :: nelem,nnode,ifpre(:,:)
 TYPE (ele08), POINTER :: heade

 INTEGER (kind=4), PARAMETER :: ndofe = 6
 INTEGER (kind=4) ielem
 TYPE (ele08), POINTER :: e

 !*** activates degrees of freedom

 e => heade
 DO ielem=1,nelem
   ifpre(1:ndofe,e%lnods(1:nnode)) = 0
   naeul(e%lnods(1:nnode)) = .TRUE.
   IF( .NOT.lumped )THEN
     gelem = gelem + 1                !update total number of elements
     gnods(1:nnode,gelem) = e%lnods(1:nnode)  !pass to global array
   END IF
   e => e%next
 END DO

 RETURN
 END SUBROUTINE acvdf8

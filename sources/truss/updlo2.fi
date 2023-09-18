SUBROUTINE updlo2(nelem,heade,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE
 INTEGER(kind=4),INTENT(IN) :: nelem,oldlb(:)
 TYPE(ele02), POINTER :: heade

 INTEGER (kind=4) ielem,i,lab,chnode
 TYPE(ele02), POINTER :: e

 !old local nodes are changed into new ones

 e => heade
 DO ielem=1,nelem
   DO i=1,nnode
     lab  = oldlb(e%lnods(i))
     e%lnods(i) = chnode(lab)
   END DO
   e => e%next
 END DO

 RETURN
END SUBROUTINE updlo2

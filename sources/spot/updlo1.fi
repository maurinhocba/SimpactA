 SUBROUTINE updlo1(nnode,nelem,heade,oldlb)
 !******************************************************************
 !
 !*** updates internal node numbering
 !
 !******************************************************************
 IMPLICIT NONE
 INTEGER(kind=4),INTENT(IN) :: nnode,nelem,oldlb(:)
 TYPE(ele01), POINTER :: heade

 INTEGER (kind=4) ielem,i,lab,chnode
 TYPE(ele01), POINTER :: e

 !old local nodes are changed into new ones

 e => heade
 DO ielem=1,nelem
   DO i=1,nnode
     lab  = oldlb(e%lnods(i))
     IF( lab > 0 )e%lnods(i) = chnode(lab)
   END DO
   e => e%next
 END DO

 RETURN
 END SUBROUTINE updlo1

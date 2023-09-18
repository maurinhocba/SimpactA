 SUBROUTINE slnd19(heade,lnods,nelem)
 !***********************************************************************
 !*** Extracts conectivities into the nodal set
 !***********************************************************************
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4), POINTER :: lnods(:,:)   !Conectivities of the set
   INTEGER(kind=4),INTENT(IN):: nelem       !number of elements
   TYPE(ele19),POINTER:: heade      !pointer to first element

   !--- Local Variables
   INTEGER(kind=4):: ielem,nn
   TYPE(ele19),POINTER:: elm


   IF( ASSOCIATED(lnods) )THEN
     nn    = SIZE(lnods,1)
     ielem = SIZE(lnods,2)
     IF( ielem /= nelem .OR. nn /= nnode)THEN
       DEALLOCATE(lnods)
       ALLOCATE(lnods(nnode,nelem))
     END IF
   ELSE
     ALLOCATE(lnods(nnode,nelem))
   END IF
   ielem = 0                             !Initiates counter
   elm => heade                          !Point to first element for each element in the set
   DO
     IF (.NOT.ASSOCIATED(elm)) EXIT      !Exit from the loop
     ielem = ielem + 1

     lnods(1:nnode,ielem) = elm%lnods(1:nnode)   !transfer data

     elm => elm%next                     !point to next element
   END DO

 RETURN
 END SUBROUTINE slnd19

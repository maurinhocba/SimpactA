 SUBROUTINE acvd11(ifpre,heade,lside,nelem,nbn,nhead)
 !******************************************************************
 !
 !*** updates internal node numbering and activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 LOGICAL, INTENT(IN OUT) :: lside !flag TRUE if topology already computed
 INTEGER (kind=4), INTENT(IN) :: nelem         !number of elements
 INTEGER (kind=4), INTENT(IN OUT) :: nbn, &    !number of branching points
                                     ifpre(:,:)
 TYPE (ele11), POINTER :: heade   !pointer to first element in the set
 TYPE (nodei), POINTER :: nhead   !pointer to first branching node

 INTEGER (kind=4) i
 TYPE (ele11), POINTER  :: e                     !element pointers


 !***  create topological array  lnods(4:6,nelem)

 IF(.NOT.lside) CALL toar11(heade,lside,nbn,nhead,nelem)

 !*** activates degrees of freedom

 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
   DO i=2,3
     IF( ifpre(1,e%lnods(i)) == 1 ) & !to avoid to overwrite RIGID constraints
       ifpre(1:2,e%lnods(i)) = 0      !release DOFs

   END DO
   e => e%next                      !point to next element
 END DO

 RETURN
 END SUBROUTINE acvd11

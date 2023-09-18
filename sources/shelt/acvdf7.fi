 SUBROUTINE acvdf7(ifpre,heade,ngaus,nnode,lside,nelem,zigzag)
 !******************************************************************
 !
 !*** activates degrees of freedom
 !
 !******************************************************************
 USE npo_db, ONLY : naeul
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: ifpre(:,:)    !(OUT) DOFs
 INTEGER (kind=4):: ngaus,nnode,nelem    !(IN) element type
 LOGICAL :: lside
 TYPE (ele07), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT(IN) :: zigzag

 TYPE (section), POINTER :: sec      !auxiliar pointer
 TYPE (ele07), POINTER  :: e                     !element pointers
 INTEGER (kind=4) :: osec,isec,secty

!***  create topological array  lnods(7:9,nelem)

 IF(nnode == 9) THEN
   IF(.NOT.lside) CALL toar07(heade,lside,nelem)
 END IF

!*** activates degrees of freedom

   e => heade                         !point to first element
   osec = -1
   DO
     IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
     isec = e%matno
     IF( isec /= osec )THEN
       sec => psecs(isec)%p
       secty = sec%secty    !section constitutive type
       osec = isec
     END IF
     IF( ngaus > 1 )THEN
       ifpre(1:5,e%lnods) = 0           !release DOFs
       naeul(e%lnods) = .TRUE.          !release DOFs
       IF( secty == 15 .AND. zigzag ) ifpre(7:8,e%lnods) = 0
     ELSE
       ifpre(1:3,e%lnods(1:3)) = 0      !release Translational DOFs
       IF( secty == 15 .AND. zigzag ) ifpre(7:8,e%lnods(1:3)) = 0
       ifpre(4:5,e%lnods(4:6)) = 0      !release Rotational DOFs
       naeul(e%lnods(4:6)) = .TRUE.     !release local systems
     END IF
     e => e%next                      !point to next element
   END DO

 RETURN
 END SUBROUTINE acvdf7

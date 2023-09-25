 SUBROUTINE acvdf6(ifpre,heade,zigzag)
 !******************************************************************
 !
 !*** activates degrees of freedom
 !
 !******************************************************************
 USE npo_db, ONLY : naeul
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4):: ifpre(:,:)    !(OUT) DOFs
 TYPE (ele06), POINTER :: heade   !pointer to first element in the set
 LOGICAL, INTENT(IN) :: zigzag

 TYPE (section), POINTER :: sec      !auxiliar pointer
 TYPE (ele06), POINTER  :: e                     !element pointers
 INTEGER (kind=4) :: i,osec,isec,secty

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
     DO i=1,nnode
       IF(ifpre(1,e%lnods(i)) == 1)THEN
         ifpre(1:5,e%lnods(i)) = 0       !release DOFs
         IF( secty == 15 .AND. zigzag ) ifpre(7:8,e%lnods(i)) = 0
       END IF  
     END DO
     naeul(e%lnods) = .TRUE.          !release DOFs
     e => e%next                      !point to next element
   END DO

 RETURN
 END SUBROUTINE acvdf6

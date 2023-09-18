 SUBROUTINE acvdf9(ifpre,heade,naeul)
 !******************************************************************
 !
 !*** activates degrees of freedom
 !
 !******************************************************************
 IMPLICIT NONE

 ! dummy arguments
 INTEGER (kind=4), INTENT(OUT) :: ifpre(:,:) !DOFs
 LOGICAL, INTENT(OUT) :: naeul(:)            !exists
 TYPE (ele09), POINTER :: heade   !pointer to first element in the set

 TYPE (ele09), POINTER  :: e                     !element pointers

 !*** activates degrees of freedom

   e => heade                         !point to first element
   DO
     IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, Exit loop
     ifpre(1:ndofe,e%lnods) = 0       !release DOFs
     naeul(e%lnods) = .TRUE.          !release DOFs
     e => e%next                      !point to next element
   END DO

 RETURN
 END SUBROUTINE acvdf9
